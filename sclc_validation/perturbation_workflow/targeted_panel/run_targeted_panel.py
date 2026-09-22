#!/usr/bin/env python3
"""Targeted, one-gene-at-a-time delete + overexpress screen for the 50-gene
panel (21 pre-registered signature genes + 29 top drivers from the prior
LUAD/LUSC/normal all-gene screen, contamination markers excluded).

Reuses the train-reference/held-out-test datasets and training-donor state
embeddings already computed for the all-gene sweep (same model, same split
-- nothing about switching to a gene-targeted approach invalidates them).

For "delete", InSilicoPerturber filters each source's cells down to only
those that actually detect the target gene before perturbing (cheap). For
"overexpress" that filter does not apply (a gene can be induced from
undetected), so every held-out cell in the source is processed for every
gene -- this is the more expensive of the two per gene, by design.

OUTPUT LOCATION (changed 2026-09-17, bf16 precision-replication task): by
default this script now writes under
sclc_validation/bf16_bench/runs/<run-tag>/targeted_panel/ instead of any
directory under this script's own location -- <run-tag> defaults to
--dtype ("fp32" or "bf16") and can be set explicitly with --run-tag so
parallel arms (e.g. an fp32 baseline vs its fp32 noise-floor repeat) don't
clobber each other. A canonical full-panel run with no bf16-specific flags
will land at .../bf16_bench/runs/fp32/targeted_panel/, NOT wherever it used
to write -- check there first if a run seems to have vanished. (Separately,
and pre-existing: check_status.sh and this README still reference
TARGETED_PANEL_RUN_DIR / ~/workspace/KD/..., which this script has never
actually written to in its current form -- that mismatch predates this
task and is unrelated to the bf16_bench redirect; flagged, not fixed here.)

PAIRED-ARM SUPPORT (added 2026-09-22, card isp-runner-paired-arms-20260922,
for Pam's pre-registered S100 design, hive/reports/s100-isp-design-20260922.md):
three gaps her CPU preflight found, all fixed here, none touching
compare_runs.py's scoring definition:

1. Delete and overexpress used to run on DIFFERENT cell sets (delete
   auto-filters to token-positive cells inside InSilicoPerturber; overexpress
   does not -- see the module docstring above, still true, that is real
   library behavior, not a bug). Both operations now run on the same
   deterministic, pre-filtered, per-donor-capped, length-sorted cell list,
   built once per (gene, source) by paired_eligible_dataset() and passed as
   input_data_file to both -- delete's own internal token-filter becomes a
   harmless no-op on an already-filtered set, and overexpress (never
   filtered by the library) now only ever sees the same set delete does.
2. Donor identity used to be reconstructed after the fact (donor_consistency.py)
   by replicating InSilicoPerturber's internal filter+sort order on a cached
   copy of the source dataset -- fragile, and silently wrong if that
   replication ever drifted from the library's real order. Fixed by
   construction instead of reconstruction: paired_eligible_dataset() sorts
   the eligible cells by "length" descending itself (the exact op
   perturber_utils.downsample_and_sort applies when max_ncells=None -- source
   verified) before saving, so pyarrow's stable sort_indices makes the
   library's own internal re-sort of an already-identically-sorted dataset
   idempotent. The manifest CSV written alongside records donor per row in
   that exact final order, so it lines up positionally with the raw
   per-cell shift list without ever needing to re-derive anything.
3. No-op (unperturbed) inference did not exist. run_noop_gene() adds it:
   two independent get_embs() forward passes over the identical eligible
   cell list (no token edits at all), fed through the same
   perturber_utils.quant_cos_sims()/cos_sim_shift() the library itself uses
   for the real delete/overexpress path, written via the same
   write_perturbation_dictionary() so InSilicoPerturberStats.get_stats()
   consumes it identically to a real arm. A genuinely nonzero result here
   reflects real forward-pass numerical noise (batching, bf16 rounding,
   kernel nondeterminism), not an algebraic identity -- a naive
   "reuse-the-same-tensor-twice" no-op would be forced to exactly 0.0 and
   would tell us nothing.

Eligibility (>=50 cells from >=3 donors) is checked once per (gene, source)
in paired_eligible_dataset() and shared across delete/overexpress/noop; a
failure is recorded as "not_estimable" in that gene's completion marker,
never silently skipped or scored as zero. Per-donor cap is 100 cells, seed
20260922 (SEED below), a donor with fewer eligible cells keeps all of them.
This is the runner's own gate for the cell lists it builds -- separate from,
and upstream of, the full per-module/per-state donor-count and classifier
gates in Pam's design, which remain a preflight/analysis-level concern.
"""

from __future__ import annotations

import argparse
import multiprocessing

# Must happen before torch/CUDA is touched anywhere below: the isp_perturb_set*
# code path (taken whenever genes_to_perturb is a specific list, as it always
# is here) calls Dataset.map(num_proc=...), which forks worker processes even
# at num_proc=1. A fork after CUDA is initialized in this process crashes
# ("Cannot re-initialize CUDA in forked subprocess"). Switching the default
# start method to spawn avoids inheriting the initialized CUDA context.
multiprocessing.set_start_method("spawn", force=True)

import hashlib
import json
import logging
import os
import pickle
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from datasets import load_from_disk

HOME = Path.home()
# Pre-existing bug fixed here (2026-09-17, discovered via the bf16-bench
# calibration run): this was `.parents[1]`, which resolves to
# perturbation_workflow/ -- but target_gene_panel.json lives next to this
# script, in targeted_panel/, i.e. `.parents[0]`. Every other ANALYSIS_ROOT
# consumer (RAW_ROOT/STATS_ROOT/TABLE_ROOT/LOG_ROOT/SOURCE_DATA_DIR) was
# already redirected under bf16_bench/ by this task, so TARGET_GENES_FILE
# is the only remaining use and this fix is isolated to it.
ANALYSIS_ROOT = Path(__file__).resolve().parents[0]
ALLGENE_ROOT = HOME / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"
# Pre-existing bug fixed here (2026-09-18, discovered via a 316M-stage
# calibration run silently loading the 104M classifier): this was a bare
# hardcoded path, so `HTAN_FINETUNE_ROOT` had no effect here even though
# run_t4_overexpression.py already honors it for exactly this purpose
# (pointing ISP at a throwaway/alternate classifier checkpoint without
# touching the default /srv/lab/KD-backed path). Same env-var name, same
# default, now consistent between both runners.
FINETUNE_ROOT = Path(
    os.environ.get("HTAN_FINETUNE_ROOT", str(HOME / "workspace/KD/sclc_luad_normal_htan_finetune"))
)
BF16_BENCH_ROOT = Path(__file__).resolve().parents[2] / "bf16_bench"

sys.path.insert(0, str(BF16_BENCH_ROOT))
from dtype_cast import DTYPES, install_dtype_cast  # noqa: E402

# Overridable so the bf16-bench arms can point at a specific pinned checkout
# (e.g. the private f45a6c7 copy) instead of the shared geneformer-workspace
# symlink target -- same GENEFORMER_ROOT override run_t4_overexpression.py
# already supports.
GENEFORMER_ROOT = Path(
    os.environ.get(
        "GENEFORMER_ROOT",
        HOME / "workspace/geneformer-uv-starter/geneformer-workspace/Geneformer",
    )
)
sys.path.insert(0, str(GENEFORMER_ROOT))
from geneformer import TOKEN_DICTIONARY_FILE, InSilicoPerturber, InSilicoPerturberStats  # noqa: E402
from geneformer import perturber_utils as pu  # noqa: E402
from geneformer.emb_extractor import get_embs  # noqa: E402


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dtype", choices=DTYPES, default="fp32",
                   help="Model dtype for the ISP forward pass (default: fp32).")
    p.add_argument("--run-tag", default=None,
                   help="Output subdirectory under sclc_validation/bf16_bench/runs/ "
                        "(default: same as --dtype). Distinguishes arms that share a "
                        "dtype, e.g. the fp32 baseline vs the fp32 noise-floor repeat, "
                        "from clobbering each other's outputs.")
    p.add_argument("--force", action="store_true")
    p.add_argument("--max-ncells", type=int, default=None,
                   help="Cap on cells perturbed per source (InSilicoPerturber's "
                        "max_ncells). Sizing lever: apply the SAME value across every "
                        "arm (fp32 baseline, fp32 repeat, bf16) so the comparison stays "
                        "valid -- the fp32-vs-fp32 noise floor then quantifies exactly "
                        "the extra noise the cap adds. Default: no cap (all cells).")
    p.add_argument("--perturb-types", nargs="+", choices=("delete", "overexpress", "noop"),
                   default=["delete", "overexpress"],
                   help="Which perturb types to run (default: both delete and "
                        "overexpress; 'noop' is opt-in, not part of the default "
                        "since it doubles forward-pass cost for a diagnostic, not "
                        "a scored arm). Second sizing lever: pass "
                        "'--perturb-types overexpress' to drop delete entirely "
                        "(see the plan's dated 2026-09-17 sizing amendment).")
    return p.parse_args()


def geneformer_provenance() -> dict:
    try:
        commit = subprocess.run(
            ["git", "-C", str(GENEFORMER_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
    except Exception as exc:  # pragma: no cover - diagnostic only
        commit = f"<unresolved: {exc}>"
    return {"geneformer_root": str(GENEFORMER_ROOT), "geneformer_commit": commit}


def runner_provenance() -> dict:
    """This script's own commit hash AND its own SHA-256, for the run
    manifest -- required by Pam's design (isp-runner-paired-arms-20260922)
    so a run can be traced back to the exact runner code that produced it,
    same spirit as geneformer_provenance() above but for this repo, not the
    Geneformer checkout.

    The SHA-256 is computed here, of this file's own bytes ON DISK AT RUN
    TIME (2026-09-23, added after three separate stale-hash incidents in one
    evening during s100-isp-execution-20260922 gate reviews -- every one of
    them was a hash written into a document *before* a later merge changed
    this file, then trusted without recomputing). A hash recorded ahead of
    a run is a prediction about a file under active development; a hash the
    runner computes of itself as it runs cannot be stale by construction --
    whatever this run actually executed is exactly what got hashed. Preflight
    documents should still record their best-known hash for planning, but
    the run manifest's copy is the one that matters, because nothing has to
    trust it after the fact."""
    try:
        commit = subprocess.run(
            ["git", "-C", str(ANALYSIS_ROOT), "rev-parse", "HEAD"],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        dirty = subprocess.run(
            ["git", "-C", str(ANALYSIS_ROOT), "status", "--porcelain", "--", str(Path(__file__).name)],
            capture_output=True, text=True, check=True,
        ).stdout.strip() != ""
    except Exception as exc:  # pragma: no cover - diagnostic only
        commit, dirty = f"<unresolved: {exc}>", None
    try:
        runner_sha256 = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    except Exception as exc:  # pragma: no cover - diagnostic only
        runner_sha256 = f"<unresolved: {exc}>"
    return {
        "runner_file": str(Path(__file__).resolve()),
        "runner_commit": commit,
        "runner_file_dirty": dirty,
        "runner_sha256": runner_sha256,
    }


ARGS = parse_args()
install_dtype_cast(ARGS.dtype)
RUN_TAG = ARGS.run_tag or ARGS.dtype
OUT_ROOT = BF16_BENCH_ROOT / "runs" / RUN_TAG / "targeted_panel"

MODEL_PATH_FILE = FINETUNE_ROOT / "runs" / "MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt"
TRAIN_DATASET = ALLGENE_ROOT / "data/train_reference.dataset"
TEST_DATASET = ALLGENE_ROOT / "data/heldout_test.dataset"
# Overridable (2026-09-18): this pickle is a fixed per-disease-state
# centroid embedding computed once from a SPECIFIC model's embedding
# space (originally the 104M classifier, 768-dim). It is not portable
# across model architectures -- loading it against a 316M model's live
# embeddings (1152-dim) crashes with a cosine_similarity shape mismatch
# (confirmed live: "size of tensor a (1152) must match ... b (768)").
# Each model size needs its own centroids file; point this at a
# 316M-specific one instead of touching ALLGENE_ROOT (shared).
STATE_EMB_FILE = Path(
    os.environ.get(
        "STATE_EMB_FILE_OVERRIDE",
        str(ALLGENE_ROOT / "state_embeddings/training_donor_disease_centroids.pkl"),
    )
)

# Overridable (2026-09-18) so calibration/subset runs never need to edit
# the shared tracked target_gene_panel.json in place -- that file lives in
# the live ts1 checkout other agents share, and a temporary in-place edit
# (even one immediately restored via `git checkout --`) is unnecessary risk
# once an override exists. Point this at a throwaway single-gene copy under
# bf16_bench/ instead.
TARGET_GENES_FILE = Path(
    os.environ.get("TARGET_GENES_FILE_OVERRIDE", str(ANALYSIS_ROOT / "target_gene_panel.json"))
)
RAW_ROOT = OUT_ROOT / "raw"
STATS_ROOT = OUT_ROOT / "stats"
TABLE_ROOT = OUT_ROOT / "tables"
LOG_ROOT = OUT_ROOT / "logs"

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
STATES = (SCLC, LUAD, NORMAL)
SLUGS = {SCLC: "sclc", LUAD: "luad", NORMAL: "normal"}
STATE_BY_SLUG = {v: k for k, v in SLUGS.items()}
DEFAULT_SOURCE_ORDER = ("normal", "sclc", "luad")
# Full valid set. ACTIVE_PERTURB_TYPES (below, from --perturb-types) is what
# ensure_dirs()/main() actually iterate -- kept separate so a run scoped to
# just "overexpress" doesn't create/expect "delete" dirs at all.
PERTURB_TYPES = ("delete", "overexpress")
ACTIVE_PERTURB_TYPES = tuple(ARGS.perturb_types)
FORWARD_BATCH_SIZE = 128
# Pre-existing bug fixed here (2026-09-17, discovered via the bf16-bench
# calibration run): this was 4. InSilicoPerturber.perturb_data() loads the
# model onto CUDA before Dataset.map(num_proc=...) runs, and Geneformer's
# map() uses the separate `multiprocess` package, which forks workers
# regardless of this file's own multiprocessing.set_start_method("spawn")
# call (that only affects stdlib multiprocessing, a different global than
# `multiprocess`). Forking after CUDA init crashes with "Cannot
# re-initialize CUDA in forked subprocess" -- confirmed by reproducing it
# live with nproc=4. run_t4_overexpression.py already documents this exact
# failure mode and sets nproc=1 for GPU runs for exactly this reason; this
# script never got the same fix, so every gene here was silently unable to
# run against a live GPU.
NPROC = 1


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def model_dir() -> Path:
    return Path(MODEL_PATH_FILE.read_text().strip())


def ensure_dirs() -> None:
    dirs = [TABLE_ROOT, LOG_ROOT]
    for ptype in ACTIVE_PERTURB_TYPES:
        for slug in SLUGS.values():
            dirs.append(RAW_ROOT / ptype / slug)
        dirs.append(STATS_ROOT / ptype)
    for d in dirs:
        d.mkdir(parents=True, exist_ok=True)


def load_target_genes() -> list[dict]:
    payload = json.loads(TARGET_GENES_FILE.read_text())
    return payload["genes"]


def canonical_states(start_state: str) -> dict:
    others = [s for s in STATES if s != start_state]
    return {"state_key": "disease", "start_state": start_state, "goal_state": others[0], "alt_states": [others[1]]}


# Dtype-independent input cache (same held-out cells for every arm), shared
# across arms rather than duplicated per --run-tag; lives under bf16_bench
# per the task boundary that new outputs stay out of the pre-existing tree.
SOURCE_DATA_DIR = BF16_BENCH_ROOT / "runs" / "_shared" / "targeted_panel_sources"

# Paired-arm design constants (2026-09-22, Pam's S100 design). Dtype- and
# run-tag-independent, same as SOURCE_DATA_DIR above: the eligible cell list
# for a given (gene, source) is a property of the data and the seed, not of
# which dtype arm is consuming it, so it is built once and shared.
SEED = 20260922
DONOR_CELL_CAP = 100
MIN_CELLS_ELIGIBLE = 50
MIN_DONORS_ELIGIBLE = 3
PAIRED_ELIGIBLE_DIR = BF16_BENCH_ROOT / "runs" / "_shared" / "targeted_panel_paired_eligible"
MANIFEST_DIR = PAIRED_ELIGIBLE_DIR / "manifests"


def source_dataset_path(source_slug: str) -> Path:
    """Materialize (once) the held-out test cells for one source disease to
    disk, since InSilicoPerturber.perturb_data requires a path, not an
    in-memory Dataset."""
    path = SOURCE_DATA_DIR / f"{source_slug}.dataset"
    if not path.exists():
        # num_proc=1: avoid forking after CUDA may already be initialized in
        # this process (see NPROC note above) -- callers should pre-materialize
        # every source's dataset before the first perturb_data() call.
        SOURCE_DATA_DIR.mkdir(parents=True, exist_ok=True)
        test = load_from_disk(str(TEST_DATASET))
        disease = STATE_BY_SLUG[source_slug]
        subset = test.filter(lambda row, d=disease: row["disease"] == d, num_proc=1)
        subset.save_to_disk(str(path))
    return path


_GENE_TOKEN_DICT: dict | None = None


def gene_token_dict() -> dict:
    """Ensembl ID -> token, loaded once. Same file InSilicoPerturber itself
    loads internally (self.token_dictionary_file defaults to this); loading
    it directly here (rather than instantiating a throwaway InSilicoPerturber
    just to read its .gene_token_dict) avoids the confusion of a perturber
    object that is never used for perturbation."""
    global _GENE_TOKEN_DICT
    if _GENE_TOKEN_DICT is None:
        with open(TOKEN_DICTIONARY_FILE, "rb") as f:
            _GENE_TOKEN_DICT = pickle.load(f)
    return _GENE_TOKEN_DICT


def _seeded_rng(*parts: str) -> np.random.Generator:
    """Deterministic per-(gene, source, donor) RNG derived from SEED and the
    given parts, via a hash rather than global generator state -- so the
    donor cap draw for one (gene, source, donor) never depends on the order
    other genes/sources/donors happen to be processed in."""
    digest = hashlib.sha256(":".join([str(SEED), *parts]).encode()).hexdigest()
    return np.random.default_rng(int(digest[:16], 16))


def paired_eligible_dataset(source_slug: str, gene: dict) -> tuple[Path | None, dict]:
    """The single deterministic, token-positive, per-donor-capped, length-
    sorted cell list for this (gene, source) -- built once and reused as-is
    for delete, overexpress, and noop, so all three run on identical cells
    (gap 1) and the saved manifest's donor column lines up positionally with
    every arm's raw per-cell output (gap 2). See the module docstring for
    why the length-descending sort here makes InSilicoPerturber's own
    internal re-sort (downsample_and_sort, always applied, a no-op reorder
    when max_ncells=None and the input is already sorted the same way)
    order-preserving rather than a second, unaccounted-for reshuffle.

    Returns (dataset_path, info). dataset_path is None if the eligibility
    gate fails (info["eligible"] is then False and info["reason"] explains
    why) -- callers must treat that as "not estimable", never as zero cells
    processed.
    """
    ensembl_id = gene["ensembl_id"]
    key = f"{source_slug}_{ensembl_id}"
    ds_path = PAIRED_ELIGIBLE_DIR / f"{key}.dataset"
    info_path = PAIRED_ELIGIBLE_DIR / f"{key}.eligibility.json"
    if info_path.exists():
        info = json.loads(info_path.read_text())
        return (ds_path if info["eligible"] else None), info

    PAIRED_ELIGIBLE_DIR.mkdir(parents=True, exist_ok=True)
    MANIFEST_DIR.mkdir(parents=True, exist_ok=True)

    token = gene_token_dict().get(ensembl_id)
    if token is None:
        info = {"eligible": False, "reason": f"{ensembl_id} not in token dictionary",
                "n_cells": 0, "n_donors": 0}
        info_path.write_text(json.dumps(info, indent=2) + "\n")
        return None, info

    full = load_from_disk(str(source_dataset_path(source_slug)))
    # Same criterion as InSilicoPerturber's own delete-path pre-filter
    # (perturber_utils.filter_data_by_tokens), applied here to BOTH
    # operations up front instead of relying on the library to apply it
    # (which it only does for delete -- see module docstring).
    token_positive = full.filter(lambda row: token in row["input_ids"], num_proc=1)

    donors = sorted(set(token_positive["individual"]))
    if len(token_positive) < MIN_CELLS_ELIGIBLE or len(donors) < MIN_DONORS_ELIGIBLE:
        info = {
            "eligible": False,
            "reason": f"{len(token_positive)} token-positive cells from {len(donors)} donors "
                      f"(need >= {MIN_CELLS_ELIGIBLE} cells from >= {MIN_DONORS_ELIGIBLE} donors)",
            "n_cells": len(token_positive), "n_donors": len(donors),
        }
        info_path.write_text(json.dumps(info, indent=2) + "\n")
        return None, info

    # Deterministic per-donor cap: keep every cell for a donor at or under
    # the cap, otherwise draw exactly DONOR_CELL_CAP without replacement
    # using a seed derived only from (SEED, source, gene, donor) -- never
    # topped up from another donor.
    keep_indices: list[int] = []
    donor_cell_counts: dict[str, int] = {}
    individuals = token_positive["individual"]
    for donor in donors:
        donor_idx = [i for i, ind in enumerate(individuals) if ind == donor]
        donor_cell_counts[donor] = len(donor_idx)
        if len(donor_idx) > DONOR_CELL_CAP:
            rng = _seeded_rng(source_slug, ensembl_id, donor)
            donor_idx = sorted(rng.choice(donor_idx, size=DONOR_CELL_CAP, replace=False).tolist())
        keep_indices.extend(donor_idx)

    capped = token_positive.select(sorted(keep_indices))
    # Sort by length descending ourselves -- see module docstring: this is
    # the exact transform perturber_utils.downsample_and_sort applies (with
    # max_ncells=None, no shuffle), so InSilicoPerturber's own internal
    # re-sort of this already-identically-sorted dataset is idempotent
    # (pyarrow's sort_indices is stable; verified against the installed
    # datasets==5.0.1 source, not assumed).
    sorted_ds = capped.sort("length", reverse=True)

    # cell_id/individual/length confirmed present on the real heldout_test.dataset
    # (checked directly, not assumed: ['input_ids', 'cell_id', 'individual',
    # 'celltype', 'disease', 'split', 'length']).
    manifest = pd.DataFrame({
        "cell_id": sorted_ds["cell_id"],
        "donor": sorted_ds["individual"],
        "length": sorted_ds["length"],
    })
    manifest_path = MANIFEST_DIR / f"{key}.csv"
    manifest.to_csv(manifest_path, index=False)

    sorted_ds.save_to_disk(str(ds_path))
    info = {
        "eligible": True,
        "n_cells": len(sorted_ds),
        "n_donors": len(donors),
        "donor_cell_counts": donor_cell_counts,
        "seed": SEED,
        "donor_cell_cap": DONOR_CELL_CAP,
        "manifest": str(manifest_path),
    }
    info_path.write_text(json.dumps(info, indent=2) + "\n")
    return ds_path, info


def state_embeddings() -> dict:
    assert STATE_EMB_FILE.exists(), f"Missing {STATE_EMB_FILE} -- run the all-gene prepare/state-embeddings stages first"
    with STATE_EMB_FILE.open("rb") as f:
        return pickle.load(f)


def _not_estimable_marker(raw_dir: Path, prefix: str, perturb_type: str, source_slug: str,
                           symbol: str, ensembl_id: str, elig: dict) -> None:
    payload = {
        "completed_utc": utc_now(),
        "perturb_type": perturb_type,
        "source": source_slug,
        "gene": symbol,
        "ensembl_id": ensembl_id,
        "dtype": ARGS.dtype,
        "status": "not_estimable",
        "eligibility": elig,
    }
    (raw_dir / f"{prefix}.complete.json").write_text(json.dumps(payload, indent=2) + "\n")
    logging.warning(
        "[%s/%s/%s] not estimable: %s -- recorded as not_estimable, NOT run, NOT a zero effect",
        perturb_type, source_slug, symbol, elig.get("reason"),
    )


def run_gene(perturb_type: str, source_slug: str, gene: dict, state_embs: dict, force: bool = False) -> None:
    disease = STATE_BY_SLUG[source_slug]
    raw_dir = RAW_ROOT / perturb_type / source_slug
    symbol, ensembl_id = gene["gene"], gene["ensembl_id"]
    prefix = f"targeted_{source_slug}_{symbol}"
    done_file = raw_dir / f"{prefix}.complete.json"
    if done_file.exists() and not force:
        logging.info("[%s/%s/%s] already complete", perturb_type, source_slug, symbol)
        return
    for partial in raw_dir.glob(f"in_silico_{perturb_type}_{prefix}_*_raw.pickle"):
        partial.unlink()

    data_path, elig = paired_eligible_dataset(source_slug, gene)
    if not elig["eligible"]:
        _not_estimable_marker(raw_dir, prefix, perturb_type, source_slug, symbol, ensembl_id, elig)
        return
    started = time.time()
    perturber = InSilicoPerturber(
        perturb_type=perturb_type,
        genes_to_perturb=[ensembl_id],
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",
        filter_data=None,
        cell_states_to_model=canonical_states(disease),
        state_embs_dict=state_embs,
        max_ncells=ARGS.max_ncells,
        emb_layer=0,
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=NPROC,
        model_version="V2",
        clear_mem_ncells=1000,
    )
    try:
        perturber.perturb_data(
            model_directory=str(model_dir()),
            input_data_file=str(data_path),
            output_directory=str(raw_dir),
            output_prefix=prefix,
        )
    except RuntimeError as exc:
        # Geneformer's filter_data_by_tokens_and_log has a bare `raise` (no
        # active exception) when zero cells in the source contain the target
        # gene -- this is a real, expected outcome for e.g. SCLC-subtype
        # transcription factors (ASCL1/NEUROD1/POU2F3) tested against normal
        # T cells, not a bug in this script. Treat it as "0 cells detected"
        # and move on rather than crashing the whole sweep.
        if "No active exception to reraise" in str(exc):
            logging.warning(
                "[%s/%s/%s] 0 cells detected -- skipping (expected for genes not "
                "biologically active in this source's cell type)",
                perturb_type, source_slug, symbol,
            )
            payload = {
                "completed_utc": utc_now(),
                "perturb_type": perturb_type,
                "source": source_slug,
                "gene": symbol,
                "ensembl_id": ensembl_id,
                "dtype": ARGS.dtype,
                "elapsed_seconds": time.time() - started,
                "n_raw_files": 0,
                "skipped_zero_cells_detected": True,
            }
            done_file.write_text(json.dumps(payload, indent=2) + "\n")
            return
        raise

    output_files = sorted(raw_dir.glob(f"in_silico_{perturb_type}_{prefix}_*_raw.pickle"))
    peak_mem_gib = None
    if torch.cuda.is_available():
        peak_mem_gib = torch.cuda.max_memory_allocated() / 2**30
        torch.cuda.reset_peak_memory_stats()
    payload = {
        "completed_utc": utc_now(),
        "perturb_type": perturb_type,
        "source": source_slug,
        "gene": symbol,
        "ensembl_id": ensembl_id,
        "dtype": ARGS.dtype,
        "elapsed_seconds": time.time() - started,
        "n_raw_files": len(output_files),
        "peak_gpu_mem_gib": peak_mem_gib,
    }
    done_file.write_text(json.dumps(payload, indent=2) + "\n")
    logging.info(
        "[%s/%s/%s] complete in %.1f s (%d raw files)",
        perturb_type, source_slug, symbol, payload["elapsed_seconds"], len(output_files),
    )


def run_noop_gene(source_slug: str, gene: dict, state_embs: dict, force: bool = False) -> None:
    """Unperturbed control: two independent forward passes over the exact
    same eligible cell list used by delete/overexpress for this gene, no
    token edits at all. Reuses perturber_utils.quant_cos_sims/cos_sim_shift
    (the library's own scoring primitive, unmodified) and
    write_perturbation_dictionary (the library's own pickle writer) so
    InSilicoPerturberStats.get_stats() and run_stats() below consume this
    identically to a real arm, with zero changes to either. See the module
    docstring for why two independent passes (not the same embedding reused
    twice, which would be forced to exactly 0.0) is the only version of this
    that actually tests numerical stability."""
    disease = STATE_BY_SLUG[source_slug]
    raw_dir = RAW_ROOT / "noop" / source_slug
    symbol, ensembl_id = gene["gene"], gene["ensembl_id"]
    prefix = f"targeted_{source_slug}_{symbol}"
    done_file = raw_dir / f"{prefix}.complete.json"
    if done_file.exists() and not force:
        logging.info("[noop/%s/%s] already complete", source_slug, symbol)
        return

    data_path, elig = paired_eligible_dataset(source_slug, gene)
    if not elig["eligible"]:
        _not_estimable_marker(raw_dir, prefix, "noop", source_slug, symbol, ensembl_id, elig)
        return
    started = time.time()

    # Constructed only to reuse its validated setup (token_gene_dict,
    # pad_token_id, tokens_to_perturb) -- perturb_data() is never called on
    # it, so perturb_type is a required-but-inert placeholder here.
    perturber = InSilicoPerturber(
        perturb_type="delete",
        genes_to_perturb=[ensembl_id],
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",
        filter_data=None,
        cell_states_to_model=canonical_states(disease),
        state_embs_dict=state_embs,
        max_ncells=None,
        emb_layer=0,
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=NPROC,
        model_version="V2",
        clear_mem_ncells=1000,
    )
    model = pu.load_model(perturber.model_type, perturber.num_classes, str(model_dir()), mode="eval")
    layer_to_quant = pu.quant_layers(model) + perturber.emb_layer

    ds = load_from_disk(str(data_path))
    token = perturber.tokens_to_perturb[0]
    possible_states = pu.get_possible_states(perturber.cell_states_to_model)
    accum = {state: [] for state in possible_states}
    for i in range(0, len(ds), FORWARD_BATCH_SIZE):
        minibatch = ds.select(range(i, min(i + FORWARD_BATCH_SIZE, len(ds))))
        pass1 = get_embs(model, minibatch, "cls", layer_to_quant, perturber.pad_token_id,
                          FORWARD_BATCH_SIZE, token_gene_dict=perturber.token_gene_dict,
                          summary_stat=None, silent=True)
        pass2 = get_embs(model, minibatch, "cls", layer_to_quant, perturber.pad_token_id,
                          FORWARD_BATCH_SIZE, token_gene_dict=perturber.token_gene_dict,
                          summary_stat=None, silent=True)
        # perturbed_emb, original_emb order matches quant_cos_sims'/cos_sim_shift's
        # own contract (see the module docstring): pass1 stands in for
        # "original", pass2 for "perturbed" -- an arbitrary but fixed
        # labeling since both passes are over identical, unperturbed input.
        cos_sims = pu.quant_cos_sims(pass2, pass1, perturber.cell_states_to_model, state_embs, emb_mode="cell")
        for state, values in cos_sims.items():
            accum[state].extend(values.tolist())

    cos_sims_dict = {state: {(token, "cell_emb"): values} for state, values in accum.items()}

    output_path_prefix = os.path.join(str(raw_dir), f"in_silico_noop_{prefix}")
    pu.write_perturbation_dictionary(cos_sims_dict, output_path_prefix)

    if torch.cuda.is_available():
        peak_mem_gib = torch.cuda.max_memory_allocated() / 2**30
        torch.cuda.reset_peak_memory_stats()
    else:
        peak_mem_gib = None
    payload = {
        "completed_utc": utc_now(),
        "perturb_type": "noop",
        "source": source_slug,
        "gene": symbol,
        "ensembl_id": ensembl_id,
        "dtype": ARGS.dtype,
        "elapsed_seconds": time.time() - started,
        "n_cells": len(ds),
        "n_raw_files": 1,
        "peak_gpu_mem_gib": peak_mem_gib,
    }
    done_file.write_text(json.dumps(payload, indent=2) + "\n")
    logging.info("[noop/%s/%s] complete in %.1f s (%d cells, 2 forward passes)",
                 source_slug, symbol, payload["elapsed_seconds"], len(ds))


def run_stats(perturb_type: str, source_slug: str, target_genes: list[dict], force: bool = False) -> None:
    source_state = STATE_BY_SLUG[source_slug]
    raw_dir = RAW_ROOT / perturb_type / source_slug
    n_genes = len(target_genes)
    completed = len(list(raw_dir.glob("targeted_*.complete.json")))
    if completed != n_genes:
        raise RuntimeError(f"Cannot run stats for {perturb_type}/{source_slug}: {completed}/{n_genes} genes complete")
    # Every gene can legitimately be not_estimable for a source (e.g. "normal"
    # has a single donor -- MIN_DONORS_ELIGIBLE=3 can never be met there,
    # regardless of gene), in which case *_raw.pickle never gets written for
    # any gene even though every gene's .complete.json exists (not_estimable
    # is itself a completion state -- see _not_estimable_marker). Without this
    # check, InSilicoPerturberStats.get_stats() -> read_dictionaries() finds
    # zero matching pickles and crashes with an opaque bare `raise`
    # ("RuntimeError: No active exception to reraise") -- a real library
    # behavior, not a bug in this script, but one this script must guard
    # against rather than let surface as a stats-phase crash after all GPU
    # work for the run is already done. Confirmed live (2026-09-22, S100 ISP
    # execution gate 2): source=normal is not_estimable for every gene by
    # construction, so this path is guaranteed to trigger on every real run
    # that includes it.
    if not any(raw_dir.glob(f"in_silico_{perturb_type}_targeted_{source_slug}_*_raw.pickle")):
        logging.warning(
            "[%s/%s] no raw output for any of %d gene(s) (all not_estimable) -- "
            "skipping stats for this source entirely, not a crash",
            perturb_type, source_slug, n_genes,
        )
        return
    targets = [s for s in STATES if s != source_state]
    stats_dir = STATS_ROOT / perturb_type
    for target in targets:
        alt = [s for s in STATES if s not in (source_state, target)]
        comparison = {"state_key": "disease", "start_state": source_state, "goal_state": target, "alt_states": alt}
        target_slug = SLUGS[target]
        output_prefix = f"targeted_{perturb_type}_{source_slug}_to_{target_slug}"
        output_file = stats_dir / f"{output_prefix}.csv"
        if output_file.exists() and not force:
            logging.info("Stats already exist: %s", output_file)
            continue
        stats = InSilicoPerturberStats(
            mode="goal_state_shift",
            genes_perturbed="all",
            combos=0,
            anchor_gene=None,
            cell_states_to_model=comparison,
            model_version="V2",
        )
        stats.get_stats(
            input_data_directory=str(raw_dir),
            null_dist_data_directory=None,
            output_directory=str(stats_dir),
            output_prefix=output_prefix,
        )
        logging.info("Wrote %s", output_file)


def main() -> None:
    logging.basicConfig(level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s")
    ensure_dirs()
    target_genes = load_target_genes()
    (TABLE_ROOT / "run_config.json").write_text(json.dumps({
        "created_utc": utc_now(),
        "model_directory": str(model_dir()),
        "n_target_genes": len(target_genes),
        "perturbation_types": list(ACTIVE_PERTURB_TYPES),
        "forward_batch_size": FORWARD_BATCH_SIZE,
        "nproc": NPROC,
        "source_order": list(DEFAULT_SOURCE_ORDER),
        "target_genes_file": str(TARGET_GENES_FILE),
        "dtype": ARGS.dtype,
        "run_tag": RUN_TAG,
        "max_ncells": ARGS.max_ncells,
        "paired_arm_seed": SEED,
        "paired_arm_donor_cell_cap": DONOR_CELL_CAP,
        "paired_arm_min_cells_eligible": MIN_CELLS_ELIGIBLE,
        "paired_arm_min_donors_eligible": MIN_DONORS_ELIGIBLE,
        **geneformer_provenance(),
        **runner_provenance(),
    }, indent=2) + "\n")

    if ARGS.max_ncells is not None and any(pt in ("delete", "overexpress", "noop") for pt in ACTIVE_PERTURB_TYPES):
        logging.warning(
            "--max-ncells=%d is set together with the paired-arm path (delete/overexpress/noop). "
            "The paired eligible dataset is already capped deterministically at %d cells/donor "
            "(seed %d). If its total size still exceeds --max-ncells, InSilicoPerturber's own "
            "downsample_and_sort will shuffle(seed=42) on top of that -- a DIFFERENT, hardcoded "
            "library seed, not this run's seed -- reintroducing non-reproducibility the paired "
            "design is meant to avoid. Leave --max-ncells unset for a design-conformant paired run; "
            "it remains available for the older bf16-bench sizing arms that do not use this path.",
            ARGS.max_ncells, DONOR_CELL_CAP, SEED,
        )

    state_embs = state_embeddings()
    logging.info("Loaded state embeddings for: %s", sorted(state_embs))

    # Materialize every source's held-out dataset before any perturb_data()
    # call touches CUDA (see source_dataset_path note).
    for source in DEFAULT_SOURCE_ORDER:
        path = source_dataset_path(source)
        logging.info("Source dataset ready: %s -> %s", source, path)

    for ptype in ACTIVE_PERTURB_TYPES:
        for source in DEFAULT_SOURCE_ORDER:
            logging.info("=== %s / %s: %d genes ===", ptype, source, len(target_genes))
            for gene in target_genes:
                if ptype == "noop":
                    run_noop_gene(source, gene, state_embs, force=ARGS.force)
                else:
                    run_gene(ptype, source, gene, state_embs, force=ARGS.force)

    for ptype in ACTIVE_PERTURB_TYPES:
        for source in DEFAULT_SOURCE_ORDER:
            run_stats(ptype, source, target_genes, force=ARGS.force)

    logging.info("Targeted panel perturbation complete.")


if __name__ == "__main__":
    main()
