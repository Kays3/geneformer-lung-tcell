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
"""

from __future__ import annotations

import multiprocessing

# Must happen before torch/CUDA is touched anywhere below: the isp_perturb_set*
# code path (taken whenever genes_to_perturb is a specific list, as it always
# is here) calls Dataset.map(num_proc=...), which forks worker processes even
# at num_proc=1. A fork after CUDA is initialized in this process crashes
# ("Cannot re-initialize CUDA in forked subprocess"). Switching the default
# start method to spawn avoids inheriting the initialized CUDA context.
multiprocessing.set_start_method("spawn", force=True)

import json
import logging
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

from datasets import load_from_disk

HOME = Path.home()
ANALYSIS_ROOT = Path(__file__).resolve().parents[1]
ALLGENE_ROOT = HOME / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"
FINETUNE_ROOT = HOME / "workspace/KD/sclc_luad_normal_htan_finetune"

sys.path.insert(0, str(HOME / "workspace/geneformer-uv-starter/geneformer-workspace/Geneformer"))
from geneformer import InSilicoPerturber, InSilicoPerturberStats

MODEL_PATH_FILE = FINETUNE_ROOT / "runs" / "MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt"
TRAIN_DATASET = ALLGENE_ROOT / "data/train_reference.dataset"
TEST_DATASET = ALLGENE_ROOT / "data/heldout_test.dataset"
STATE_EMB_FILE = ALLGENE_ROOT / "state_embeddings/training_donor_disease_centroids.pkl"

TARGET_GENES_FILE = ANALYSIS_ROOT / "target_gene_panel.json"
RAW_ROOT = ANALYSIS_ROOT / "raw"
STATS_ROOT = ANALYSIS_ROOT / "stats"
TABLE_ROOT = ANALYSIS_ROOT / "tables"
LOG_ROOT = ANALYSIS_ROOT / "logs"

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
STATES = (SCLC, LUAD, NORMAL)
SLUGS = {SCLC: "sclc", LUAD: "luad", NORMAL: "normal"}
STATE_BY_SLUG = {v: k for k, v in SLUGS.items()}
DEFAULT_SOURCE_ORDER = ("normal", "sclc", "luad")
PERTURB_TYPES = ("delete", "overexpress")
FORWARD_BATCH_SIZE = 128
NPROC = 4


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def model_dir() -> Path:
    return Path(MODEL_PATH_FILE.read_text().strip())


def ensure_dirs() -> None:
    dirs = [TABLE_ROOT, LOG_ROOT]
    for ptype in PERTURB_TYPES:
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


SOURCE_DATA_DIR = ANALYSIS_ROOT / "data" / "sources"


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


def state_embeddings() -> dict:
    assert STATE_EMB_FILE.exists(), f"Missing {STATE_EMB_FILE} -- run the all-gene prepare/state-embeddings stages first"
    with STATE_EMB_FILE.open("rb") as f:
        return pickle.load(f)


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

    data_path = source_dataset_path(source_slug)
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
        max_ncells=None,
        emb_layer=0,
        forward_batch_size=FORWARD_BATCH_SIZE,
        nproc=NPROC,
        model_version="V2",
        clear_mem_ncells=1000,
    )
    perturber.perturb_data(
        model_directory=str(model_dir()),
        input_data_file=str(data_path),
        output_directory=str(raw_dir),
        output_prefix=prefix,
    )
    output_files = sorted(raw_dir.glob(f"in_silico_{perturb_type}_{prefix}_*_raw.pickle"))
    payload = {
        "completed_utc": utc_now(),
        "perturb_type": perturb_type,
        "source": source_slug,
        "gene": symbol,
        "ensembl_id": ensembl_id,
        "elapsed_seconds": time.time() - started,
        "n_raw_files": len(output_files),
    }
    done_file.write_text(json.dumps(payload, indent=2) + "\n")
    logging.info(
        "[%s/%s/%s] complete in %.1f s (%d raw files)",
        perturb_type, source_slug, symbol, payload["elapsed_seconds"], len(output_files),
    )


def run_stats(perturb_type: str, source_slug: str, target_genes: list[dict], force: bool = False) -> None:
    source_state = STATE_BY_SLUG[source_slug]
    raw_dir = RAW_ROOT / perturb_type / source_slug
    n_genes = len(target_genes)
    completed = len(list(raw_dir.glob("targeted_*.complete.json")))
    if completed != n_genes:
        raise RuntimeError(f"Cannot run stats for {perturb_type}/{source_slug}: {completed}/{n_genes} genes complete")
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
        "perturbation_types": list(PERTURB_TYPES),
        "forward_batch_size": FORWARD_BATCH_SIZE,
        "nproc": NPROC,
        "source_order": list(DEFAULT_SOURCE_ORDER),
        "target_genes_file": str(TARGET_GENES_FILE),
    }, indent=2) + "\n")

    state_embs = state_embeddings()
    logging.info("Loaded state embeddings for: %s", sorted(state_embs))

    # Materialize every source's held-out dataset before any perturb_data()
    # call touches CUDA (see source_dataset_path note).
    for source in DEFAULT_SOURCE_ORDER:
        path = source_dataset_path(source)
        logging.info("Source dataset ready: %s -> %s", source, path)

    for ptype in PERTURB_TYPES:
        for source in DEFAULT_SOURCE_ORDER:
            logging.info("=== %s / %s: %d genes ===", ptype, source, len(target_genes))
            for gene in target_genes:
                run_gene(ptype, source, gene, state_embs)

    for ptype in PERTURB_TYPES:
        for source in DEFAULT_SOURCE_ORDER:
            run_stats(ptype, source, target_genes)

    logging.info("Targeted panel perturbation complete.")


if __name__ == "__main__":
    main()
