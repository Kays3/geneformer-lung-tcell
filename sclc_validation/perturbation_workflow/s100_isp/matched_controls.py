"""Matched-control table construction for the S100 ISP design
(hive/reports/s100-isp-design-20260922.md, "Controls and analysis", item 3).

STATUS: UNBLOCKED 2026-09-23 (human ruling, s100-isp-execution-20260922).
The two definitional questions that blocked this module are now ruled --
see retained_rows_spec_20260922.md's second dated amendment for the full
text. Summary:

1. "Median token rank" is not in ambient_risk_all_genes.csv (that file has
   detect_frac, not a rank column, and detect_frac must not stand in for
   it -- they are two independent matching axes). It is computed directly
   from the tokenized held-out dataset (ALLGENE_ROOT/data/heldout_test.dataset,
   the same source the eligibility counts come from): for each gene, the
   median, over cells in which that gene is token-positive, of the gene's
   0-based position in that cell's rank-ordered `input_ids`. See
   `median_token_rank_and_detection()` below.
2. Matching scope is PER-SOURCE-STATE, not global. With the SCLC source arm
   dropped (see the same amendment), this means: LUAD held-out cells only.
   Both matching quantities (detection fraction AND median token rank) are
   computed from that same LUAD-only cell set -- ambient_risk_all_genes.csv's
   detect_frac is a global (cross-source) figure and is NOT used here for
   that reason, not because detection fraction itself changed meaning.

`synthetic_control_table()` below remains for ambient_stats.py's own tests
(fabricated ids that cannot collide with anything real) -- it was never a
guess at the real answer, and nothing here retires it.
"""
from __future__ import annotations

import hashlib

import numpy as np
import pandas as pd
from scipy.stats import rankdata

from provenance_utils import sha256_file
from retained_rows import STATUS_NOT_ESTIMABLE_CONTROL_STRATUM

MIN_COMMON_CONTROLS = 20
SEED = 20260922

# The six fixed strata, exactly as registered in the design doc and as
# present in the "stratum" column of s100_gene_panel_20260922.json. Keyed
# by gene SYMBOL, matching the panel file's own grouping.
STRATA = {
    "clean_high": ("S100A4", "S100A6", "S100A10", "S100A11"),
    "low_p_a16": ("S100P", "S100A16"),
    "low_a2_b": ("S100A2", "S100B"),
    "singleton_a13": ("S100A13",),
    "singleton_a8": ("S100A8",),
    "low_a9_pbp": ("S100A9", "S100PBP"),
}


# ---------------------------------------------------------------------------
# Step 1: per-gene detection fraction + median token rank, from LUAD-only
# tokenized held-out cells.
# ---------------------------------------------------------------------------


def median_token_rank_and_detection(
    cells: list[list[int]],
    gene_token_dict: dict[str, int],
    rank_convention: str = "0_based",
) -> pd.DataFrame:
    """Detection fraction and median token rank per gene, over the given
    cells only -- caller filters to one source-state (LUAD) before calling
    this; this function has no notion of "source" at all, by design, so it
    cannot silently mix sources.

    `cells` is each cell's `input_ids`, already rank-ordered by the
    tokenizer (position 0 = the most highly expressed gene in that cell --
    real Geneformer convention, not re-derived here). One pass over all
    cells: for every token seen, record how many cells contain it (the
    detection-fraction numerator) and every position it was found at (the
    median-token-rank sample). Bounded by the total number of
    (cell, detected-gene) pairs in the corpus -- the same quantity
    ambient_risk_all_genes.csv's own detect_frac computation already
    processes over the whole corpus, just restricted to LUAD cells and
    extended with position.

    `rank_convention`: "0_based" (index in the list, Python's own
    convention -- the default) or "1_based" (index + 1). This choice does
    NOT affect matching (every tolerance check is a difference between two
    genes' values under the SAME convention), but it must be stated once
    and frozen, per the ruling -- an unstated convention is a future
    argument, not a live one.
    """
    if rank_convention not in ("0_based", "1_based"):
        raise ValueError(f"rank_convention must be '0_based' or '1_based', got {rank_convention!r}")
    offset = 1 if rank_convention == "1_based" else 0

    n_cells = len(cells)
    detect_count: dict[int, int] = {}
    positions: dict[int, list[int]] = {}
    for input_ids in cells:
        for pos, token in enumerate(input_ids):
            detect_count[token] = detect_count.get(token, 0) + 1
            positions.setdefault(token, []).append(pos + offset)

    rows = []
    for ensembl_id, token in gene_token_dict.items():
        count = detect_count.get(token, 0)
        rows.append({
            "ensembl_id": ensembl_id,
            "detect_frac": (count / n_cells) if n_cells else np.nan,
            "median_token_rank": float(np.median(positions[token])) if token in positions else np.nan,
            "n_positive_cells": count,
        })
    return pd.DataFrame(rows)


def load_and_freeze_luad_gene_stats(
    dataset_path,
    gene_token_dict: dict[str, int],
    out_path,
    *,
    disease_column: str = "disease",
    luad_disease_value: str = "lung adenocarcinoma",
    rank_convention: str = "0_based",
) -> tuple[pd.DataFrame, str]:
    """Load the real tokenized held-out dataset, filter to LUAD cells only,
    compute detect_frac/median_token_rank for every gene in
    `gene_token_dict`, write the frozen table to `out_path`, and hash it
    from its own bytes on disk (provenance_utils.sha256_file -- same
    principle as run_targeted_panel.py's runner/input hashing: a hash
    computed of the table as written cannot be stale by construction).

    Local import of `datasets` -- this is the one function in this module
    that touches the real corpus; everything else (median_token_rank_and_
    detection, build_matched_control_table) is plain pandas/numpy and
    testable without it.
    """
    from datasets import load_from_disk

    ds = load_from_disk(str(dataset_path))
    luad = ds.filter(lambda row: row[disease_column] == luad_disease_value, num_proc=1)
    stats = median_token_rank_and_detection(luad["input_ids"], gene_token_dict, rank_convention=rank_convention)
    stats.to_csv(out_path, index=False)
    return stats, sha256_file(out_path)


# ---------------------------------------------------------------------------
# Step 2: rank-percentile transform (for the 5-percentile-point tolerance)
# ---------------------------------------------------------------------------


def rank_percentile_transform(values: pd.Series) -> pd.Series:
    """Each gene's median_token_rank expressed as a percentile among the
    population of genes passed in (ties averaged) -- raw rank positions
    span thousands of tokens, so "5 median-rank-percentile-points" (the
    design doc's own tolerance unit) is meaningless without first
    normalizing to a percentile. Same midrank convention as
    ambient_stats.midrank_percentile(): rank via rankdata(method="average"),
    map rank r among N points to (r - 0.5) / N * 100."""
    ranks = rankdata(values, method="average")
    n = len(values)
    return pd.Series(100.0 * (ranks - 0.5) / n, index=values.index)


# ---------------------------------------------------------------------------
# Step 3: the matched-control table itself
# ---------------------------------------------------------------------------


def _seeded_rng(seed: int, *parts: str) -> np.random.Generator:
    """Same construction as run_targeted_panel.py's own _seeded_rng: a
    deterministic per-stratum draw from one shared seed via a hash, so the
    stratum draw order never depends on dict iteration order."""
    digest = hashlib.sha256(":".join([str(seed), *parts]).encode()).hexdigest()
    return np.random.default_rng(int(digest[:16], 16))


def build_matched_control_table(
    panel_genes: list[dict],
    luad_gene_stats: pd.DataFrame,
    *,
    strata: dict[str, tuple[str, ...]] = STRATA,
    detection_tolerance_log2: float = 0.5,
    rank_tolerance_percentile: float = 5.0,
    min_common_controls: int = MIN_COMMON_CONTROLS,
    seed: int = SEED,
) -> dict[str, dict]:
    """The six fixed 20-control matched strata tables, matched on LUAD-only
    detection fraction and median token rank (both computed by
    `median_token_rank_and_detection` / `load_and_freeze_luad_gene_stats`
    upstream of this function -- this function is pure pandas/numpy and
    takes the already-computed table as `luad_gene_stats`).

    `luad_gene_stats` must have columns `ensembl_id`, `detect_frac`,
    `median_token_rank`, covering every panel gene plus the full candidate
    universe (non-S100, non-ambient-anchor genes) -- one row per gene,
    LUAD-only, one fixed table shared by every stratum.

    A candidate qualifies for a stratum only if it is within
    `detection_tolerance_log2` log2 detection units AND
    `rank_tolerance_percentile` rank-percentile points of EVERY member of
    that stratum (not a stratum centroid) -- the design's own wording.
    Candidates that are themselves S100 panel genes or ambient-anchor genes
    are excluded from every stratum's candidate pool.

    A stratum with fewer than `min_common_controls` qualifying candidates
    is reported as not_estimable -- never widened, split, or given private
    controls (design doc's own rule, already enforced identically by
    ambient_stats.compute_Q_for_contrast's 20-count gate downstream).
    Otherwise exactly `min_common_controls` are drawn once, without
    replacement, from a generator seeded deterministically from `seed` and
    the stratum name -- never re-sampled.

    Returns {stratum_name: {"status": ..., "controls": tuple[str, ...] | None,
    "n_candidates": int}}. "controls" holds ensembl ids.
    """
    panel_ensembl_ids = {g["ensembl_id"] for g in panel_genes}
    anchor_ensembl_ids = {g["ensembl_id"] for g in panel_genes if g.get("role") == "ambient_anchor"}
    symbol_to_ensembl = {g["gene"]: g["ensembl_id"] for g in panel_genes}

    stats = luad_gene_stats.set_index("ensembl_id")
    estimable = stats[stats["detect_frac"] > 0].dropna(subset=["median_token_rank"])
    log2_detect = np.log2(estimable["detect_frac"])
    rank_pctile = rank_percentile_transform(estimable["median_token_rank"])

    candidate_pool = [
        eid for eid in estimable.index
        if eid not in panel_ensembl_ids and eid not in anchor_ensembl_ids
    ]

    result: dict[str, dict] = {}
    for stratum_name, member_symbols in strata.items():
        member_ids = [symbol_to_ensembl[s] for s in member_symbols]
        missing = [eid for eid in member_ids if eid not in estimable.index]
        if missing:
            result[stratum_name] = {
                "status": STATUS_NOT_ESTIMABLE_CONTROL_STRATUM, "controls": None,
                "n_candidates": 0,
                "reason": f"stratum member(s) {missing} not estimable in luad_gene_stats "
                          f"(zero detection or no token-positive cells in LUAD)",
            }
            continue

        member_log2 = log2_detect.loc[member_ids]
        member_pct = rank_pctile.loc[member_ids]

        qualifying = []
        for cand in candidate_pool:
            cand_log2 = log2_detect.loc[cand]
            cand_pct = rank_pctile.loc[cand]
            if (member_log2.sub(cand_log2).abs() <= detection_tolerance_log2).all() and \
               (member_pct.sub(cand_pct).abs() <= rank_tolerance_percentile).all():
                qualifying.append(cand)
        qualifying.sort()

        if len(qualifying) < min_common_controls:
            result[stratum_name] = {
                "status": STATUS_NOT_ESTIMABLE_CONTROL_STRATUM, "controls": None,
                "n_candidates": len(qualifying), "reason": None,
            }
            continue

        rng = _seeded_rng(seed, stratum_name)
        drawn = rng.choice(qualifying, size=min_common_controls, replace=False)
        result[stratum_name] = {
            "status": "eligible", "controls": tuple(sorted(drawn.tolist())),
            "n_candidates": len(qualifying), "reason": None,
        }
    return result


def synthetic_control_table(
    strata: dict[str, tuple[str, ...]] = STRATA,
    n_controls: int = MIN_COMMON_CONTROLS,
) -> dict[str, tuple[str, ...]]:
    """A structurally-valid FAKE control table for testing ambient_stats.py
    only. Fabricated ensembl-shaped ids ("ENSG9SYNTH...") that cannot
    collide with any real gene, so a caller can never mistake this for real
    controls by accident. Deterministic (no randomness needed -- the ids
    are just enumerated), since determinism of the *test* isn't the thing
    seed 20260922 is protecting.

    Still used by ambient_stats.py's own tests, unaffected by this module
    unblocking -- proving the E/Q/stratum-rollup/primary-test/stability
    code is correct never depended on the real controls existing.
    """
    table: dict[str, tuple[str, ...]] = {}
    for stratum_name in strata:
        table[stratum_name] = tuple(
            f"ENSG9SYNTH{stratum_name.upper()}{i:03d}" for i in range(n_controls)
        )
    return table
