#!/usr/bin/env python3
"""Donor-level consistency for the all-gene held-out perturbation results.

This closes the gap left open by the primary all-gene report: every concordant hit
there is a cell-level statistic, so a candidate driven by one patient is
indistinguishable from one that reproduces across patients.

Reconstruction
--------------
The all-gene run stores raw output per cell-batch, not per gene as the targeted panel
did. Each file is named `..._shard{S}_dict_cell_embs_{C}batch{B}_raw.pickle` and holds
`{goal_state: {(gene_token, 'cell_emb'): [one shift value]}}` for a single cell.

`C` indexes cells in **length-descending order within the shard**, because
`perturber_utils.downsample_and_sort` always sorts. This is verified empirically at
startup: the number of perturbed genes for cell `C` must equal `length - 2` for the
`C`-th longest cell in that shard. Combined with `heldout_shard_manifest.csv`, which
carries `individual` per cell, that recovers the donor behind every shift value.

Validation
----------
Before any consistency claim, the per-cell values are re-aggregated per gene and
compared against `Shift_to_goal_end` in the Geneformer stats CSVs. If the
reconstruction were misaligned, those means would not agree.
"""
from __future__ import annotations

import os
import pickle
import re
import sys
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

HOME = Path.home()
ROOT = Path(os.environ.get("SCLC_PERTURBATION_ROOT", HOME / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"))
TOKEN_DICT = Path(os.environ.get(
    "GENEFORMER_TOKEN_DICT",
    HOME / "workspace/geneformer-uv-starter/geneformer-workspace/Geneformer/geneformer/token_dictionary_gc104M.pkl",
))
OUT_DIR = Path(os.environ.get("DONOR_OUT_DIR", Path(__file__).resolve().parents[1] / "tables"))
HITS_CSV = Path(os.environ.get("PRIMARY_HITS", Path(__file__).resolve().parents[1] / "tables" / "primary_concordant_hits.csv"))

SOURCES = ["normal", "sclc", "luad"]
ARMS = ["delete", "overexpress"]
STATE_NAMES = {"sclc": "small cell lung carcinoma", "luad": "lung adenocarcinoma", "normal": "normal"}
NAME_TO_SLUG = {v: k for k, v in STATE_NAMES.items()}
FILE_RE = re.compile(r"shard(\d+)_dict_cell_embs_(\d+)batch(-?\d+)_raw\.pickle$")
N_WORKERS = int(os.environ.get("N_WORKERS", 8))


def shard_cell_donors(manifest: pd.DataFrame) -> dict[tuple[str, int], list[str]]:
    """Map (source, shard) -> donor per cell index, in length-descending order."""
    mapping = {}
    for (source, shard), group in manifest.groupby(["source", "shard"]):
        ordered = group.sort_values("length", ascending=False, kind="stable")
        mapping[(source, int(shard))] = list(ordered["individual"])
    return mapping


def shard_cell_lengths(manifest: pd.DataFrame) -> dict[tuple[str, int], list[int]]:
    mapping = {}
    for (source, shard), group in manifest.groupby(["source", "shard"]):
        ordered = group.sort_values("length", ascending=False, kind="stable")
        mapping[(source, int(shard))] = list(ordered["length"])
    return mapping


def verify_ordering(arm: str, source: str, donors: dict, lengths: dict, n_shards: int = 3) -> None:
    """Assert that cell index C corresponds to the C-th longest cell in the shard."""
    raw_dir = ROOT / "raw" / arm / source
    checked = 0
    for shard in range(n_shards):
        per_cell: dict[int, set] = defaultdict(set)
        for path in raw_dir.glob(f"*shard{shard:04d}_dict_cell_embs_*_raw.pickle"):
            match = FILE_RE.search(path.name)
            if not match:
                continue
            cell_index = int(match.group(2))
            with open(path, "rb") as handle:
                payload = pickle.load(handle)
            any_state = next(iter(payload.values()))
            per_cell[cell_index].update(token for token, _ in any_state)
        if not per_cell:
            continue
        expected = lengths[(source, shard)]
        for cell_index, tokens in sorted(per_cell.items()):
            if len(tokens) != expected[cell_index] - 2:
                raise AssertionError(
                    f"ordering mismatch {arm}/{source} shard{shard} cell{cell_index}: "
                    f"{len(tokens)} genes vs expected {expected[cell_index] - 2}"
                )
        checked += 1
    print(f"  [verify] {arm}/{source}: cell ordering confirmed on {checked} shards", flush=True)


def accumulate_shard(args: tuple[str, str, int, list[str]]) -> dict:
    """Sum shift values per (token, goal_slug, donor) for one shard."""
    arm, source, shard, donors = args
    raw_dir = ROOT / "raw" / arm / source
    totals: dict[tuple, list[float]] = defaultdict(lambda: [0.0, 0])
    for path in raw_dir.glob(f"*shard{shard:04d}_dict_cell_embs_*_raw.pickle"):
        match = FILE_RE.search(path.name)
        if not match:
            continue
        cell_index = int(match.group(2))
        if cell_index >= len(donors):
            continue
        donor = donors[cell_index]
        with open(path, "rb") as handle:
            payload = pickle.load(handle)
        for goal_name, inner in payload.items():
            goal = NAME_TO_SLUG.get(goal_name)
            if goal is None or goal == source:
                continue
            for (token, _), values in inner.items():
                entry = totals[(token, goal, donor)]
                entry[0] += float(np.sum(values))
                entry[1] += len(values)
    return dict(totals)


def collect(arm: str, source: str, donors_by_shard: dict) -> pd.DataFrame:
    shards = sorted(s for (src, s) in donors_by_shard if src == source)
    jobs = [(arm, source, shard, donors_by_shard[(source, shard)]) for shard in shards]
    merged: dict[tuple, list[float]] = defaultdict(lambda: [0.0, 0])
    done = 0
    with ProcessPoolExecutor(max_workers=N_WORKERS) as pool:
        for partial in pool.map(accumulate_shard, jobs, chunksize=4):
            for key, (total, count) in partial.items():
                entry = merged[key]
                entry[0] += total
                entry[1] += count
            done += 1
            if done % 50 == 0:
                print(f"    {arm}/{source}: {done}/{len(jobs)} shards", flush=True)
    rows = [
        {"token": token, "goal_state": goal, "individual": donor, "sum_shift": total, "n_cells": count}
        for (token, goal, donor), (total, count) in merged.items()
    ]
    frame = pd.DataFrame(rows)
    frame["source_state"] = source
    frame["arm"] = arm
    frame["donor_mean_shift"] = frame.sum_shift / frame.n_cells
    return frame


def validate_against_stats(per_donor: pd.DataFrame, token_to_ensembl: dict) -> pd.DataFrame:
    """Re-aggregate to gene level and compare with the Geneformer stats CSVs."""
    overall = (
        per_donor.groupby(["arm", "source_state", "goal_state", "token"], as_index=False)
        .agg(sum_shift=("sum_shift", "sum"), n_cells=("n_cells", "sum"))
    )
    overall["reconstructed_shift"] = overall.sum_shift / overall.n_cells
    overall["Ensembl_ID"] = overall.token.map(token_to_ensembl)

    checks = []
    for (arm, source, goal), group in overall.groupby(["arm", "source_state", "goal_state"]):
        stats_path = ROOT / "stats" / arm / f"heldout_allgene_{arm}_{source}_to_{goal}.csv"
        if not stats_path.exists():
            continue
        stats = pd.read_csv(stats_path)[["Ensembl_ID", "Shift_to_goal_end", "N_Detections"]]
        merged = group.merge(stats, on="Ensembl_ID", how="inner")
        if merged.empty:
            continue
        delta = (merged.reconstructed_shift - merged.Shift_to_goal_end).abs()
        checks.append({
            "arm": arm, "comparison": f"{source}_to_{goal}", "n_genes_compared": len(merged),
            "max_abs_diff": float(delta.max()), "median_abs_diff": float(delta.median()),
            "n_detection_mismatch": int((merged.n_cells != merged.N_Detections).sum()),
            "correlation": float(np.corrcoef(merged.reconstructed_shift, merged.Shift_to_goal_end)[0, 1]),
        })
    return pd.DataFrame(checks)


def classify(frac_delete: float, frac_over: float, donors_delete: int, donors_over: int) -> str:
    if min(donors_delete, donors_over) < 2:
        return "single_donor_only"
    worst = min(frac_delete, frac_over)
    if worst >= 1.0:
        return "fully_consistent"
    if worst >= 0.5:
        return "majority_consistent"
    return "inconsistent"


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(ROOT / "tables" / "heldout_shard_manifest.csv")
    donors_by_shard = shard_cell_donors(manifest)
    lengths_by_shard = shard_cell_lengths(manifest)

    print("Verifying cell-index ordering assumption")
    for arm in ARMS:
        for source in SOURCES:
            verify_ordering(arm, source, donors_by_shard, lengths_by_shard)

    with open(TOKEN_DICT, "rb") as handle:
        token_dict = pickle.load(handle)
    token_to_ensembl = {v: k for k, v in token_dict.items()}

    print("\nAccumulating per-donor shifts")
    frames = []
    for arm in ARMS:
        for source in SOURCES:
            print(f"  {arm}/{source}", flush=True)
            frames.append(collect(arm, source, donors_by_shard))
    per_donor = pd.concat(frames, ignore_index=True)
    per_donor["Ensembl_ID"] = per_donor.token.map(token_to_ensembl)
    per_donor["comparison"] = per_donor.source_state + "_to_" + per_donor.goal_state
    per_donor.to_parquet(OUT_DIR / "allgene_per_donor_shifts.parquet", index=False)
    print(f"  per-donor rows: {len(per_donor):,}")

    print("\nValidating reconstruction against Geneformer stats")
    validation = validate_against_stats(per_donor, token_to_ensembl)
    validation.to_csv(OUT_DIR / "allgene_donor_reconstruction_validation.csv", index=False)
    print(validation.to_string(index=False))
    worst = validation.max_abs_diff.max() if not validation.empty else float("nan")
    if not np.isfinite(worst) or worst > 1e-6:
        print(f"\n  [WARN] reconstruction differs from stats by up to {worst:.3g}", file=sys.stderr)
    else:
        print(f"\n  reconstruction matches stats exactly (max abs diff {worst:.3g})")

    print("\nScoring donor consistency for concordant hits")
    hits = pd.read_csv(HITS_CSV)
    hits["source_state"] = hits.comparison.str.split("_to_").str[0]

    per_donor["sign"] = np.sign(per_donor.donor_mean_shift)
    gene_overall = (
        per_donor.groupby(["arm", "comparison", "Ensembl_ID"], as_index=False)
        .agg(sum_shift=("sum_shift", "sum"), n_cells=("n_cells", "sum"))
    )
    gene_overall["overall_sign"] = np.sign(gene_overall.sum_shift / gene_overall.n_cells)

    joined = per_donor.merge(
        gene_overall[["arm", "comparison", "Ensembl_ID", "overall_sign"]],
        on=["arm", "comparison", "Ensembl_ID"], how="left",
    )
    joined["agrees"] = (joined.sign == joined.overall_sign) | (joined.donor_mean_shift == 0)
    agreement = (
        joined.groupby(["arm", "comparison", "Ensembl_ID"], as_index=False)
        .agg(n_donors=("individual", "nunique"), frac_same_sign=("agrees", "mean"))
    )

    wide = agreement.pivot_table(index=["comparison", "Ensembl_ID"], columns="arm",
                                 values=["n_donors", "frac_same_sign"])
    wide.columns = [f"{a}_{b}" for a, b in wide.columns]
    wide = wide.reset_index()

    result = hits.merge(wide, on=["comparison", "Ensembl_ID"], how="left")
    for column in ("n_donors_delete", "n_donors_overexpress", "frac_same_sign_delete", "frac_same_sign_overexpress"):
        if column not in result:
            result[column] = np.nan
    result["donor_robustness"] = [
        classify(row.frac_same_sign_delete, row.frac_same_sign_overexpress,
                 row.n_donors_delete if pd.notna(row.n_donors_delete) else 0,
                 row.n_donors_overexpress if pd.notna(row.n_donors_overexpress) else 0)
        if pd.notna(row.frac_same_sign_delete) and pd.notna(row.frac_same_sign_overexpress)
        else "not_recovered"
        for row in result.itertuples()
    ]
    result.to_csv(OUT_DIR / "allgene_concordant_hits_with_donor_robustness.csv", index=False)

    summary = result.donor_robustness.value_counts().rename_axis("donor_robustness").reset_index(name="n_hits")
    summary["pct"] = (100 * summary.n_hits / len(result)).round(1)
    summary.to_csv(OUT_DIR / "allgene_donor_robustness_summary.csv", index=False)
    print(summary.to_string(index=False))

    by_comparison = (
        result.groupby(["comparison", "donor_robustness"]).size().unstack(fill_value=0)
    )
    by_comparison.to_csv(OUT_DIR / "allgene_donor_robustness_by_comparison.csv")
    print()
    print(by_comparison.to_string())
    print(f"\nWrote outputs to {OUT_DIR}")


if __name__ == "__main__":
    main()
