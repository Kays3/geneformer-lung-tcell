#!/usr/bin/env python3
"""Donor-level consistency check for the targeted-panel concordant hits.

The raw perturbation pickles store, per gene/source/type, a flat list of
per-cell shift values per possible goal state -- in the exact order
InSilicoPerturber processes cells: filtered to the target gene's token
(delete only; overexpress has no such filter) then sorted by token sequence
length descending (perturber_utils.downsample_and_sort always sorts,
regardless of max_ncells). Replicating that same filter+sort on the cached
per-source dataset recovers, in order, which donor each shift value
belongs to.
"""
from __future__ import annotations

import json
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from datasets import load_from_disk

HOME = Path.home()
RUN_DIR = HOME / "workspace/KD/sclc_luad_normal_htan_targeted_panel_perturbation"
RAW_ROOT = RUN_DIR / "raw"
SOURCE_DATA_DIR = RUN_DIR / "data" / "sources"
CONCORDANT_CSV = RUN_DIR / "results_import" / "targeted_panel_concordant_hits.csv"
TOKEN_DICT_FILE = HOME / "workspace/geneformer-uv-starter/geneformer-workspace/Geneformer/geneformer/token_dictionary_gc104M.pkl"
OUT_CSV = RUN_DIR / "results_import" / "targeted_panel_donor_consistency.csv"

STATE_NAMES = {
    "sclc": "small cell lung carcinoma",
    "luad": "lung adenocarcinoma",
    "normal": "normal",
}


def ordered_cells(source: str, gene_ensembl: str | None) -> pd.DataFrame:
    """Cells in the exact order InSilicoPerturber processed them for this
    source: (delete only) filtered to cells containing the gene token, then
    always sorted by length descending."""
    ds = load_from_disk(str(SOURCE_DATA_DIR / f"{source}.dataset"))
    if gene_ensembl is not None:
        with open(TOKEN_DICT_FILE, "rb") as f:
            token_dict = pickle.load(f)
        token = token_dict.get(gene_ensembl)
        ds = ds.filter(lambda row: token in row["input_ids"], num_proc=1)
    df = pd.DataFrame({"individual": ds["individual"], "length": ds["length"]})
    df = df.sort_values("length", ascending=False, kind="stable").reset_index(drop=True)
    return df


def load_shift_list(perturb_type: str, source: str, gene: str, goal_state_name: str) -> list[float] | None:
    raw_dir = RAW_ROOT / perturb_type / source
    matches = sorted(raw_dir.glob(f"in_silico_{perturb_type}_targeted_{source}_{gene}_cell_embs_dict_*_raw.pickle"))
    if not matches:
        return None
    with open(matches[0], "rb") as f:
        d = pickle.load(f)
    inner = d.get(goal_state_name)
    if not inner:
        return None
    # exactly one (token, 'cell_emb') key per single-gene run
    (key, values), = inner.items()
    return values


def donor_summary(source: str, gene: str, gene_ensembl: str, perturb_type: str, goal_state_slug: str) -> dict | None:
    goal_state_name = STATE_NAMES[goal_state_slug]
    shifts = load_shift_list(perturb_type, source, gene, goal_state_name)
    if shifts is None:
        return {"status": "no_pickle"}

    cells = ordered_cells(source, gene_ensembl if perturb_type == "delete" else None)
    if len(cells) != len(shifts):
        return {"status": f"length_mismatch cells={len(cells)} shifts={len(shifts)}"}

    cells = cells.copy()
    cells["shift"] = shifts
    by_donor = cells.groupby("individual")["shift"].agg(["mean", "count"]).reset_index()
    overall_mean = float(np.mean(shifts))
    overall_sign = np.sign(overall_mean)
    same_sign = (np.sign(by_donor["mean"]) == overall_sign) | (by_donor["mean"] == 0)

    return {
        "status": "ok",
        "n_cells": len(shifts),
        "n_donors": len(by_donor),
        "overall_mean_shift": overall_mean,
        "n_donors_same_sign": int(same_sign.sum()),
        "frac_donors_same_sign": float(same_sign.mean()),
        "donor_means_min": float(by_donor["mean"].min()),
        "donor_means_max": float(by_donor["mean"].max()),
        "donor_cell_counts": ";".join(f"{r.individual}:{r['count']}" for _, r in by_donor.iterrows()),
    }


def main() -> None:
    hits = pd.read_csv(CONCORDANT_CSV)
    rows = []
    for i, r in hits.iterrows():
        source = r["source_state"]
        gene = r["Gene_name"]
        ensembl = r["Ensembl_ID"]
        for ptype in ("delete", "overexpress"):
            res = donor_summary(source, gene, ensembl, ptype, r["goal_state"])
            row = {
                "comparison": r["comparison"],
                "gene": gene,
                "gene_set": r["gene_set"],
                "perturb_type": ptype,
            }
            row.update(res or {"status": "error"})
            rows.append(row)
            print(f"[{i+1}/{len(hits)}] {ptype}/{r['comparison']}/{gene}: {res.get('status')}", file=sys.stderr)

    out = pd.DataFrame(rows)
    out.to_csv(OUT_CSV, index=False)
    print(f"Wrote {OUT_CSV}: {len(out)} rows")


if __name__ == "__main__":
    main()
