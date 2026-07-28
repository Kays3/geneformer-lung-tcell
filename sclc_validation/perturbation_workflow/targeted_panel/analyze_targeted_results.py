#!/usr/bin/env python3
"""Merge delete + overexpress stats per comparison, flag concordant genes
(opposite-sign shift between delete and overexpress, both FDR<0.05), and
summarize results scoped to the 21-gene panel vs. the 29 top-driver genes.
"""
import json
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
STATS_DIR = HERE / "results"
PANEL_JSON = HERE / "target_gene_panel.json"
OUT_DIR = HERE / "results"

COMPARISONS = [
    ("sclc", "luad"), ("sclc", "normal"),
    ("luad", "sclc"), ("luad", "normal"),
    ("normal", "sclc"), ("normal", "luad"),
]
FDR_THRESH = 0.05


def main():
    gene_meta = {g["gene"]: g["source"] for g in json.loads(PANEL_JSON.read_text())["genes"]}

    all_rows = []
    for src, tgt in COMPARISONS:
        d = pd.read_csv(STATS_DIR / "delete" / f"targeted_delete_{src}_to_{tgt}.csv")
        o = pd.read_csv(STATS_DIR / "overexpress" / f"targeted_overexpress_{src}_to_{tgt}.csv")
        d = d[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections", "Sig"]].rename(
            columns={"Shift_to_goal_end": "delete_shift", "Goal_end_FDR": "delete_fdr", "N_Detections": "delete_n", "Sig": "delete_sig"}
        )
        o = o[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections", "Sig"]].rename(
            columns={"Shift_to_goal_end": "overexpress_shift", "Goal_end_FDR": "overexpress_fdr", "N_Detections": "overexpress_n", "Sig": "overexpress_sig"}
        )
        m = d.merge(o, on=["Gene_name", "Ensembl_ID"], how="outer")
        m["comparison"] = f"{src}_to_{tgt}"
        m["source_state"] = src
        m["goal_state"] = tgt
        m["gene_set"] = m["Gene_name"].map(gene_meta)
        all_rows.append(m)

    full = pd.concat(all_rows, ignore_index=True)

    both_sig = (full["delete_fdr"] < FDR_THRESH) & (full["overexpress_fdr"] < FDR_THRESH)
    opposite_sign = (full["delete_shift"] * full["overexpress_shift"]) < 0
    full["concordant"] = both_sig & opposite_sign

    full = full.sort_values(["comparison", "concordant"], ascending=[True, False])
    full.to_csv(OUT_DIR / "targeted_panel_delete_overexpress_merged.csv", index=False)

    concordant = full[full["concordant"]].copy()
    concordant["abs_delete_shift"] = concordant["delete_shift"].abs()
    concordant = concordant.sort_values("abs_delete_shift", ascending=False)
    concordant.to_csv(OUT_DIR / "targeted_panel_concordant_hits.csv", index=False)

    print(f"Total gene x comparison rows: {len(full)}")
    print(f"Concordant (both FDR<{FDR_THRESH}, opposite sign): {full['concordant'].sum()}")
    print()
    print("=== Concordant hits, panel genes only ===")
    panel_hits = concordant[concordant["gene_set"] == "panel"]
    print(panel_hits[["comparison", "Gene_name", "delete_shift", "overexpress_shift", "delete_fdr", "overexpress_fdr", "delete_n"]].to_string(index=False))
    print()
    print("=== Concordant hits, top-driver genes only ===")
    driver_hits = concordant[concordant["gene_set"] == "top_driver_luad_lusc_normal"]
    print(driver_hits[["comparison", "Gene_name", "delete_shift", "overexpress_shift", "delete_fdr", "overexpress_fdr", "delete_n"]].head(20).to_string(index=False))

    summary = full.groupby(["gene_set", "comparison"]).agg(
        n_genes=("Gene_name", "nunique"),
        n_concordant=("concordant", "sum"),
    ).reset_index()
    summary.to_csv(OUT_DIR / "targeted_panel_concordance_summary_by_comparison.csv", index=False)
    print()
    print("=== Concordance rate by gene set x comparison ===")
    print(summary.to_string(index=False))


if __name__ == "__main__":
    main()
