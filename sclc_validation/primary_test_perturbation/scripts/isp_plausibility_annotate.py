#!/usr/bin/env python3
"""Annotate the top candidates of the complete in-silico perturbation (ISP) screen
with the evidence needed to judge biological plausibility.

This script MEASURES. It does not decide. Every column it writes is derived from a
committed table in this repository; the prose verdicts live in
`../reports/isp_plausibility_20260917.md` and are labelled there as judgment.

SELECTION RULE  (fixed before any result was looked at; see the report)
----------------------------------------------------------------------
Universe      every gene x comparison row of allgene_delete_overexpress_shift.csv
Gate          delete_fdr < 0.05  AND  overexpress_fdr < 0.05  AND  delete_n >= 25
              -- the significance and detection-adequacy halves of this repo's own
              `concordant` criterion, applied WITHOUT its sign requirement, so that
              same-sign rows stay in the table and can be counted rather than
              silently filtered out.
Rank          within each comparison x direction (direction = sign of delete_shift),
              by |delete_shift| descending. Effect size, not p-value: ranking by FDR
              ranks by detection count, which is the very confound under assessment.
              The FDR-based rank is carried alongside as `rank_in_stratum_by_fdr` so
              the divergence is visible instead of argued.
Top-N         N = 10 per comparison x direction  ->  up to 120 selected rows.
Always-in     all rows for S100A8 and S100A9, whether or not they pass the gate,
              flagged by `selected_by`. They were named in the work order.

WHY THE SIGN REQUIREMENT IS DROPPED FROM THE GATE, NOT FROM THE ANALYSIS
-----------------------------------------------------------------------
This repo defines a concordant hit as delete and overexpress shifts that are
significant, adequately detected, and OPPOSITE in sign -- overexpression pushing one
way and deletion undoing it. A same-sign pair means the model moved the same
direction whether the gene was added or removed, which is a response to perturbing
the token at all rather than to its dose. That is diagnostic, so it is recorded per
row (`same_sign_delete_overexpress`) rather than used to drop rows.

Usage:  python3 isp_plausibility_annotate.py [--top-n 10] [--out <csv>]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PTP = HERE.parent                                    # primary_test_perturbation/
SCLC = PTP.parent                                    # sclc_validation/

SHIFT_CSV = PTP / "tables" / "allgene_delete_overexpress_shift.csv"
AMBIENT_CSV = PTP / "tables" / "ambient_risk_all_genes.csv"
AMBIENT_MANIFEST = PTP / "tables" / "ambient_risk_manifest.json"
PAM_SUMMARY = SCLC / "immune_axis_test" / "results" / "axis_sensitivity_summary.json"
PAM_AUDIT = SCLC / "immune_axis_test" / "results" / "axis_sensitivity_exclusion_audit.csv"

FOCUS_GENES = ("S100A8", "S100A9")

# ---------------------------------------------------------------------------
# LINEAGE CALL -- THIS BLOCK IS JUDGMENT, NOT MEASUREMENT.
# Every other column in the output is computed from a committed table. This one
# is a literature-based assignment of the canonical cell type that expresses each
# gene, made by the analyst and cited in the report. It is here, in code, rather
# than in prose so that the counts derived from it are reproducible and so that a
# reader can disagree with a specific call by editing one line.
# Sources for the assignments: Travaglini et al. 2020 Nature 587:619-625 (lung
# epithelial and stromal identities); Edgeworth et al. 1991 J Biol Chem
# 266:7706-13 (S100A8/A9 myeloid restriction); standard marker usage for the
# neuroendocrine, mast, plasma-cell and myeloid entries.
# "t_cell_or_broad" means the gene is expressed by T cells or is broadly
# expressed across cell types -- NOT that it is T-cell specific.
LINEAGE_CALL = {
    # lung epithelial (alveolar type 2, club, secretory)
    "SFTPC": "lung_epithelial", "SFTPB": "lung_epithelial", "SFTPA1": "lung_epithelial",
    "SFTPA2": "lung_epithelial", "SCGB1A1": "lung_epithelial", "SCGB3A1": "lung_epithelial",
    "SCGB3A2": "lung_epithelial", "SLPI": "lung_epithelial", "TFF3": "lung_epithelial",
    "EPCAM": "lung_epithelial", "NAPSA": "lung_epithelial", "WFDC2": "lung_epithelial",
    "MUC1": "lung_epithelial", "AGER": "lung_epithelial", "CLDN18": "lung_epithelial",
    "KRT8": "lung_epithelial", "KRT18": "lung_epithelial", "KRT19": "lung_epithelial",
    "SMIM22": "lung_epithelial", "RBP1": "lung_epithelial", "MEST": "lung_epithelial",
    # myeloid
    "S100A8": "myeloid", "S100A9": "myeloid", "S100A12": "myeloid", "LYZ": "myeloid",
    "MARCO": "myeloid", "MSR1": "myeloid", "CD68": "myeloid", "FABP4": "myeloid",
    "CCL18": "myeloid", "MCEMP1": "myeloid", "APOC1": "myeloid", "CFD": "myeloid",
    "SERPINA1": "myeloid", "FCN1": "myeloid", "VCAN": "myeloid", "MNDA": "myeloid",
    # stromal / endothelial
    "COL1A1": "stromal", "COL1A2": "stromal", "COL3A1": "stromal", "DCN": "stromal",
    "LUM": "stromal", "ACTA2": "stromal", "PECAM1": "endothelial", "VWF": "endothelial",
    "MGP": "stromal",
    # neuroendocrine / SCLC tumour programme
    "CHGA": "neuroendocrine", "INSM1": "neuroendocrine", "GRP": "neuroendocrine",
    "UCHL1": "neuroendocrine", "INA": "neuroendocrine", "CALB2": "neuroendocrine",
    "BEX1": "neuroendocrine", "TAGLN3": "neuroendocrine", "CD24": "neuroendocrine",
    "MDK": "neuroendocrine", "TUBB2B": "neuroendocrine", "ASCL1": "neuroendocrine",
    # other non-T immune
    "JCHAIN": "plasma_cell", "TPSB2": "mast_cell",
    # erythroid
    "HBA1": "erythroid", "HBA2": "erythroid", "HBB": "erythroid", "HBD": "erythroid",
}


# Values asserted from the source files, so a changed input fails loudly here
# instead of quietly changing the report's numbers.
EXPECTED_ROWS = 82_702
EXPECTED_COMPARISONS = 6


def load_sources() -> tuple[pd.DataFrame, pd.DataFrame, dict, pd.DataFrame]:
    for path in (SHIFT_CSV, AMBIENT_CSV, AMBIENT_MANIFEST, PAM_SUMMARY, PAM_AUDIT):
        if not path.exists():
            raise SystemExit(f"FAIL: missing source {path}")

    shift = pd.read_csv(SHIFT_CSV)
    if len(shift) != EXPECTED_ROWS:
        raise SystemExit(f"FAIL: {SHIFT_CSV.name} has {len(shift)} rows, expected {EXPECTED_ROWS}")
    if shift["comparison"].nunique() != EXPECTED_COMPARISONS:
        raise SystemExit(f"FAIL: {shift['comparison'].nunique()} comparisons, expected {EXPECTED_COMPARISONS}")

    ambient = pd.read_csv(AMBIENT_CSV)
    manifest = json.loads(AMBIENT_MANIFEST.read_text())
    pam = json.loads(PAM_SUMMARY.read_text())
    audit = pd.read_csv(PAM_AUDIT)

    # The focus genes must exist in the screen, or every statement about them below
    # would be an assertion over an empty set rather than a measurement.
    for gene in FOCUS_GENES:
        n = int((shift["Gene_name"] == gene).sum())
        if n != EXPECTED_COMPARISONS:
            raise SystemExit(f"FAIL: {gene} has {n} rows in the screen, expected {EXPECTED_COMPARISONS}")
    return shift, ambient, manifest, pam, audit


def annotate(shift, ambient, manifest, pam, audit) -> pd.DataFrame:
    df = shift.copy()

    df["direction"] = df["delete_shift"].apply(lambda v: "toward_goal" if v > 0 else "away_from_goal")
    df["abs_delete_shift"] = df["delete_shift"].abs()
    df["same_sign_delete_overexpress"] = (df["delete_shift"] * df["overexpress_shift"]) > 0
    df["passes_gate"] = (
        (df["delete_fdr"] < 0.05) & (df["overexpress_fdr"] < 0.05) & (df["delete_n"] >= 25)
    )

    # --- ambient-RNA diagnostic (this repo's own, tables/ambient_risk_all_genes.csv)
    amb = ambient[["gene", "detect_frac", "mean_expr", "ambient_risk", "ambient_pct"]].rename(
        columns={"gene": "Gene_name"}
    )
    df = df.merge(amb, on="Gene_name", how="left")
    threshold = float(manifest["flag_threshold"])
    df["ambient_flag"] = df["ambient_risk"] > threshold
    df["ambient_scored"] = df["ambient_risk"].notna()

    anchors_amb = set(manifest["anchors_ambient"])
    anchors_tc = set(manifest["anchors_tcell"])
    df["is_ambient_anchor"] = df["Gene_name"].isin(anchors_amb)
    df["is_tcell_anchor"] = df["Gene_name"].isin(anchors_tc)

    # --- Pam's T1 axis sensitivity exclusions
    excluded = set(pam["excluded_genes"])
    in_panel_of_50 = set(audit["Gene_name"])
    df["pam_excluded"] = df["Gene_name"].isin(excluded)
    df["pam_assessed"] = df["Gene_name"].isin(in_panel_of_50)

    # judgment column, see LINEAGE_CALL
    df["lineage_call"] = df["Gene_name"].map(LINEAGE_CALL).fillna("t_cell_or_broad_or_unassigned")
    df["non_t_lineage"] = df["lineage_call"].isin(
        {"lung_epithelial", "myeloid", "stromal", "endothelial", "neuroendocrine",
         "plasma_cell", "mast_cell", "erythroid"}
    )

    # --- ranks WITHIN THE GATED SUBSET of each comparison x direction stratum.
    # Ranking over all rows instead would let a row ranked 62nd appear in a "top 10"
    # list, because selection happens after the gate: rank and selection have to be
    # computed over the same set or the table contradicts itself.
    gated = df[df["passes_gate"]]
    grp = gated.groupby(["comparison", "direction"])
    df["rank_gated_by_effect"] = grp["abs_delete_shift"].rank(ascending=False, method="min")
    df["rank_gated_by_fdr"] = grp["delete_fdr"].rank(ascending=True, method="min")
    df["n_gated_in_stratum"] = grp["abs_delete_shift"].transform("size")
    return df, threshold


def select(df: pd.DataFrame, top_n: int) -> pd.DataFrame:
    gated = df[df["passes_gate"]].copy()
    top = (
        gated.sort_values(["comparison", "direction", "abs_delete_shift"], ascending=[True, True, False])
        .groupby(["comparison", "direction"], as_index=False)
        .head(top_n)
    )
    top["selected_by"] = "top_n_by_effect"

    focus = df[df["Gene_name"].isin(FOCUS_GENES)].copy()
    focus["selected_by"] = "named_in_work_order"

    key = ["Gene_name", "comparison"]
    focus_only = focus.merge(top[key], on=key, how="left", indicator=True)
    focus_only = focus_only[focus_only["_merge"] == "left_only"].drop(columns="_merge")

    out = pd.concat([top, focus_only], ignore_index=True)
    if out.empty:
        raise SystemExit("FAIL: selection produced zero rows")
    if not set(FOCUS_GENES).issubset(set(out["Gene_name"])):
        raise SystemExit("FAIL: focus genes absent from the selected table")
    return out.sort_values(
        ["comparison", "direction", "abs_delete_shift"], ascending=[True, True, False]
    ).reset_index(drop=True)


COLUMNS = [
    "Gene_name", "comparison", "direction", "selected_by",
    "delete_shift", "delete_fdr", "delete_n",
    "overexpress_shift", "overexpress_fdr", "overexpress_n",
    "concordant", "same_sign_delete_overexpress", "passes_gate",
    "abs_delete_shift", "rank_gated_by_effect", "rank_gated_by_fdr", "n_gated_in_stratum",
    "detect_frac", "mean_expr", "ambient_risk", "ambient_pct", "ambient_flag",
    "ambient_scored", "is_ambient_anchor", "is_tcell_anchor",
    "pam_excluded", "pam_assessed", "lineage_call", "non_t_lineage",
]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--top-n", type=int, default=10)
    ap.add_argument("--out", type=Path, default=PTP / "tables" / "isp_plausibility_top_candidates.csv")
    args = ap.parse_args()

    shift, ambient, manifest, pam, audit = load_sources()
    df, threshold = annotate(shift, ambient, manifest, pam, audit)
    out = select(df, args.top_n)
    out[COLUMNS].to_csv(args.out, index=False)

    # ---- everything the report quotes is printed here, so the prose can be checked
    # against this run rather than against a remembered number.
    print(f"SELECTION RULE: gate delete_fdr<0.05 & overexpress_fdr<0.05 & delete_n>=25;")
    print(f"                top {args.top_n} per comparison x direction by |delete_shift|;")
    print(f"                plus all rows for {', '.join(FOCUS_GENES)}.")
    print(f"ambient flag threshold (from manifest): {threshold:.6f}")
    print()
    print(f"screen rows                     : {len(df)}")
    print(f"rows passing the gate           : {int(df['passes_gate'].sum())}")
    print(f"  of those, same-sign           : {int((df['passes_gate'] & df['same_sign_delete_overexpress']).sum())}")
    print(f"  of those, concordant (repo)   : {int((df['passes_gate'] & (df['concordant'] == True)).sum())}")
    print(f"concordant rows in whole screen : {int((df['concordant'] == True).sum())}")
    print(f"selected rows written           : {len(out)}  -> {args.out}")
    print(f"  top_n_by_effect               : {int((out['selected_by'] == 'top_n_by_effect').sum())}")
    print(f"  named_in_work_order (extra)   : {int((out['selected_by'] == 'named_in_work_order').sum())}")
    print()
    sel = out[out["selected_by"] == "top_n_by_effect"]
    print(f"selected top-N: ambient-flagged  : {int(sel['ambient_flag'].sum())} / {len(sel)}")
    print(f"selected top-N: ambient anchors  : {int(sel['is_ambient_anchor'].sum())}")
    print(f"selected top-N: T-cell anchors   : {int(sel['is_tcell_anchor'].sum())}")
    print(f"selected top-N: Pam-excluded     : {int(sel['pam_excluded'].sum())}")
    print(f"selected top-N: same-sign        : {int(sel['same_sign_delete_overexpress'].sum())} / {len(sel)}")
    print(f"selected top-N: distinct genes   : {sel['Gene_name'].nunique()}")
    print(f"selected top-N: median delete_n  : {sel['delete_n'].median():.0f}")
    print()
    print("selected top-N by LINEAGE CALL (judgment column, see LINEAGE_CALL):")
    for lin, n in sel["lineage_call"].value_counts().items():
        genes = sorted(sel.loc[sel["lineage_call"] == lin, "Gene_name"].unique())
        print(f"  {lin:32} rows={n:3d} genes={len(genes):2d}  {', '.join(genes)}")
    print(f"  -> rows on a NON-T lineage      : {int(sel['non_t_lineage'].sum())} / {len(sel)}")
    print()
    print("recurring genes (top-N in >=2 comparison x direction strata):")
    rec = (sel.groupby("Gene_name")
              .agg(strata=("comparison", "size"), concordant_rows=("concordant", "sum"),
                   same_sign_rows=("same_sign_delete_overexpress", "sum"),
                   max_abs_shift=("abs_delete_shift", "max"), max_n=("delete_n", "max"),
                   detect_frac=("detect_frac", "max"), ambient_risk=("ambient_risk", "max"),
                   ambient_flag=("ambient_flag", "max"), lineage=("lineage_call", "first"))
              .query("strata >= 2").sort_values(["strata", "max_abs_shift"], ascending=False))
    with pd.option_context("display.width", 220, "display.max_rows", 80):
        print(rec.to_string())
    print()
    print("--- focus genes, every comparison (all columns that the report quotes)")
    cols = ["comparison", "Gene_name", "delete_shift", "delete_fdr", "delete_n",
            "overexpress_shift", "overexpress_fdr", "concordant",
            "same_sign_delete_overexpress", "passes_gate",
            "rank_gated_by_effect", "rank_gated_by_fdr", "n_gated_in_stratum"]
    focus = df[df["Gene_name"].isin(FOCUS_GENES)].sort_values(["Gene_name", "comparison"])
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print(focus[cols].to_string(index=False))
    print()
    print("--- focus genes, ambient diagnostic (one row per gene)")
    amb_cols = ["Gene_name", "detect_frac", "mean_expr", "ambient_risk", "ambient_pct",
                "ambient_flag", "is_ambient_anchor", "is_tcell_anchor", "pam_excluded"]
    print(focus[amb_cols].drop_duplicates("Gene_name").to_string(index=False))
    print()
    print("--- S100 family context: which S100A genes look ambient and which do not")
    s100 = df[df["Gene_name"].str.match(r"^S100A", na=False)]
    s100 = s100[amb_cols].drop_duplicates("Gene_name").sort_values("ambient_risk", ascending=False)
    print(s100.to_string(index=False))


if __name__ == "__main__":
    main()
