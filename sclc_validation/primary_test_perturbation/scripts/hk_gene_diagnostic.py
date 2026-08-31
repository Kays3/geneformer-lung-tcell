#!/usr/bin/env python3
"""Housekeeping/ubiquitous-gene enrichment diagnostic for T-cell ISP hits.

Motivating question (external review comment): a substantial fraction of genes
ranked as top in-silico-perturbation (ISP) hits are housekeeping (HK) or
broadly-expressed genes. Two explanations are on the table -- (a) perturbing an
essential/ubiquitous gene can genuinely produce a large transcriptomic shift,
or (b) this is a systematic property of Geneformer's rank-value encoding and
ISP procedure, not disease-specific T-cell biology. This script does not
resolve that question; it quantifies how big the HK footprint actually is in
this project's own hit lists, and applies the one check the Geneformer authors
recommend before treating a candidate HK gene as informative: is it
differentially detected between the source and goal T-cell states.

Two gene-level flags, computed independently and reported separately because
they are different mechanisms:

  1. Ambient / lineage-foreign (already diagnosed elsewhere in this repo, see
     ../METHODS_ambient_risk.md): a gene a T cell should not transcribe at
     all -- surfactant, haemoglobin, myeloid, epithelial, stromal markers.
     Reused here by joining `tables/ambient_risk_all_genes.csv` (11,047 genes,
     calibrated logistic score, cross-validated AUC = 1.0 on held-out anchors)
     against this project's hit lists. Not recomputed.

  2. Housekeeping/ubiquitous (new in this script): genes genuinely transcribed
     by T cells but constitutively, at similar levels regardless of cell
     state -- ribosomal proteins, translation factors, heat-shock proteins,
     proteasome subunits, cytoskeleton, and a short list of classic reference
     genes. This is a coarse, standard family list (the kind used for
     "percent.ribo"/"percent.hsp" QC in single-cell pipelines), not an
     exhaustive HK-gene catalogue, and that limitation is reported plainly.

For every hit list -- targeted 50-gene panel and whole-genome, both the
delete-vs-overexpress concordant hits and the goal-vs-alt significant movers --
this script reports the HK/ambient fraction among hits against the background
rate among all tested genes (Fisher's exact test), and, for whole-genome HK
hits, the per-gene detection-fraction gap between the source and goal T-cell
state (delete_n / total held-out cells for that state) as a coarse stand-in
for differential expression -- coarse because detection rate is not the same
statistic as measured fold-change, which is why this is a proxy, not the T2
analysis the immune-axis PLAN still has scheduled.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd
from scipy.stats import fisher_exact

ROOT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parents[1]
TARGETED_DIR = ROOT / "sclc_validation/perturbation_workflow/targeted_panel"
TARGETED_MERGED = TARGETED_DIR / "results/targeted_panel_delete_overexpress_merged.csv"
TARGETED_PANEL_JSON = TARGETED_DIR / "target_gene_panel.json"
ALLGENE_DVO = HERE / "tables/allgene_delete_overexpress_shift.csv"
ALLGENE_GOAL_ALT = HERE / "tables/allgene_goal_vs_alt_shift.csv"
AMBIENT_SCORES = HERE / "tables/ambient_risk_all_genes.csv"
AMBIENT_MANIFEST = HERE / "tables/ambient_risk_manifest.json"
OUT = HERE / "tables/hk_gene_diagnostic"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]

# Authoritative held-out T-cell counts per source state, from
# cohort_and_workload_summary.csv (all-gene perturbation run, thinkstation2).
# overexpress_n is NOT a usable stand-in for this at genome scale: it tracks
# baseline detection, same as delete_n, not a census-complete perturbed count
# (that property held only for the targeted 50-gene panel's curated genes).
TOTAL_CELLS = {"luad": 6387, "sclc": 2424, "normal": 566}

# Classic, low-copy-number reference genes used across the literature as
# qPCR/normalization controls -- constitutively expressed, not disease-linked.
HK_REFERENCE = {
    "ACTB", "GAPDH", "B2M", "TBP", "PPIA", "YWHAZ", "HPRT1", "SDHA", "PGK1",
    "UBC", "RPLP0", "EEF1A1", "EEF2", "POLR2A", "TUBB", "PSMB4",
}

# Broad constitutive gene families (regex on symbol prefix). Coarse by
# construction -- a handful of ribosomal or heat-shock paralogs are known to
# be regulated, but as a family these are the standard single-cell QC classes.
HK_FAMILY_PATTERNS = [
    r"^RPL\d", r"^RPS\d", r"^MRPL\d", r"^MRPS\d",       # ribosomal proteins
    r"^EEF\d", r"^EIF\d",                                # translation factors
    r"^HSPA\d", r"^HSPB\d", r"^HSP90", r"^DNAJ", r"^HSPD1$", r"^HSPE1$",  # heat shock
    r"^PSMA\d", r"^PSMB\d", r"^PSMC\d", r"^PSMD\d",      # proteasome
    r"^ACTB$", r"^ACTG1$", r"^TUBA\d", r"^TUBB\d",       # cytoskeleton
    r"^NDUFA\d", r"^NDUFB\d", r"^COX4", r"^COX5", r"^ATP5",  # oxidative phosphorylation
]
_HK_FAMILY_RE = re.compile("|".join(HK_FAMILY_PATTERNS))


def is_hk(gene: str) -> bool:
    return gene in HK_REFERENCE or bool(_HK_FAMILY_RE.match(gene))


def load_ambient() -> pd.DataFrame:
    df = pd.read_csv(AMBIENT_SCORES)
    threshold = json.loads(AMBIENT_MANIFEST.read_text())["flag_threshold"]
    df["ambient_flag"] = df["ambient_risk"] >= threshold
    return df[["gene", "ambient_risk", "ambient_flag"]].rename(columns={"gene": "Gene_name"})


def enrichment_table(df: pd.DataFrame, hit_col: str, group_col: str = "comparison") -> pd.DataFrame:
    """Fisher's exact test: HK fraction among hits vs among all tested genes, per group."""
    rows = []
    for group, sub in df.groupby(group_col):
        hk = sub["is_hk"]
        hit = sub[hit_col].astype(bool)
        table = [
            [int((hk & hit).sum()), int((~hk & hit).sum())],
            [int((hk & ~hit).sum()), int((~hk & ~hit).sum())],
        ]
        odds_ratio, p_value = fisher_exact(table)
        n_hit = int(hit.sum())
        rows.append({
            group_col: group,
            "n_genes": len(sub),
            "n_hit": n_hit,
            "n_hk_total": int(hk.sum()),
            "n_hk_among_hit": int((hk & hit).sum()),
            "pct_hk_background": round(100 * hk.mean(), 2),
            "pct_hk_among_hit": round(100 * (hk & hit).sum() / n_hit, 2) if n_hit else float("nan"),
            "odds_ratio": round(odds_ratio, 3),
            "fisher_p": p_value,
        })
    return pd.DataFrame(rows)


def targeted_panel_analysis(ambient: pd.DataFrame) -> pd.DataFrame:
    panel_meta = json.loads(TARGETED_PANEL_JSON.read_text())
    gene_source = {g["gene"]: g["source"] for g in panel_meta["genes"]}
    df = pd.read_csv(TARGETED_MERGED)
    df["is_hk"] = df["Gene_name"].map(is_hk)
    df["gene_source"] = df["Gene_name"].map(gene_source)
    df = df.merge(ambient, on="Gene_name", how="left")
    df["ambient_flag"] = df["ambient_flag"].fillna(False)

    per_gene = df.drop_duplicates("Gene_name")[["Gene_name", "gene_source", "is_hk", "ambient_flag"]]
    per_gene.to_csv(OUT / "targeted_panel_gene_flags.csv", index=False)

    enrichment = enrichment_table(df, "concordant")
    enrichment.to_csv(OUT / "targeted_panel_hk_enrichment.csv", index=False)
    return df


def allgene_analysis(ambient: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    dvo = pd.read_csv(ALLGENE_DVO)
    dvo["is_hk"] = dvo["Gene_name"].map(is_hk)
    dvo = dvo.merge(ambient, on="Gene_name", how="left")
    dvo["ambient_flag"] = dvo["ambient_flag"].fillna(False)
    dvo_enrichment = enrichment_table(dvo, "concordant")
    dvo_enrichment.to_csv(OUT / "allgene_dvo_hk_enrichment.csv", index=False)

    ga = pd.read_csv(ALLGENE_GOAL_ALT)
    ga = ga.rename(columns={"Gene_name": "Gene_name"})
    ga["is_hk"] = ga["Gene_name"].map(is_hk)
    ga = ga.merge(ambient, on="Gene_name", how="left")
    ga["ambient_flag"] = ga["ambient_flag"].fillna(False)
    ga_enrichment_rows = []
    for arm in ga["arm"].unique():
        sub = ga[ga["arm"] == arm]
        e = enrichment_table(sub, "sig")
        e.insert(0, "arm", arm)
        ga_enrichment_rows.append(e)
    ga_enrichment = pd.concat(ga_enrichment_rows, ignore_index=True)
    ga_enrichment.to_csv(OUT / "allgene_goal_alt_hk_enrichment.csv", index=False)

    return dvo, ga


def detection_fraction_de_proxy(dvo: pd.DataFrame) -> pd.DataFrame:
    """For each gene, its delete-arm detection fraction as a source in each of
    the three states (delete_n / total held-out cells for that state) -- a
    coarse proxy for 'is this gene's transcript prevalence different between
    the source and goal T-cell state', per the Geneformer authors' recommended
    check before trusting an HK gene as a real perturbation-response driver.
    """
    dvo = dvo.copy()
    dvo["source_state"] = dvo["comparison"].str.split("_to_").str[0]
    dvo["goal_state"] = dvo["comparison"].str.split("_to_").str[1]
    dvo["source_detect_frac"] = dvo["delete_n"] / dvo["source_state"].map(TOTAL_CELLS)

    # detection fraction of the same gene when the GOAL state is itself the
    # source of a different comparison row. Two rows share each (gene,
    # source_state) pair (one per alternate goal) with an identical
    # source_detect_frac, so drop duplicates before the join or it fans out.
    goal_frac = dvo[["Gene_name", "source_state", "source_detect_frac"]].drop_duplicates().rename(
        columns={"source_state": "goal_state", "source_detect_frac": "goal_detect_frac"})
    merged = dvo.merge(goal_frac, on=["Gene_name", "goal_state"], how="left")
    merged["detect_frac_gap"] = (merged["goal_detect_frac"] - merged["source_detect_frac"]).abs()

    hk_concordant = merged[merged["is_hk"] & merged["concordant"]].copy()
    hk_concordant = hk_concordant.sort_values("detect_frac_gap")
    cols = ["Gene_name", "comparison", "source_state", "goal_state", "source_detect_frac",
            "goal_detect_frac", "detect_frac_gap", "delete_shift", "overexpress_shift",
            "ambient_flag"]
    hk_concordant[cols].to_csv(OUT / "allgene_hk_concordant_detection_gap.csv", index=False)
    return hk_concordant[cols]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    ambient = load_ambient()

    targeted = targeted_panel_analysis(ambient)
    dvo, ga = allgene_analysis(ambient)
    gap = detection_fraction_de_proxy(dvo)

    n_panel_hk = targeted.drop_duplicates("Gene_name")["is_hk"].sum()
    print(f"Targeted panel: {n_panel_hk}/50 genes flagged HK/ubiquitous by family or reference list")
    print(f"  of which top_driver-sourced: "
          f"{targeted.drop_duplicates('Gene_name').query('is_hk and gene_source == \"top_driver_luad_lusc_normal\"').shape[0]}"
          f" / {targeted.drop_duplicates('Gene_name').query('gene_source == \"top_driver_luad_lusc_normal\"').shape[0]}")

    print("\nWhole-genome delete-vs-overexpress: HK enrichment among concordant hits")
    print(pd.read_csv(OUT / "allgene_dvo_hk_enrichment.csv")[
        ["comparison", "n_hit", "pct_hk_background", "pct_hk_among_hit", "odds_ratio", "fisher_p"]
    ].to_string(index=False))

    print(f"\nHK-flagged, concordant, near-zero detection gap (<0.02) genes: "
          f"{(gap['detect_frac_gap'] < 0.02).sum()} / {len(gap)}")
    print(f"Wrote tables to {OUT.relative_to(ROOT)}/")


if __name__ == "__main__":
    main()
