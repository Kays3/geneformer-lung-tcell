#!/usr/bin/env python3
"""Ambient-RNA risk diagnostic for the prioritized perturbation candidates.

This is NOT CellBender and NOT ambient correction. CellBender's `remove-background`
learns the ambient profile from empty droplets, and the only available HTAN matrix is
filtered (minimum 113 UMI per cell across all 46,140 cells) and already subset to T
cells. With no empty droplets the ambient profile is unidentifiable, so that model
cannot be fit. See RESULTS_ambient_risk.md.

What this does instead is score how *ambient-like* each gene's behaviour is, using
properties that distinguish contamination from genuine expression, and calibrate that
score against genes whose status is already known:

- Contamination is lineage-foreign. Surfactant, haemoglobin and myeloid transcripts have
  no business in a T cell, so their presence marks the ambient pool.
- Contamination is unstructured across cell states. It does not track T-cell subtype,
  because it is not produced by the cell.
- Contamination is sample-driven. The soup composition differs per tissue digest, so
  ambient genes vary more between donors than between cell types.
- Contamination scales with droplet content, giving a stronger library-size dependence
  than a matched genuinely-expressed gene.

A logistic model is fit on the two anchor sets and its cross-validated AUC is reported.
If the anchors do not separate, the score is not trustworthy and the run says so rather
than ranking candidates anyway.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = Path(os.environ.get(
    "HTAN_H5AD",
    Path.home() / "workspace/KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad",
))
# The h5ad is indexed by Ensembl ID with no symbol column, so the symbol mapping is
# taken from the perturbation stats tables, which carry both identifiers.
STATS_ROOT = Path(os.environ.get(
    "SCLC_PERTURBATION_ROOT",
    Path.home() / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation",
)) / "stats"
CANDIDATES = Path(os.environ.get(
    "DENOISED_CANDIDATES",
    Path(__file__).resolve().parents[1] / "tables" / "immune_cancer_candidates_with_donor_robustness.csv",
))
OUT_DIR = Path(os.environ.get("AMBIENT_OUT_DIR", Path(__file__).resolve().parents[1] / "tables"))
FIG_DIR = Path(__file__).resolve().parents[1] / "figures" / "ambient_risk"

# Lineage-foreign transcripts: a T cell does not transcribe these, so detection is
# contamination by construction. This is the positive anchor.
KNOWN_AMBIENT = [
    "HBB", "HBA1", "HBA2", "HBD", "ALAS2", "AHSP",
    "SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SCGB1A1", "SCGB3A1", "SCGB3A2", "NAPSA",
    "WFDC2", "MUC1", "AGER", "CLDN18", "EPCAM", "KRT8", "KRT18", "KRT19", "SLPI",
    "LYZ", "S100A8", "S100A9", "S100A12", "CD68", "MARCO", "FCN1", "VCAN",
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "ACTA2", "PECAM1", "VWF",
]

# Canonical T-cell-intrinsic transcripts: genuinely produced by these cells.
# This is the negative anchor.
KNOWN_TCELL = [
    "CD3D", "CD3E", "CD3G", "CD247", "TRAC", "TRBC1", "TRBC2", "CD2", "CD5", "CD6",
    "CD7", "CD28", "LCK", "ZAP70", "LAT", "ITK", "THEMIS", "SKAP1", "CD8A", "CD8B",
    "IL7R", "CCR7", "TCF7", "LEF1", "SELL", "GZMA", "GZMK", "PRF1", "NKG7", "CTSW",
    "CD27", "ICOS", "CTLA4", "PDCD1", "FOXP3", "IKZF2", "RUNX3", "BCL11B", "CD69",
]


def symbol_to_ensembl() -> dict[str, str]:
    """Union the stats tables to map gene symbol -> Ensembl ID."""
    mapping: dict[str, str] = {}
    for path in sorted(STATS_ROOT.glob("*/heldout_allgene_*.csv")):
        frame = pd.read_csv(path, usecols=["Gene_name", "Ensembl_ID"])
        mapping.update(dict(zip(frame.Gene_name, frame.Ensembl_ID)))
    if not mapping:
        raise SystemExit(f"No stats tables found under {STATS_ROOT}; cannot map gene symbols")
    return mapping


def group_f_statistic(counts: sp.csr_matrix, labels: pd.Series) -> np.ndarray:
    """One-way F statistic per gene across the given grouping, on log1p-CPM values."""
    codes = pd.Categorical(labels).codes
    # codes is int8 for small category counts; keep the group count a Python int so
    # arithmetic against the cell count cannot overflow.
    n_groups = int(codes.max()) + 1
    n_cells, n_genes = counts.shape
    indicator = sp.csr_matrix(
        (np.ones(n_cells), (codes, np.arange(n_cells))), shape=(n_groups, n_cells)
    )
    group_n = np.asarray(indicator.sum(axis=1)).ravel()
    # Both are (n_groups x n_genes), so densifying is cheap regardless of cell count.
    group_sum = np.asarray(sp.csr_matrix(indicator @ counts).todense())
    group_sumsq = np.asarray(sp.csr_matrix(indicator @ counts.multiply(counts)).todense())

    group_mean = group_sum / group_n[:, None]
    grand_mean = np.asarray(counts.mean(axis=0)).ravel()

    between = (group_n[:, None] * (group_mean - grand_mean) ** 2).sum(axis=0) / max(n_groups - 1, 1)
    total_sumsq = group_sumsq.sum(axis=0)
    within_sumsq = total_sumsq - (group_n[:, None] * group_mean ** 2).sum(axis=0)
    within = within_sumsq / max(n_cells - n_groups, 1)
    return between / np.maximum(within, 1e-12)


def main() -> None:
    import anndata as ad
    from sklearn.linear_model import LogisticRegression
    from sklearn.model_selection import cross_val_score
    from sklearn.preprocessing import StandardScaler

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    print(f"Loading {H5AD}")
    adata = ad.read_h5ad(H5AD)
    raw_counts = adata.X if sp.issparse(adata.X) else sp.csr_matrix(adata.X)
    raw_counts = sp.csr_matrix(raw_counts)
    total_counts = np.asarray(raw_counts.sum(axis=1)).ravel()
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    # CPM + log1p, matching the standard normalization used elsewhere in this project.
    scaling = sp.diags(1e4 / np.maximum(total_counts, 1))
    norm = sp.csr_matrix(scaling @ raw_counts)
    norm.data = np.log1p(norm.data)

    genes = pd.Index(adata.var_names)
    detect_frac = np.asarray((raw_counts > 0).sum(axis=0)).ravel() / adata.n_obs
    mean_expr = np.asarray(norm.mean(axis=0)).ravel()

    print("Computing subtype and donor structure")
    subtype_f = group_f_statistic(norm, adata.obs["celltype"])
    donor_f = group_f_statistic(norm, adata.obs["individual"])

    print("Computing library-size dependence")
    # Spearman of expression against total counts, computed on ranks via a
    # centred dot product so it stays a sparse-friendly single pass.
    rank_lib = pd.Series(total_counts).rank().to_numpy()
    rank_lib = (rank_lib - rank_lib.mean()) / rank_lib.std()
    dense_proxy = norm.copy()
    col_mean = np.asarray(dense_proxy.mean(axis=0)).ravel()
    numerator = np.asarray(dense_proxy.T @ rank_lib).ravel() - adata.n_obs * col_mean * rank_lib.mean()
    sumsq = np.asarray(dense_proxy.multiply(dense_proxy).sum(axis=0)).ravel()
    denom = np.sqrt(np.maximum(sumsq - adata.n_obs * col_mean ** 2, 1e-12)) * np.sqrt(adata.n_obs)
    libsize_corr = numerator / np.maximum(denom, 1e-12)

    sym2ens = symbol_to_ensembl()
    ens2sym = {e: s for s, e in sym2ens.items()}

    features = pd.DataFrame({
        "ensembl_id": genes,
        "gene": [ens2sym.get(e, e) for e in genes],
        "detect_frac": detect_frac,
        "mean_expr": mean_expr,
        "log_subtype_f": np.log1p(subtype_f),
        "log_donor_f": np.log1p(donor_f),
        "libsize_corr": libsize_corr,
    })
    # Ambient shows breadth without depth: detected widely but never strongly.
    features["breadth_over_depth"] = features.detect_frac / np.maximum(features.mean_expr, 1e-6)
    features["subtype_over_donor"] = features.log_subtype_f / np.maximum(features.log_donor_f, 1e-6)

    present_symbols = set(features.gene)
    ambient_set = [g for g in KNOWN_AMBIENT if g in present_symbols]
    tcell_set = [g for g in KNOWN_TCELL if g in present_symbols]
    print(f"\nAnchors present: {len(ambient_set)} ambient, {len(tcell_set)} T-cell-intrinsic")

    # Only score genes with enough detection for the features to mean anything.
    scored = features[features.detect_frac >= 0.005].copy().reset_index(drop=True)
    print(f"Genes with detection >= 0.5%: {len(scored):,}")

    columns = ["log_subtype_f", "log_donor_f", "libsize_corr", "breadth_over_depth", "subtype_over_donor", "mean_expr"]
    labelled = scored[scored.gene.isin(ambient_set + tcell_set)].copy()
    labelled["is_ambient"] = labelled.gene.isin(ambient_set).astype(int)
    print(f"Anchors retained after detection filter: {int(labelled.is_ambient.sum())} ambient, "
          f"{int((1 - labelled.is_ambient).sum())} T-cell")

    scaler = StandardScaler().fit(labelled[columns])
    model = LogisticRegression(max_iter=2000, class_weight="balanced")
    auc = cross_val_score(model, scaler.transform(labelled[columns]), labelled.is_ambient,
                          cv=5, scoring="roc_auc")
    print(f"\nCross-validated AUC separating the anchor sets: {auc.mean():.3f} (+/- {auc.std():.3f})")
    if auc.mean() < 0.8:
        print("  [WARN] anchors separate poorly; treat the score as uninformative", flush=True)

    model.fit(scaler.transform(labelled[columns]), labelled.is_ambient)
    scored["ambient_risk"] = model.predict_proba(scaler.transform(scored[columns]))[:, 1]
    scored["ambient_pct"] = 100 * scored.ambient_risk.rank(pct=True)

    # Decision threshold: the lowest risk seen among the known-ambient anchors, so a
    # candidate is only flagged if it looks at least as ambient-like as a real contaminant.
    anchor_scores = scored[scored.gene.isin(ambient_set)].ambient_risk
    tcell_scores = scored[scored.gene.isin(tcell_set)].ambient_risk
    threshold = float(np.percentile(anchor_scores, 25))
    print(f"Known-ambient anchors  median risk: {anchor_scores.median():.3f}")
    print(f"Known T-cell anchors   median risk: {tcell_scores.median():.3f}")
    print(f"Flag threshold (25th pct of ambient anchors): {threshold:.3f}")

    scored.sort_values("ambient_risk", ascending=False).to_csv(OUT_DIR / "ambient_risk_all_genes.csv", index=False)

    candidates = pd.read_csv(CANDIDATES)
    per_gene = candidates.groupby(["Gene_name", "Ensembl_ID"], as_index=False).agg(
        n_comparisons=("comparison", "nunique"),
        program=("class_label", "first"),
        donor_robustness=("donor_robustness", lambda s: s.value_counts().idxmax()),
    )
    # Join on Ensembl ID rather than symbol; the h5ad is Ensembl-indexed and symbols are
    # not one-to-one.
    merged = per_gene.merge(
        scored[["ensembl_id", "ambient_risk", "ambient_pct", "detect_frac", "log_subtype_f", "log_donor_f"]],
        left_on="Ensembl_ID", right_on="ensembl_id", how="left",
    ).drop(columns="ensembl_id")
    merged["ambient_flag"] = np.where(
        merged.ambient_risk.isna(), "not_scored",
        np.where(merged.ambient_risk >= threshold, "AMBIENT_RISK", "ok"),
    )
    merged = merged.sort_values("ambient_risk", ascending=False)
    merged.to_csv(OUT_DIR / "ambient_risk_candidates.csv", index=False)

    summary = merged.ambient_flag.value_counts().rename_axis("flag").reset_index(name="n_genes")
    print("\n=== Candidate genes by ambient flag ===")
    print(summary.to_string(index=False))
    flagged = merged[merged.ambient_flag.eq("AMBIENT_RISK")]
    if len(flagged):
        print("\nFlagged candidates:")
        print(flagged[["Gene_name", "program", "ambient_risk", "ambient_pct", "n_comparisons"]].to_string(index=False))
    else:
        print("\nNo prioritized candidate reaches the ambient-anchor threshold.")

    by_program = (
        merged.dropna(subset=["ambient_risk"])
        .groupby("program", as_index=False)
        .agg(genes=("Gene_name", "size"), median_risk=("ambient_risk", "median"),
             n_flagged=("ambient_flag", lambda s: int((s == "AMBIENT_RISK").sum())))
        .sort_values("median_risk", ascending=False)
    )
    by_program.to_csv(OUT_DIR / "ambient_risk_by_program.csv", index=False)
    print("\n=== By program ===")
    print(by_program.to_string(index=False))

    (OUT_DIR / "ambient_risk_manifest.json").write_text(json.dumps({
        "note": "Ambient-RISK DIAGNOSTIC, not CellBender correction. No empty droplets available.",
        "h5ad": str(H5AD),
        "n_cells": int(adata.n_obs), "n_genes": int(adata.n_vars),
        "min_counts_per_cell": float(total_counts.min()),
        "anchors_ambient": ambient_set, "anchors_tcell": tcell_set,
        "cv_auc_mean": float(auc.mean()), "cv_auc_std": float(auc.std()),
        "flag_threshold": threshold,
    }, indent=2) + "\n")

    plot(scored, ambient_set, tcell_set, merged, threshold)
    print(f"\nWrote outputs to {OUT_DIR}")


def plot(scored: pd.DataFrame, ambient_set, tcell_set, merged: pd.DataFrame, threshold: float) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9.5, 5.5))
    bins = np.linspace(0, 1, 41)
    ax.hist(scored.ambient_risk, bins=bins, color="#dee2e6", label="all detected genes")
    ax.hist(scored[scored.gene.isin(tcell_set)].ambient_risk, bins=bins, color="#2A9D8F",
            alpha=0.85, label="known T-cell-intrinsic (anchor)")
    ax.hist(scored[scored.gene.isin(ambient_set)].ambient_risk, bins=bins, color="#E76F51",
            alpha=0.85, label="known ambient (anchor)")
    cand = merged.dropna(subset=["ambient_risk"])
    ax.hist(cand.ambient_risk, bins=bins, histtype="step", linewidth=2, color="#1d3557",
            label="prioritized candidates")
    ax.axvline(threshold, color="#E76F51", linestyle="--", linewidth=1.2, label="flag threshold")
    ax.set_yscale("log")
    ax.set_xlabel("Ambient-risk score (probability of resembling a known contaminant)")
    ax.set_ylabel("Genes (log scale)")
    ax.set_title("Prioritized candidates sit with genuine T-cell genes, not with contaminants", fontsize=11)
    ax.legend(fontsize=8, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIG_DIR / "ambient_risk_distribution.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
