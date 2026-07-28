#!/usr/bin/env python3
"""GSE263196 spatial validation: test whether the pre-registered T-cell
dysfunction signature is enriched in T-cell-rich regions of SCLC Visium
spots.

No cell-type labels are bundled with this dataset, so T-cell abundance is
quantified with a marker-gene score per spot (not formal deconvolution --
an explicit scope decision, since the audit's design allows either
"deconvolve or annotate"). Effect sizes (Spearman rho with 95% CI) are
reported per sample and combined across the 5 patient samples via a
Fisher z-weighted meta-analysis, per the audit's "patient-level estimates,
not spot-level P values alone" guidance.
"""
from __future__ import annotations

import gzip
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
import scipy.stats as st
import anndata as ad

HERE = Path(__file__).resolve().parent
RAW_DIR = HERE.parent / "audit" / "source_metadata" / "GSE263196_RAW"
RESULTS_DIR = HERE / "results"
FIGURE_DIR = HERE / "figures"

SAMPLES = {
    "GSM8187469": "SCLC3",
    "GSM8187470": "SCLC4",
    "GSM8187471": "SCLC8",
    "GSM8187472": "SCLC9",
    "GSM8187473": "SCLC12",
}

# Pan-T-cell markers (abundance / identity), deliberately distinct from the
# dysfunction panel below -- this is "is a T cell here", not "how exhausted".
TCELL_MARKERS = ["CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD28", "TRBC1", "TRBC2", "IL7R", "CD8A", "CD8B", "CD4"]

# Exhaustion/dysfunction subset of the audit's pre-registered 21-gene panel.
DYSFUNCTION_MARKERS = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]

MIN_COUNTS_PER_SPOT = 200


def load_sample(gsm: str, label: str) -> ad.AnnData:
    matrix = scipy.io.mmread(RAW_DIR / f"{gsm}_{label}_matrix.mtx.gz").tocsr().T.tocsr()
    features = pd.read_csv(
        RAW_DIR / f"{gsm}_{label}_features.tsv.gz", sep="\t", header=None,
        names=["ensembl_id", "gene_symbol", "feature_type"],
    )
    barcodes = pd.read_csv(RAW_DIR / f"{gsm}_{label}_barcodes.tsv.gz", header=None, names=["barcode"])
    positions = pd.read_csv(
        RAW_DIR / f"{gsm}_{label}_tissue_positions_list.csv.gz", header=None,
        names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"],
    )

    obs = barcodes.merge(positions, on="barcode", how="left").set_index("barcode")
    var = features.set_index("gene_symbol")
    var["ensembl_id"] = var["ensembl_id"].astype(str)
    var_names_dedup = pd.Index(var.index).where(~pd.Index(var.index).duplicated(), var.index + "_dup")
    var.index = var_names_dedup

    a = ad.AnnData(X=matrix, obs=obs, var=var)
    a.obs["sample_gsm"] = gsm
    a.obs["sample_label"] = label
    a.obs_names_make_unique()

    a = a[a.obs["in_tissue"] == 1].copy()
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None)
    a = a[a.obs["total_counts"] >= MIN_COUNTS_PER_SPOT].copy()
    return a


def score_sample(a: ad.AnnData) -> ad.AnnData:
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)

    tcell_present = [g for g in TCELL_MARKERS if g in a.var_names]
    dysfunction_present = [g for g in DYSFUNCTION_MARKERS if g in a.var_names]
    missing_tcell = sorted(set(TCELL_MARKERS) - set(tcell_present))
    missing_dys = sorted(set(DYSFUNCTION_MARKERS) - set(dysfunction_present))
    if missing_tcell:
        print(f"  [warn] missing T-cell markers: {missing_tcell}")
    if missing_dys:
        print(f"  [warn] missing dysfunction markers: {missing_dys}")

    sc.tl.score_genes(a, tcell_present, score_name="tcell_score")
    sc.tl.score_genes(a, dysfunction_present, score_name="dysfunction_score")
    return a


def spearman_ci(x: np.ndarray, y: np.ndarray, alpha: float = 0.05) -> dict:
    rho, p = st.spearmanr(x, y)
    n = len(x)
    z = np.arctanh(rho)
    se = 1 / np.sqrt(n - 3)
    z_crit = st.norm.ppf(1 - alpha / 2)
    lo, hi = np.tanh(z - z_crit * se), np.tanh(z + z_crit * se)
    return {"rho": rho, "p_value": p, "ci_low": lo, "ci_high": hi, "n_spots": n, "fisher_z": z, "se": se}


def meta_analyze(per_sample: list[dict]) -> dict:
    weights = np.array([1 / r["se"] ** 2 for r in per_sample])
    z_values = np.array([r["fisher_z"] for r in per_sample])
    z_pooled = np.sum(weights * z_values) / np.sum(weights)
    se_pooled = 1 / np.sqrt(np.sum(weights))
    z_crit = st.norm.ppf(0.975)
    rho_pooled = np.tanh(z_pooled)
    lo, hi = np.tanh(z_pooled - z_crit * se_pooled), np.tanh(z_pooled + z_crit * se_pooled)
    p_pooled = 2 * (1 - st.norm.cdf(abs(z_pooled / se_pooled)))
    return {"rho_pooled": rho_pooled, "ci_low": lo, "ci_high": hi, "p_value": p_pooled, "n_samples": len(per_sample)}


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)

    per_sample_rows = []
    all_spot_rows = []
    for gsm, label in SAMPLES.items():
        print(f"=== {gsm} ({label}) ===")
        a = load_sample(gsm, label)
        print(f"  {a.n_obs} in-tissue spots (>= {MIN_COUNTS_PER_SPOT} counts), {a.n_vars} genes")
        a = score_sample(a)

        stats = spearman_ci(a.obs["tcell_score"].values, a.obs["dysfunction_score"].values)
        stats.update({"sample_gsm": gsm, "sample_label": label})
        per_sample_rows.append(stats)
        print(f"  Spearman rho={stats['rho']:.3f} [{stats['ci_low']:.3f}, {stats['ci_high']:.3f}], p={stats['p_value']:.3e}, n={stats['n_spots']}")

        spot_df = a.obs[["sample_gsm", "sample_label", "total_counts", "tcell_score", "dysfunction_score"]].copy()
        spot_df["barcode"] = a.obs_names
        all_spot_rows.append(spot_df)

    per_sample_df = pd.DataFrame(per_sample_rows)[
        ["sample_gsm", "sample_label", "n_spots", "rho", "ci_low", "ci_high", "p_value"]
    ]
    per_sample_df.to_csv(RESULTS_DIR / "spatial_tcell_dysfunction_correlation_by_sample.csv", index=False)

    pooled = meta_analyze(per_sample_rows)
    pooled_df = pd.DataFrame([pooled])
    pooled_df.to_csv(RESULTS_DIR / "spatial_tcell_dysfunction_correlation_pooled.csv", index=False)

    spot_level = pd.concat(all_spot_rows, ignore_index=True)
    spot_level.to_csv(RESULTS_DIR / "spatial_spot_level_scores.csv.gz", index=False, compression="gzip")

    print("\n=== Per-sample results ===")
    print(per_sample_df.to_string(index=False))
    print("\n=== Pooled (5-sample meta-analysis) ===")
    print(pooled_df.to_string(index=False))

    plot_forest(per_sample_df, pooled)


def plot_forest(per_sample_df: pd.DataFrame, pooled: dict) -> None:
    rows = per_sample_df.sort_values("rho").to_dict(orient="records")
    labels = [f"{r['sample_label']} (n={r['n_spots']:,})" for r in rows] + ["Pooled (5 samples)"]
    rhos = [r["rho"] for r in rows] + [pooled["rho_pooled"]]
    los = [r["ci_low"] for r in rows] + [pooled["ci_low"]]
    his = [r["ci_high"] for r in rows] + [pooled["ci_high"]]
    y = np.arange(len(labels))

    fig, ax = plt.subplots(figsize=(8.5, 3.8))
    ax.axvline(0, color="#9aa5a3", linewidth=1, linestyle="--", zorder=1)
    colors = ["#457B9D"] * (len(labels) - 1) + ["#E76F51"]
    for i, (r, lo, hi, c) in enumerate(zip(rhos, los, his, colors)):
        ax.plot([lo, hi], [i, i], color=c, linewidth=2, zorder=2)
        ax.scatter([r], [i], color=c, s=60, zorder=3, edgecolor="white", linewidth=0.8)
    ax.set_yticks(y, labels)
    ax.set_xlabel("Spearman ρ (T-cell abundance score vs. dysfunction score)")
    ax.set_title("Dysfunction-signature enrichment in T-cell-rich SCLC Visium spots", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    plt.tight_layout()
    fig.savefig(FIGURE_DIR / "tcell_dysfunction_correlation_forest.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
