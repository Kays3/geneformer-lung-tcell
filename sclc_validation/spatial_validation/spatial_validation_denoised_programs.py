#!/usr/bin/env python3
"""GSE263196 spatial validation of the denoised immune/cancer perturbation programs.

This reuses the pipeline in `spatial_validation.py` -- same spot QC, same
`scanpy.tl.score_genes` scoring, same Spearman + Fisher-z meta-analysis across the
5 patient samples -- but tests the gene programs prioritized by the all-gene
perturbation screen after technical-noise removal, rather than the single
pre-registered 7-gene dysfunction panel.

Three guards make this a real test rather than a guaranteed positive:

1. Circularity guard. Genes shared with the T-cell abundance marker set are removed
   from every tested program. Without this, the "T-cell identity / TCR" program would
   trivially correlate with T-cell abundance because it is largely the same genes.

2. Library-size control. MHC-I genes are expressed by nearly every cell, so a raw
   correlation can reflect spot cellularity rather than immune co-localization. A
   partial Spearman controlling for log10 total counts is reported alongside the raw one.

3. Expression-matched random null. For each program, random gene sets matched on mean
   expression decile are scored the same way. With ~15,000 spots almost any gene set
   reaches significance, so the informative question is whether a program is more
   correlated with T-cell abundance than equally-expressed random genes.

The original 7-gene dysfunction panel is re-run as a benchmark; reproducing its
published pooled rho of 0.161 is the pipeline-equivalence check.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.stats as st
import anndata as ad

HERE = Path(__file__).resolve().parent
RAW_DIR = Path(os.environ.get("GSE263196_RAW_DIR", HERE.parent / "audit" / "source_metadata" / "GSE263196_RAW"))
CANDIDATES = Path(
    os.environ.get(
        "DENOISED_CANDIDATES",
        HERE.parent / "primary_test_perturbation" / "tables" / "immune_cancer_candidates.csv",
    )
)
RESULTS_DIR = HERE / "results_denoised_programs"
FIGURE_DIR = HERE / "figures"

SAMPLES = {
    "GSM8187469": "SCLC3",
    "GSM8187470": "SCLC4",
    "GSM8187471": "SCLC8",
    "GSM8187472": "SCLC9",
    "GSM8187473": "SCLC12",
}

# Unchanged from spatial_validation.py -- this defines the x-axis (is a T cell here).
TCELL_MARKERS = ["CD3D", "CD3E", "CD3G", "CD2", "CD5", "CD28", "TRBC1", "TRBC2", "IL7R", "CD8A", "CD8B", "CD4"]

# Benchmark: the pre-registered panel used in the published spatial result.
DYSFUNCTION_MARKERS = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]

MIN_COUNTS_PER_SPOT = 200
MIN_GENES_PER_PROGRAM = 3
MIN_COMPARISONS = 2
N_NULL = int(os.environ.get("N_NULL", 100))
RANDOM_SEED = 43


def load_programs() -> tuple[dict[str, list[str]], pd.DataFrame]:
    """Build gene programs from the denoised perturbation candidates."""
    candidates = pd.read_csv(CANDIDATES)
    recurrent = (
        candidates.groupby(["Gene_name", "class_label"], as_index=False)
        .agg(n_comparisons=("comparison", "nunique"))
        .query("n_comparisons >= @MIN_COMPARISONS")
    )
    programs: dict[str, list[str]] = {}
    audit_rows = []
    for label, group in recurrent.groupby("class_label"):
        genes = sorted(group.Gene_name)
        dropped = sorted(set(genes) & set(TCELL_MARKERS))
        kept = [g for g in genes if g not in TCELL_MARKERS]
        audit_rows.append(
            {
                "program": label,
                "n_recurrent_genes": len(genes),
                "n_dropped_circular": len(dropped),
                "dropped_genes": "; ".join(dropped),
                "n_tested": len(kept),
                "tested_genes": "; ".join(kept),
                "status": "tested" if len(kept) >= MIN_GENES_PER_PROGRAM else "skipped_too_few_genes",
            }
        )
        if len(kept) >= MIN_GENES_PER_PROGRAM:
            programs[label] = kept
    programs["Dysfunction panel (benchmark)"] = DYSFUNCTION_MARKERS
    audit_rows.append(
        {
            "program": "Dysfunction panel (benchmark)",
            "n_recurrent_genes": len(DYSFUNCTION_MARKERS),
            "n_dropped_circular": 0,
            "dropped_genes": "",
            "n_tested": len(DYSFUNCTION_MARKERS),
            "tested_genes": "; ".join(DYSFUNCTION_MARKERS),
            "status": "benchmark",
        }
    )
    return programs, pd.DataFrame(audit_rows)


def load_sample(gsm: str, label: str) -> ad.AnnData:
    """Identical spot loading and QC to spatial_validation.py."""
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
    var.index = pd.Index(var.index).where(~pd.Index(var.index).duplicated(), var.index + "_dup")

    a = ad.AnnData(X=matrix, obs=obs, var=var)
    a.obs["sample_gsm"] = gsm
    a.obs["sample_label"] = label
    a.obs_names_make_unique()
    a = a[a.obs["in_tissue"] == 1].copy()
    sc.pp.calculate_qc_metrics(a, inplace=True, percent_top=None)
    a = a[a.obs["total_counts"] >= MIN_COUNTS_PER_SPOT].copy()

    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    return a


def spearman_ci(x: np.ndarray, y: np.ndarray, n_controls: int = 0, alpha: float = 0.05) -> dict:
    rho, p = st.spearmanr(x, y)
    n = len(x)
    z = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1 / np.sqrt(n - 3 - n_controls)
    z_crit = st.norm.ppf(1 - alpha / 2)
    return {
        "rho": rho, "p_value": p, "ci_low": np.tanh(z - z_crit * se), "ci_high": np.tanh(z + z_crit * se),
        "n_spots": n, "fisher_z": z, "se": se,
    }


def partial_spearman(x: np.ndarray, y: np.ndarray, z: np.ndarray) -> dict:
    """Spearman of x and y with z partialled out, computed on ranks."""
    rx, ry, rz = st.rankdata(x), st.rankdata(y), st.rankdata(z)
    rz_design = np.column_stack([np.ones_like(rz), rz])
    res_x = rx - rz_design @ np.linalg.lstsq(rz_design, rx, rcond=None)[0]
    res_y = ry - rz_design @ np.linalg.lstsq(rz_design, ry, rcond=None)[0]
    rho = float(np.corrcoef(res_x, res_y)[0, 1])
    n = len(x)
    z_t = np.arctanh(np.clip(rho, -0.999999, 0.999999))
    se = 1 / np.sqrt(n - 4)
    z_crit = st.norm.ppf(0.975)
    p = 2 * (1 - st.norm.cdf(abs(z_t / se)))
    return {
        "rho": rho, "p_value": p, "ci_low": np.tanh(z_t - z_crit * se), "ci_high": np.tanh(z_t + z_crit * se),
        "n_spots": n, "fisher_z": z_t, "se": se,
    }


def meta_analyze(per_sample: list[dict]) -> dict:
    weights = np.array([1 / r["se"] ** 2 for r in per_sample])
    z_values = np.array([r["fisher_z"] for r in per_sample])
    z_pooled = float(np.sum(weights * z_values) / np.sum(weights))
    se_pooled = float(1 / np.sqrt(np.sum(weights)))
    z_crit = st.norm.ppf(0.975)
    # Cochran's Q quantifies between-patient heterogeneity, which the earlier
    # single-panel result flagged qualitatively but never tested.
    q = float(np.sum(weights * (z_values - z_pooled) ** 2))
    df = len(per_sample) - 1
    return {
        "rho_pooled": float(np.tanh(z_pooled)),
        "ci_low": float(np.tanh(z_pooled - z_crit * se_pooled)),
        "ci_high": float(np.tanh(z_pooled + z_crit * se_pooled)),
        "p_value": float(2 * (1 - st.norm.cdf(abs(z_pooled / se_pooled)))),
        "n_samples": len(per_sample),
        "q_statistic": q,
        "q_p_value": float(1 - st.chi2.cdf(q, df)) if df > 0 else np.nan,
        "i_squared": float(max(0.0, (q - df) / q) * 100) if q > 0 else 0.0,
    }


def expression_bins(adatas: dict[str, ad.AnnData], n_bins: int = 10) -> dict[str, pd.Series]:
    """Assign each gene to a mean-expression bin, per sample, for matched null sampling."""
    bins = {}
    for label, a in adatas.items():
        mean_expression = np.asarray(a.X.mean(axis=0)).ravel()
        series = pd.Series(mean_expression, index=a.var_names)
        ranked = series.rank(method="first")
        bins[label] = pd.cut(ranked, bins=n_bins, labels=False).astype(int)
    return bins


def sample_matched_genes(program: list[str], gene_bins: pd.Series, rng: np.random.Generator) -> list[str]:
    """Draw a random gene set matched to the program's expression-bin composition."""
    present = [g for g in program if g in gene_bins.index]
    picks = []
    for gene in present:
        pool = gene_bins.index[gene_bins.values == gene_bins[gene]]
        picks.append(str(rng.choice(pool)))
    return picks


def score(a: ad.AnnData, genes: list[str], name: str) -> np.ndarray | None:
    present = [g for g in genes if g in a.var_names]
    if len(present) < MIN_GENES_PER_PROGRAM:
        return None
    sc.tl.score_genes(a, present, score_name=name)
    return a.obs[name].to_numpy()


def main() -> None:
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(RANDOM_SEED)

    programs, program_audit = load_programs()
    program_audit.to_csv(RESULTS_DIR / "program_definitions.csv", index=False)
    print(f"Testing {len(programs)} programs (>= {MIN_GENES_PER_PROGRAM} genes after circularity guard)")
    for name, genes in programs.items():
        print(f"  {name}: {len(genes)} genes")

    print("\nLoading samples...")
    adatas: dict[str, ad.AnnData] = {}
    for gsm, label in SAMPLES.items():
        a = load_sample(gsm, label)
        tcell_present = [g for g in TCELL_MARKERS if g in a.var_names]
        sc.tl.score_genes(a, tcell_present, score_name="tcell_score")
        adatas[label] = a
        print(f"  {label}: {a.n_obs} spots, {a.n_vars} genes, {len(tcell_present)}/{len(TCELL_MARKERS)} T-cell markers")

    gene_bins = expression_bins(adatas)

    per_sample_rows, pooled_rows, coverage_rows = [], [], []
    for name, genes in programs.items():
        print(f"\n=== {name} ===")
        raw_stats, partial_stats = [], []
        for label, a in adatas.items():
            values = score(a, genes, "program_score")
            present = [g for g in genes if g in a.var_names]
            coverage_rows.append(
                {"program": name, "sample_label": label, "n_requested": len(genes), "n_present": len(present),
                 "missing_genes": "; ".join(sorted(set(genes) - set(present)))}
            )
            if values is None:
                print(f"  [warn] {label}: too few genes present, skipped")
                continue
            tcell = a.obs["tcell_score"].to_numpy()
            log_counts = np.log10(a.obs["total_counts"].to_numpy())

            raw = spearman_ci(tcell, values)
            partial = partial_spearman(tcell, values, log_counts)
            raw.update({"program": name, "sample_label": label, "model": "raw"})
            partial.update({"program": name, "sample_label": label, "model": "partial_logcounts"})
            raw_stats.append(raw)
            partial_stats.append(partial)
            per_sample_rows.extend([raw, partial])
            print(f"  {label}: raw rho={raw['rho']:+.3f}  partial rho={partial['rho']:+.3f}  n={raw['n_spots']}")

        if not raw_stats:
            continue

        # Expression-matched null on the raw statistic.
        null_rhos = []
        for _ in range(N_NULL):
            null_sample_stats = []
            for label, a in adatas.items():
                null_genes = sample_matched_genes(genes, gene_bins[label], rng)
                values = score(a, null_genes, "null_score")
                if values is None:
                    continue
                null_sample_stats.append(spearman_ci(a.obs["tcell_score"].to_numpy(), values))
            if null_sample_stats:
                null_rhos.append(meta_analyze(null_sample_stats)["rho_pooled"])
        null_rhos = np.array(null_rhos)

        pooled_raw = meta_analyze(raw_stats)
        pooled_partial = meta_analyze(partial_stats)
        observed = pooled_raw["rho_pooled"]
        # One-sided empirical p: how often does a matched random set reach this rho?
        empirical_p = float((np.sum(null_rhos >= observed) + 1) / (len(null_rhos) + 1))
        null_z = float((observed - null_rhos.mean()) / null_rhos.std()) if null_rhos.std() > 0 else np.nan

        row = {
            "program": name,
            "n_genes_tested": len(genes),
            "rho_pooled_raw": observed,
            "ci_low_raw": pooled_raw["ci_low"],
            "ci_high_raw": pooled_raw["ci_high"],
            "p_value_raw": pooled_raw["p_value"],
            "rho_pooled_partial": pooled_partial["rho_pooled"],
            "ci_low_partial": pooled_partial["ci_low"],
            "ci_high_partial": pooled_partial["ci_high"],
            "i_squared": pooled_raw["i_squared"],
            "q_p_value": pooled_raw["q_p_value"],
            "null_mean_rho": float(null_rhos.mean()),
            "null_sd_rho": float(null_rhos.std()),
            "null_z": null_z,
            "empirical_p_vs_null": empirical_p,
            "n_null_iterations": len(null_rhos),
            "exceeds_null": bool(empirical_p < 0.05),
        }
        pooled_rows.append(row)
        print(
            f"  POOLED raw={observed:+.3f} [{pooled_raw['ci_low']:+.3f}, {pooled_raw['ci_high']:+.3f}]  "
            f"partial={pooled_partial['rho_pooled']:+.3f}  "
            f"null={null_rhos.mean():+.3f}±{null_rhos.std():.3f}  z={null_z:+.2f}  emp_p={empirical_p:.3f}  "
            f"I2={pooled_raw['i_squared']:.0f}%"
        )

    per_sample_df = pd.DataFrame(per_sample_rows)[
        ["program", "sample_label", "model", "n_spots", "rho", "ci_low", "ci_high", "p_value"]
    ]
    per_sample_df.to_csv(RESULTS_DIR / "denoised_programs_by_sample.csv", index=False)
    pooled_df = pd.DataFrame(pooled_rows).sort_values("null_z", ascending=False)
    pooled_df.to_csv(RESULTS_DIR / "denoised_programs_pooled.csv", index=False)
    pd.DataFrame(coverage_rows).to_csv(RESULTS_DIR / "denoised_programs_gene_coverage.csv", index=False)

    manifest = {
        "source_cohort": "GSE263196 (5 fresh-frozen SCLC Visium samples)",
        "candidates_table": str(CANDIDATES),
        "min_counts_per_spot": MIN_COUNTS_PER_SPOT,
        "min_genes_per_program": MIN_GENES_PER_PROGRAM,
        "min_comparisons_for_recurrence": MIN_COMPARISONS,
        "n_null_iterations": N_NULL,
        "random_seed": RANDOM_SEED,
        "circularity_guard": "Genes overlapping the T-cell abundance marker set are removed from every tested program.",
        "programs_tested": {k: v for k, v in programs.items()},
    }
    (RESULTS_DIR / "denoised_programs_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")

    print("\n=== Pooled summary ===")
    print(
        pooled_df[
            ["program", "n_genes_tested", "rho_pooled_raw", "rho_pooled_partial", "null_mean_rho", "null_z",
             "empirical_p_vs_null", "exceeds_null", "i_squared"]
        ].to_string(index=False)
    )
    plot_program_forest(pooled_df)
    plot_null_comparison(pooled_df)


def plot_program_forest(pooled_df: pd.DataFrame) -> None:
    rows = pooled_df.sort_values("rho_pooled_raw").to_dict(orient="records")
    labels = [f"{r['program']} ({r['n_genes_tested']}g)" for r in rows]
    y = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(11.5, 0.55 * len(labels) + 2.4))
    ax.axvline(0, color="#9aa5a3", linewidth=1, linestyle="--", zorder=1)
    for i, r in enumerate(rows):
        colour = "#E76F51" if r["program"].startswith("Dysfunction") else ("#2A9D8F" if r["exceeds_null"] else "#8D99AE")
        ax.plot([r["ci_low_raw"], r["ci_high_raw"]], [i, i], color=colour, linewidth=2.4, zorder=3)
        ax.scatter([r["rho_pooled_raw"]], [i], color=colour, s=70, zorder=4, edgecolor="white", linewidth=0.8)
        ax.scatter([r["rho_pooled_partial"]], [i], color=colour, s=46, marker="D", zorder=4, alpha=0.65)
        ax.scatter([r["null_mean_rho"]], [i], color="#4a4e69", s=34, marker="|", zorder=5)
    ax.set_yticks(y, labels, fontsize=9)
    ax.set_xlabel("Pooled Spearman ρ vs. T-cell abundance score (5-sample Fisher-z meta-analysis)")
    ax.set_title(
        "Denoised perturbation programs in SCLC Visium tissue\n"
        "circle = raw ρ · diamond = partial ρ (counts controlled) · tick = random-null mean",
        fontsize=10.5,
    )
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.2)
    plt.tight_layout()
    fig.savefig(FIGURE_DIR / "denoised_programs_forest.png", dpi=200)
    plt.close(fig)


def plot_null_comparison(pooled_df: pd.DataFrame) -> None:
    df = pooled_df.sort_values("null_z")
    fig, ax = plt.subplots(figsize=(9.5, 0.5 * len(df) + 2.2))
    colours = ["#2A9D8F" if e else "#8D99AE" for e in df.exceeds_null]
    ax.barh(df.program, df.null_z, color=colours)
    ax.axvline(0, color="#4a4e69", linewidth=1)
    ax.axvline(1.96, color="#E76F51", linewidth=1, linestyle="--", label="z = 1.96")
    ax.set_xlabel("Standard deviations above the expression-matched random-gene-set null")
    ax.set_title("Is each program more T-cell-colocalized than equally-expressed random genes?", fontsize=11)
    ax.legend(frameon=False, fontsize=8)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", alpha=0.2)
    plt.tight_layout()
    fig.savefig(FIGURE_DIR / "denoised_programs_null_z.png", dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    main()
