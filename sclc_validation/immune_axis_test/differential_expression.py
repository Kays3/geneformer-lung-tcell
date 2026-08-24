#!/usr/bin/env python3
"""T5 -- genome-wide differential expression (PLAN.md section 6, T5a/T5b).

Two variants of the same script, selected by the POPULATION environment
variable, because they answer different questions and deliberately use
different cell populations (see PLAN.md's "Why two variants" for the full
rationale -- in short: this is a raw-counts test with no dependency on the
fine-tuned model, so it does not inherit the ISP pipeline's train/test
decision boundary, except for the one goal that specifically needs population
parity with the ISP hit lists):

  POPULATION=complete  (T5a, default) -- all 46,140 T cells (train+eval+test).
    Used for: where the 7 pre-registered exhaustion genes rank among all
    genome-wide DE genes; building a data-driven dysfunction score.
  POPULATION=test_only (T5b) -- the same 9,377 held-out test cells T2 used.
    Used only for: cross-referencing DE-significant genes against the
    whole-genome ISP hit lists, which were themselves computed from
    test-only cells.

Pairwise Wilcoxon DE (scanpy, the standard, well-tested implementation --
this is not hand-rolled) on log1p-CP10k for all three state pairs. Also
writes a per-donor pseudobulk mean-expression table (every gene, every donor)
so a downstream analysis can check whether an apparent DE signal is
donor-consistent or dominated by one donor, the same spirit as
`donor_consistency_allgene.py`'s check for ISP hits -- formal pseudobulk
significance testing is not attempted here since Normal has only 1 donor in
either population and SCLC/LUAD have too few for a parametric test to mean
much; the per-donor table is left for `analyze_differential_expression.py`
(and a human) to inspect directly.

Runs on the compute host against the prepared H5AD; no GPU, no perturbation
output needed. Given T2's measured runtime for a much smaller gene/cell
subset (PLAN.md section 9), expect seconds to low minutes here, not hours.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

H5AD = Path(os.environ.get(
    "HTAN_H5AD",
    Path.home() / "workspace/KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad",
))
POPULATION = os.environ.get("POPULATION", "complete")
OUT_DIR = Path(__file__).resolve().parent / "results"

DISEASE_LABEL = {"normal": "normal", "small cell lung carcinoma": "sclc", "lung adenocarcinoma": "luad"}
CD4CD8_GROUP = {
    "CD4-positive helper T cell": "CD4",
    "CD8-positive, alpha-beta memory T cell": "CD8",
    "effector CD8-positive, alpha-beta T cell": "CD8",
    "regulatory T cell": "CD4 (Treg)",
    "exhausted T cell": "other",
    "gamma-delta T cell": "other",
}
PAIRS = [("sclc", "luad"), ("sclc", "normal"), ("luad", "normal")]

if POPULATION not in ("complete", "test_only"):
    raise SystemExit(f"POPULATION must be 'complete' or 'test_only', got {POPULATION!r}")


def load_population() -> ad.AnnData:
    adata = ad.read_h5ad(H5AD, backed="r")
    mask = (adata.obs["split"] == "test") if POPULATION == "test_only" else pd.Series(True, index=adata.obs.index)
    sub = adata[mask].to_memory()
    sub.obs["state"] = sub.obs["disease"].map(DISEASE_LABEL).astype("category")
    sub.obs["cd4cd8"] = sub.obs["celltype"].map(CD4CD8_GROUP).fillna("other").astype("category")
    return sub


def run_de(sub: ad.AnnData) -> pd.DataFrame:
    sc.pp.normalize_total(sub, target_sum=1e4)
    sc.pp.log1p(sub)

    rows = []
    for group, reference in PAIRS:
        sc.tl.rank_genes_groups(sub, groupby="state", groups=[group], reference=reference, method="wilcoxon")
        res = sc.get.rank_genes_groups_df(sub, group=group)
        res["comparison"] = f"{group}_vs_{reference}"
        rows.append(res)
    de = pd.concat(rows, ignore_index=True)
    return de.rename(columns={
        "names": "gene", "logfoldchanges": "log2fc", "pvals": "pval",
        "pvals_adj": "padj", "scores": "wilcoxon_score",
    })[["comparison", "gene", "log2fc", "pval", "padj", "wilcoxon_score"]]


def per_donor_pseudobulk(sub: ad.AnnData) -> pd.DataFrame:
    """Mean log1p-CP10k expression per gene, per donor -- for inspecting
    whether a DE hit is donor-consistent, not a formal significance test."""
    rows = []
    for (state, donor), idx in sub.obs.groupby(["state", "individual"], observed=True).indices.items():
        idx = np.asarray(idx)
        block = sub.X[idx, :]
        n_cells = len(idx)
        mean_expr = np.asarray(block.mean(axis=0)).ravel()
        rows.append(pd.DataFrame({
            "state": state, "donor": donor, "n_cells": n_cells,
            "gene": sub.var_names, "mean_log1p_cp10k": mean_expr,
        }))
    return pd.concat(rows, ignore_index=True)


def state_pooled_expression(sub: ad.AnnData) -> pd.DataFrame:
    """Genome-wide state-level mean expression and detection rate, pooling all
    cells of that state directly (not an average of donor means) -- same
    convention as measure_baseline_expression.py's pooled table, extended here
    to every gene instead of the 176-gene program/HK subset."""
    rows = []
    for state, idx in sub.obs.groupby("state", observed=True).indices.items():
        idx = np.asarray(idx)
        block = sub.X[idx, :]
        n_cells = len(idx)
        n_donors = sub.obs.loc[sub.obs["state"] == state, "individual"].nunique()
        detect_rate = np.asarray((block > 0).sum(axis=0)).ravel() / n_cells
        mean_expr = np.asarray(block.mean(axis=0)).ravel()
        rows.append(pd.DataFrame({
            "state": state, "n_cells": n_cells, "n_donors": n_donors,
            "gene": sub.var_names, "detect_rate": detect_rate,
            "mean_log1p_cp10k": mean_expr,
        }))
    return pd.concat(rows, ignore_index=True)


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    sub = load_population()
    print(f"POPULATION={POPULATION}: {sub.n_obs} cells x {sub.n_vars} genes")
    print(sub.obs.groupby("state", observed=True).agg(n_cells=("individual", "size"), n_donors=("individual", "nunique")))

    de = run_de(sub)
    de_path = OUT_DIR / f"differential_expression_{POPULATION}.csv"
    de.to_csv(de_path, index=False)
    print(f"\nWrote {de_path} ({len(de)} rows)")

    pseudobulk = per_donor_pseudobulk(sub)
    pb_path = OUT_DIR / f"pseudobulk_per_donor_{POPULATION}.csv"
    pseudobulk.to_csv(pb_path, index=False)
    print(f"Wrote {pb_path} ({len(pseudobulk)} rows)")

    pooled = state_pooled_expression(sub)
    pooled_path = OUT_DIR / f"state_pooled_expression_{POPULATION}.csv"
    pooled.to_csv(pooled_path, index=False)
    print(f"Wrote {pooled_path} ({len(pooled)} rows)")

    sig = de[de["padj"] < 0.05]
    print(f"\nSignificant (padj<0.05) genes per comparison:")
    print(sig.groupby("comparison").size())


if __name__ == "__main__":
    main()
