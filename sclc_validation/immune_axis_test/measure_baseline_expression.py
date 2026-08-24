#!/usr/bin/env python3
"""T2 -- measured baseline expression (PLAN.md section 6).

The control the axis claim depends on and which had never been run: what is
the actual per-cell checkpoint expression in SCLC vs LUAD vs Normal T cells in
this atlas? Per-donor mean expression and detection rate for each of the four
pre-registered immune-axis programs, on the held-out test cells, stratified by
donor and by CD4/CD8.

Extended here (same atlas load, marginal cost) to the whole-genome
housekeeping-gene diagnostic's HK-flagged concordant hits
(../primary_test_perturbation/tables/hk_gene_diagnostic/allgene_hk_concordant_detection_gap.csv),
answering that report's own "run T2" follow-up in the same pass rather than a
second one.

This decides whether the model is describing the data or distorting it: if a
gene's detection/expression genuinely differs between states, an ISP shift for
it has support from the input; if it does not, the shift has none.

Runs on the compute host against the prepared H5AD (raw counts, `split`,
`individual`, `celltype`, `disease` in obs) -- no perturbation output needed,
no GPU.
"""
from __future__ import annotations

import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = Path(os.environ.get(
    "HTAN_H5AD",
    Path.home() / "workspace/KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad",
))
STATS_ROOT = Path(os.environ.get(
    "SCLC_PERTURBATION_ROOT",
    Path.home() / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation",
)) / "stats"
HK_GENE_LIST = Path(os.environ.get(
    "HK_GENE_LIST",
    Path(__file__).resolve().parent / "hk_concordant_gene_list.csv",
))
OUT_DIR = Path(__file__).resolve().parent / "results"

DISEASE_LABEL = {"normal": "normal", "small cell lung carcinoma": "sclc", "lung adenocarcinoma": "luad"}

# Same four programs and gene sets as axis_consistency.py, so T1 and T2 talk
# about the same genes.
PROGRAMS = {
    "exhaustion": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"],
    "cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG"],
    "progenitor": ["TCF7", "SLAMF6", "IL7R", "CCR7"],
    "sclc_subtype_tf": ["ASCL1", "NEUROD1", "POU2F3", "YAP1"],
}

CD4CD8_GROUP = {
    "CD4-positive helper T cell": "CD4",
    "CD8-positive, alpha-beta memory T cell": "CD8",
    "effector CD8-positive, alpha-beta T cell": "CD8",
    "regulatory T cell": "CD4 (Treg)",
    "exhausted T cell": "other",
    "gamma-delta T cell": "other",
}


def symbol_to_ensembl() -> dict[str, str]:
    """Same mapping ambient_risk_diagnostic.py uses: union the stats tables."""
    mapping: dict[str, str] = {}
    for path in sorted(STATS_ROOT.glob("*/heldout_allgene_*.csv")):
        frame = pd.read_csv(path, usecols=["Gene_name", "Ensembl_ID"])
        mapping.update(dict(zip(frame.Gene_name, frame.Ensembl_ID)))
    if not mapping:
        raise SystemExit(f"No stats tables found under {STATS_ROOT}; cannot map gene symbols")
    return mapping


def gene_list() -> tuple[list[str], dict[str, str]]:
    program_genes = sorted({g for genes in PROGRAMS.values() for g in genes})
    hk_genes: list[str] = []
    if HK_GENE_LIST.exists():
        hk_genes = sorted(pd.read_csv(HK_GENE_LIST)["Gene_name"].unique().tolist())
    else:
        print(f"[warn] {HK_GENE_LIST} not found; running program genes only")
    all_genes = sorted(set(program_genes) | set(hk_genes))
    gene_to_program = {g: name for name, genes in PROGRAMS.items() for g in genes}
    for g in hk_genes:
        gene_to_program.setdefault(g, "hk_concordant_whole_genome")
    return all_genes, gene_to_program


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    genes, gene_to_program = gene_list()
    sym2ens = symbol_to_ensembl()

    missing_mapping = [g for g in genes if g not in sym2ens]
    if missing_mapping:
        print(f"[warn] {len(missing_mapping)} genes have no Ensembl mapping, skipped: {missing_mapping}")
    genes = [g for g in genes if g in sym2ens]
    ensembl_ids = [sym2ens[g] for g in genes]

    adata = ad.read_h5ad(H5AD, backed="r")
    missing_var = [e for e in ensembl_ids if e not in adata.var_names]
    if missing_var:
        print(f"[warn] {len(missing_var)} Ensembl IDs not in H5AD var, skipped")
    keep = [(g, e) for g, e in zip(genes, ensembl_ids) if e in adata.var_names]
    genes, ensembl_ids = zip(*keep)

    test_mask = adata.obs["split"] == "test"
    sub = adata[test_mask, list(ensembl_ids)].to_memory()
    print(f"Subset: {sub.n_obs} held-out test cells x {sub.n_vars} genes")

    obs = sub.obs.copy()
    obs["state"] = obs["disease"].map(DISEASE_LABEL)
    obs["cd4cd8"] = obs["celltype"].map(CD4CD8_GROUP).fillna("other")
    n_counts = obs["n_counts"].to_numpy()

    X = sub.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsc()
    cp10k = X.multiply(1.0 / n_counts[:, None] * 1e4).tocsc()
    log1p_cp10k = cp10k.copy()
    log1p_cp10k.data = np.log1p(log1p_cp10k.data)

    rows = []
    group_cols = ["state", "individual", "cd4cd8"]
    for (state, donor, cd4cd8), idx in obs.groupby(group_cols, observed=True).indices.items():
        if len(idx) == 0:
            continue
        idx = np.asarray(idx)
        block = X[idx, :]
        block_norm = log1p_cp10k[idx, :]
        n_cells = len(idx)
        detect_rate = np.asarray((block > 0).sum(axis=0)).ravel() / n_cells
        mean_raw = np.asarray(block.mean(axis=0)).ravel()
        mean_log1p_cp10k = np.asarray(block_norm.mean(axis=0)).ravel()
        for j, gene in enumerate(genes):
            rows.append({
                "gene": gene, "program": gene_to_program[gene], "state": state,
                "donor": donor, "cd4cd8": cd4cd8, "n_cells": n_cells,
                "detect_rate": detect_rate[j], "mean_raw_count": mean_raw[j],
                "mean_log1p_cp10k": mean_log1p_cp10k[j],
            })
    per_donor = pd.DataFrame(rows)
    per_donor.to_csv(OUT_DIR / "baseline_expression_per_donor.csv", index=False)

    # State-level rollup pooling all held-out test cells of that state directly
    # (not an average of donor means), since donor cell counts are wildly
    # unequal (Normal is a single donor).
    pooled_rows = []
    for state, idx in obs.groupby("state", observed=True).indices.items():
        idx = np.asarray(idx)
        block = X[idx, :]
        block_norm = log1p_cp10k[idx, :]
        n_cells = len(idx)
        detect_rate = np.asarray((block > 0).sum(axis=0)).ravel() / n_cells
        mean_raw = np.asarray(block.mean(axis=0)).ravel()
        mean_log1p_cp10k = np.asarray(block_norm.mean(axis=0)).ravel()
        n_donors = obs.loc[obs["state"] == state, "individual"].nunique()
        for j, gene in enumerate(genes):
            pooled_rows.append({
                "gene": gene, "program": gene_to_program[gene], "state": state,
                "n_cells": n_cells, "n_donors": n_donors,
                "detect_rate": detect_rate[j], "mean_raw_count": mean_raw[j],
                "mean_log1p_cp10k": mean_log1p_cp10k[j],
            })
    pooled = pd.DataFrame(pooled_rows)
    pooled.to_csv(OUT_DIR / "baseline_expression_pooled.csv", index=False)

    print(f"\nWrote {OUT_DIR / 'baseline_expression_per_donor.csv'} ({len(per_donor)} rows)")
    print(f"Wrote {OUT_DIR / 'baseline_expression_pooled.csv'} ({len(pooled)} rows)")
    print("\nPooled detection rate, exhaustion program:")
    exh = pooled[pooled["program"] == "exhaustion"].pivot(index="gene", columns="state", values="detect_rate")
    print(exh[["normal", "sclc", "luad"]].round(3).to_string())


if __name__ == "__main__":
    main()
