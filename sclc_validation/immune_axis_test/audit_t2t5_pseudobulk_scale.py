#!/usr/bin/env python3
"""Audit the T2/T5 pseudobulk normalization scale without changing outputs.

T2's ``measure_baseline_expression.py`` divides counts by ``obs['n_counts']``.
T5's ``differential_expression.py`` uses Scanpy ``normalize_total``, which divides
by the row sum of ``X``.  These are equivalent only if the stored ``n_counts``
equals the current matrix row sum.  This audit measures that assumption on the
same held-out cells and seven exhaustion genes, then reproduces both published
state means from the raw object before making any source or table change.

It writes a compact JSON report only.  It never overwrites T2/T5 artifacts.
"""
from __future__ import annotations

import json
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

LAB_ROOT = Path(os.environ.get("LAB_ROOT", "/srv/lab"))
H5AD = Path(os.environ.get(
    "HTAN_H5AD",
    LAB_ROOT / "KD/sclc_luad_normal_htan_finetune/data/htan_sclc_luad_normal_tcells_prepared.h5ad",
))
STATS_ROOT = Path(os.environ.get(
    "SCLC_PERTURBATION_ROOT",
    LAB_ROOT / "KD/sclc_luad_normal_htan_heldout_allgene_perturbation",
)) / "stats"
OUT = Path(os.environ.get(
    "T2T5_SCALE_AUDIT_OUT",
    Path(__file__).resolve().parent / "results" / "t2t5_pseudobulk_scale_audit.json",
))

DISEASE_LABEL = {
    "normal": "normal",
    "small cell lung carcinoma": "sclc",
    "lung adenocarcinoma": "luad",
}
EXHAUSTION = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]


def symbol_to_ensembl(stats_root: Path = STATS_ROOT) -> dict[str, str]:
    """Read the same symbol map T2 uses."""
    mapping: dict[str, str] = {}
    for path in sorted(stats_root.glob("*/heldout_allgene_*.csv")):
        frame = pd.read_csv(path, usecols=["Gene_name", "Ensembl_ID"])
        mapping.update(dict(zip(frame.Gene_name, frame.Ensembl_ID)))
    if not mapping:
        raise SystemExit(f"No stats tables found under {stats_root}")
    return mapping


def as_csr(matrix) -> sp.csr_matrix:
    return matrix.tocsr() if sp.issparse(matrix) else sp.csr_matrix(matrix)


def log1p_cp10k(matrix: sp.csr_matrix, denominator: np.ndarray) -> sp.csr_matrix:
    """T2/T5 shared transform, parameterized only by the library denominator."""
    if np.any(denominator <= 0):
        raise ValueError("Library denominators must be positive")
    scaled = matrix.multiply(1e4 / denominator[:, None]).tocsr()
    scaled.data = np.log1p(scaled.data)
    return scaled


def state_program_means(norm: sp.csr_matrix, state: pd.Series) -> dict[str, float]:
    """Cell-weighted mean over the seven exhaustion genes, by disease state."""
    values: dict[str, float] = {}
    for label in ("normal", "sclc", "luad"):
        idx = np.flatnonzero(state.to_numpy() == label)
        if not len(idx):
            raise ValueError(f"No held-out cells for {label}")
        values[label] = float(norm[idx, :].mean())
    return values


def summarize_ratio(n_counts: np.ndarray, row_sums: np.ndarray) -> dict[str, float]:
    """Describe current-matrix counts relative to the stored count metadata."""
    ratio = row_sums / n_counts
    return {
        "min": float(ratio.min()),
        "median": float(np.median(ratio)),
        "mean": float(ratio.mean()),
        "max": float(ratio.max()),
        "fraction_exactly_one": float(np.mean(np.isclose(ratio, 1.0, rtol=0, atol=1e-10))),
    }


def main() -> None:
    adata = ad.read_h5ad(H5AD, backed="r")
    mask = adata.obs["split"] == "test"
    symbols = symbol_to_ensembl()
    missing = [gene for gene in EXHAUSTION if symbols.get(gene) not in adata.var_names]
    if missing:
        raise SystemExit(f"Exhaustion genes missing from H5AD: {missing}")
    ensembl = [symbols[gene] for gene in EXHAUSTION]

    # Keep a full-gene view for the T5 denominator and a seven-gene view for
    # the output score.  T2's numerator is the latter but its denominator is
    # `obs['n_counts']`; T5's denominator is the former row sum.
    full = adata[mask].to_memory()
    genes = full[:, ensembl].copy()
    full_x = as_csr(full.X)
    gene_x = as_csr(genes.X)
    stored = full.obs["n_counts"].to_numpy(dtype=float)
    row_sums = np.asarray(full_x.sum(axis=1)).ravel().astype(float)
    state = full.obs["disease"].map(DISEASE_LABEL)
    if state.isna().any():
        raise SystemExit("Unmapped disease labels in held-out cells")

    t2_norm = log1p_cp10k(gene_x, stored)
    t5_norm = log1p_cp10k(gene_x, row_sums)
    t2_means = state_program_means(t2_norm, state)
    t5_means = state_program_means(t5_norm, state)
    ratios = {key: t5_means[key] / t2_means[key] for key in t2_means}
    count_ratio = summarize_ratio(stored, row_sums)

    report = {
        "h5ad": str(H5AD),
        "population": "test_only",
        "n_cells": int(full.n_obs),
        "genes": EXHAUSTION,
        "t2_denominator": "obs['n_counts']",
        "t5_denominator": "row sum of X (Scanpy normalize_total default)",
        "current_matrix_over_stored_n_counts": count_ratio,
        "t2_formula_state_means": t2_means,
        "t5_formula_state_means": t5_means,
        "t5_over_t2_state_mean_ratio": ratios,
        "root_cause_supported": not np.isclose(count_ratio["median"], 1.0, rtol=0, atol=1e-10),
        "interpretation": (
            "A non-unit current_matrix_over_stored_n_counts ratio proves the two pipelines use "
            "different library-size denominators on otherwise identical cells and genes. "
            "No result tables are changed by this audit."
        ),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
