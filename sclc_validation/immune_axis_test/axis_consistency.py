#!/usr/bin/env python3
"""Test whether the poster's "Normal < SCLC < LUAD checkpoint axis" is a
self-consistent one-dimensional ordering, using only the overexpression arm
of the targeted 50-gene panel.

The deletion arm cannot answer this: it is detection-limited, so a different
subset of genes is testable in each source state and power varies by three
orders of magnitude between genes. The overexpression arm is census-complete
-- every held-out cell of the source is perturbed for every gene -- so the
50 genes x 3 sources x 2 goals design is balanced by construction and the
three source states are directly comparable.

Three tests, all computable from the committed stats CSVs:

  T1a collinearity  If the three states lie on one line in the order
                    Normal < SCLC < LUAD, then a gene that raises the
                    coordinate must, from Normal, move the cell toward both
                    SCLC and LUAD (same sign); from LUAD, away from both
                    (same sign); and from SCLC, toward one and away from the
                    other (opposite sign). All three must hold together.

  T1b reciprocity   A perturbation that displaces cells along a fixed
                    direction must raise similarity to one state exactly as
                    it lowers similarity to the other: shift(A->B) and
                    shift(B->A) must carry opposite signs. Same-sign pairs
                    mean the displacement is largely orthogonal to the A-B
                    contrast, so neither number is evidence about where A and
                    B sit relative to each other.

  T1c complement    If the shift toward the second goal is predictable from
                    the shift toward the first by a single gene-independent
                    slope, it is a geometric consequence of the first, not an
                    independent measurement. High R^2 means the "away from
                    Normal" term carries no information of its own.

Writes results/ and prints the summary. No GPU, no raw pickles -- runs on the
laptop from the tracked result tables.
"""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PANEL_DIR = HERE.parent / "perturbation_workflow" / "targeted_panel"
MERGED = PANEL_DIR / "results" / "targeted_panel_delete_overexpress_merged.csv"
PANEL_JSON = PANEL_DIR / "target_gene_panel.json"
OUT_DIR = HERE / "results"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]

# Reciprocal pairs: the same two states, scored from each side.
RECIPROCAL = [
    ("sclc_to_luad", "luad_to_sclc"),
    ("sclc_to_normal", "normal_to_sclc"),
    ("luad_to_normal", "normal_to_luad"),
]

# Source state -> (first goal, second goal). The complement test asks how much
# of the second is predictable from the first.
COMPLEMENT = [
    ("sclc", "sclc_to_luad", "sclc_to_normal"),
    ("normal", "normal_to_luad", "normal_to_sclc"),
    ("luad", "luad_to_normal", "luad_to_sclc"),
]

# The 21-gene pre-registered panel splits into four disjoint programs. The
# axis claim is about the exhaustion program specifically, so every test is
# reported for it separately as well as for all 50 genes.
PROGRAMS = {
    "exhaustion": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"],
    "cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG"],
    "progenitor": ["TCF7", "SLAMF6", "IL7R", "CCR7"],
    "sclc_subtype_tf": ["ASCL1", "NEUROD1", "POU2F3", "YAP1"],
}


def load_shifts() -> pd.DataFrame:
    merged = pd.read_csv(MERGED)
    wide = merged.pivot_table(
        index="Gene_name", columns="comparison", values="overexpress_shift"
    )[COMPARISONS]
    assert not wide.isna().any().any(), "overexpression arm should be complete for all 50 genes"
    return wide


def gene_groups(index: pd.Index) -> pd.DataFrame:
    meta = {g["gene"]: g["source"] for g in json.loads(PANEL_JSON.read_text())["genes"]}
    program = {gene: name for name, genes in PROGRAMS.items() for gene in genes}
    return pd.DataFrame(
        {
            "gene_set": [meta.get(g, "unknown") for g in index],
            "program": [program.get(g, "none") for g in index],
        },
        index=index,
    )


def collinearity(wide: pd.DataFrame) -> pd.DataFrame:
    """Per-gene truth table for the three sign predictions of a 1-D ordering."""
    sign = np.sign(wide)
    out = pd.DataFrame(index=wide.index)
    # From Normal, both other states lie on the same side of the axis.
    out["from_normal_same_sign"] = sign["normal_to_sclc"] == sign["normal_to_luad"]
    # From LUAD, likewise -- it is the far end.
    out["from_luad_same_sign"] = sign["luad_to_sclc"] == sign["luad_to_normal"]
    # From SCLC, the two goals lie on opposite sides -- it is the middle.
    out["from_sclc_opposite_sign"] = sign["sclc_to_luad"] != sign["sclc_to_normal"]
    out["consistent_1d"] = out.all(axis=1)
    return out


def reciprocity(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for a, b in RECIPROCAL:
        rows.append(
            {
                "pair": f"{a} / {b}",
                "frac_opposite_sign": float((np.sign(wide[a]) != np.sign(wide[b])).mean()),
                "pearson_r": float(np.corrcoef(wide[a], wide[b])[0, 1]),
                "n_genes": int(len(wide)),
            }
        )
    return pd.DataFrame(rows)


def complement(wide: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for source, first, second in COMPLEMENT:
        slope, intercept = np.polyfit(wide[first], wide[second], 1)
        r = float(np.corrcoef(wide[first], wide[second])[0, 1])
        rows.append(
            {
                "source_state": source,
                "predictor": first,
                "response": second,
                "slope": float(slope),
                "intercept": float(intercept),
                "r_squared": r ** 2,
            }
        )
    return pd.DataFrame(rows)


def response_dimensionality(wide: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    """How many independent directions does the 6-comparison response occupy?

    A genuine 1-D axis would put nearly all variance on one component whose
    loadings are antisymmetric across reciprocal pairs.
    """
    centered = wide.values - wide.values.mean(axis=0)
    _, singular, right = np.linalg.svd(centered, full_matrices=False)
    variance = singular ** 2 / (singular ** 2).sum()
    summary = pd.DataFrame(
        {"component": np.arange(1, len(variance) + 1), "variance_explained": variance}
    )
    loadings = pd.DataFrame(
        right[:3].T, index=wide.columns, columns=["PC1", "PC2", "PC3"]
    )
    return summary, loadings


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    wide = load_shifts()
    groups = gene_groups(wide.index)

    signs = collinearity(wide).join(groups)
    signs.to_csv(OUT_DIR / "axis_collinearity_by_gene.csv")

    recip = reciprocity(wide)
    recip.to_csv(OUT_DIR / "axis_reciprocity.csv", index=False)

    comp = complement(wide)
    comp.to_csv(OUT_DIR / "axis_goal_complement_regression.csv", index=False)

    variance, loadings = response_dimensionality(wide)
    variance.to_csv(OUT_DIR / "axis_response_variance.csv", index=False)
    loadings.to_csv(OUT_DIR / "axis_response_loadings.csv")

    checks = ["from_normal_same_sign", "from_luad_same_sign", "from_sclc_opposite_sign", "consistent_1d"]
    by_program = signs.groupby("program")[checks].mean()
    by_program.loc["all_50_genes"] = signs[checks].mean()
    by_program.to_csv(OUT_DIR / "axis_collinearity_by_program.csv")

    print("T1a  Collinearity -- fraction of genes satisfying each 1-D prediction")
    print(by_program.round(3).to_string())
    print(f"\n     Genes satisfying all three: {int(signs['consistent_1d'].sum())} / {len(signs)}")

    print("\nT1b  Reciprocity -- fixed-direction displacement predicts frac_opposite_sign = 1.0")
    print(recip.round(3).to_string(index=False))

    print("\nT1c  Second goal as a geometric complement of the first")
    print(comp.round(3).to_string(index=False))

    print("\nT1d  Dimensionality of the 6-comparison response")
    print(variance.round(3).to_string(index=False))
    print()
    print(loadings.round(3).to_string())

    print(f"\nWrote 6 tables to {OUT_DIR.relative_to(HERE.parent.parent)}/")


if __name__ == "__main__":
    main()
