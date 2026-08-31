#!/usr/bin/env python3
"""Re-run T1 after excluding flagged ambient/stress/ribosomal genes."""
from __future__ import annotations

import json
import re
from pathlib import Path

import pandas as pd

import axis_consistency as axis

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
FLAGS = (
    ROOT
    / "sclc_validation/primary_test_perturbation/tables/hk_gene_diagnostic/targeted_panel_gene_flags.csv"
)
OUT = HERE / "results"
EXPLICIT_PATTERN = re.compile(r"^(?:HBA|HBB$|HBD$|HSPA1|RPL|RPS|MRPL|MRPS|S100A|SFTPA|EEF1)")


def flagged_genes(path: Path = FLAGS) -> tuple[set[str], pd.DataFrame]:
    frame = pd.read_csv(path)
    required = {"Gene_name", "is_hk", "ambient_flag"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Flag table is missing {sorted(missing)}")
    frame["explicit_class_flag"] = frame["Gene_name"].str.match(EXPLICIT_PATTERN)
    frame["excluded"] = frame[["is_hk", "ambient_flag", "explicit_class_flag"]].any(axis=1)
    excluded = set(frame.loc[frame["excluded"], "Gene_name"])
    return excluded, frame


def scenario_summary(name: str, wide: pd.DataFrame) -> tuple[dict, dict[str, pd.DataFrame]]:
    collinearity = axis.collinearity(wide)
    reciprocity = axis.reciprocity(wide)
    complement = axis.complement(wide)
    variance, loadings = axis.response_dimensionality(wide)
    row = {
        "scenario": name,
        "n_genes": len(wide),
        "n_consistent_1d": int(collinearity["consistent_1d"].sum()),
        "frac_consistent_1d": float(collinearity["consistent_1d"].mean()),
        "frac_from_normal_same_sign": float(collinearity["from_normal_same_sign"].mean()),
        "frac_from_luad_same_sign": float(collinearity["from_luad_same_sign"].mean()),
        "frac_from_sclc_opposite_sign": float(collinearity["from_sclc_opposite_sign"].mean()),
        "pc1_variance_explained": float(variance.iloc[0]["variance_explained"]),
        "pc2_variance_explained": float(variance.iloc[1]["variance_explained"]),
    }
    tables = {
        "collinearity": collinearity,
        "reciprocity": reciprocity,
        "complement": complement,
        "variance": variance,
        "loadings": loadings,
    }
    return row, tables


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    wide = axis.load_shifts()
    excluded, audit = flagged_genes()
    present_excluded = sorted(excluded & set(wide.index))
    clean = wide.drop(index=present_excluded)
    if len(clean) < 3:
        raise ValueError("Too few genes remain for sensitivity analysis")

    summaries = []
    all_tables = {}
    for name, frame in [("all_50_genes", wide), ("exclude_ambient_stress_ribosomal", clean)]:
        row, tables = scenario_summary(name, frame)
        summaries.append(row)
        all_tables[name] = tables
        if name != "all_50_genes":
            for table_name, table in tables.items():
                table.to_csv(OUT / f"axis_sensitivity_{name}_{table_name}.csv")

    summary = pd.DataFrame(summaries)
    summary.to_csv(OUT / "axis_sensitivity_summary.csv", index=False)
    audit.to_csv(OUT / "axis_sensitivity_exclusion_audit.csv", index=False)

    baseline = summary.set_index("scenario").loc["all_50_genes"]
    cleaned = summary.set_index("scenario").loc["exclude_ambient_stress_ribosomal"]
    result = {
        "status": "complete",
        "exclusion_pattern": EXPLICIT_PATTERN.pattern,
        "n_excluded": len(present_excluded),
        "excluded_genes": present_excluded,
        "n_remaining": len(clean),
        "consistent_1d_before": int(baseline["n_consistent_1d"]),
        "consistent_1d_after": int(cleaned["n_consistent_1d"]),
        "pc1_variance_before": float(baseline["pc1_variance_explained"]),
        "pc1_variance_after": float(cleaned["pc1_variance_explained"]),
        "conclusion_reversed": bool(
            cleaned["frac_consistent_1d"] >= 0.5
            and baseline["frac_consistent_1d"] < 0.5
        ),
        "decision_rule": "reversal requires at least half of remaining genes to satisfy all three 1-D sign predictions",
    }
    (OUT / "axis_sensitivity_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(summary.round(3).to_string(index=False))
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
