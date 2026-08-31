#!/usr/bin/env python3
"""Build the deterministic T4 program, titration, and matched-null manifest."""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
EXPRESSION = HERE / "results/state_pooled_expression_test_only.csv"
AMBIENT_TABLE = HERE.parent / "primary_test_perturbation/tables/ambient_risk_all_genes.csv"
AMBIENT_MANIFEST = HERE.parent / "primary_test_perturbation/tables/ambient_risk_manifest.json"
TARGET_PANEL = HERE.parent / "perturbation_workflow/targeted_panel/target_gene_panel.json"
OUTPUT = HERE / "t4_program_manifest.json"

SEED = 43
N_BINS = 20
N_NULL = 20
STATES = ("normal", "sclc", "luad")
PROGRAMS = {
    "exhaustion": ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"],
    "cytotoxicity": ["NKG7", "GNLY", "PRF1", "GZMB", "GZMH", "IFNG"],
    "progenitor": ["TCF7", "SLAMF6", "IL7R", "CCR7"],
    "sclc_subtype_tf": ["ASCL1", "NEUROD1", "POU2F3", "YAP1"],
}
EXCLUDE_PATTERN = re.compile(r"^(?:RPL|RPS|MRPL|MRPS|MT-|HBA|HBB$|HBD$|HSPA1)")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_gene_map(path: Path) -> dict[str, str]:
    payload = json.loads(path.read_text())
    mapping = {row["gene"]: row["ensembl_id"] for row in payload["genes"]}
    missing = set(gene for genes in PROGRAMS.values() for gene in genes) - set(mapping)
    if missing:
        raise ValueError(f"Target panel lacks T4 genes: {sorted(missing)}")
    return mapping


def pooled_expression(path: Path, n_bins: int) -> pd.DataFrame:
    frame = pd.read_csv(path)
    required = {"state", "n_cells", "gene", "gene_symbol", "detect_rate", "mean_log1p_cp10k"}
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Expression table is missing {sorted(missing)}")
    if set(frame["state"]) != set(STATES):
        raise ValueError(f"Expected states {STATES}, found {sorted(frame['state'].unique())}")
    frame["weighted_expression"] = frame["n_cells"] * frame["mean_log1p_cp10k"]
    frame["weighted_detection"] = frame["n_cells"] * frame["detect_rate"]
    pooled = (
        frame.groupby(["gene", "gene_symbol"], as_index=False)
        .agg(
            weighted_expression=("weighted_expression", "sum"),
            weighted_detection=("weighted_detection", "sum"),
            total_cells=("n_cells", "sum"),
        )
    )
    pooled["mean_log1p_cp10k"] = pooled["weighted_expression"] / pooled["total_cells"]
    pooled["detect_rate"] = pooled["weighted_detection"] / pooled["total_cells"]
    # A real symbol indicates this Ensembl ID was present in the ISP-derived
    # symbol map used by T5, which is the available tokenizability proxy here.
    pooled = pooled[~pooled["gene_symbol"].str.startswith("ENSG")].copy()
    pooled = pooled.drop_duplicates("gene_symbol", keep=False)
    ranks = pooled["mean_log1p_cp10k"].rank(method="first", pct=True)
    pooled["expression_bin"] = np.minimum((ranks * n_bins).astype(int), n_bins - 1)
    return pooled.sort_values("gene_symbol").reset_index(drop=True)


def null_pool(
    expression: pd.DataFrame,
    target_genes: set[str],
    ambient_table: Path,
    ambient_manifest: Path,
) -> tuple[pd.DataFrame, dict]:
    known_ambient: set[str] = set()
    threshold = np.nan
    if ambient_manifest.exists():
        payload = json.loads(ambient_manifest.read_text())
        known_ambient = set(payload.get("anchors_ambient", []))
        threshold = float(payload.get("flag_threshold", np.nan))

    pool = expression[~expression["gene_symbol"].isin(target_genes | known_ambient)].copy()
    pool = pool[~pool["gene_symbol"].str.match(EXCLUDE_PATTERN)]
    scored_excluded = 0
    if ambient_table.exists() and np.isfinite(threshold):
        ambient = pd.read_csv(ambient_table, usecols=["ensembl_id", "ambient_risk"])
        pool = pool.merge(ambient, left_on="gene", right_on="ensembl_id", how="left")
        flagged = pool["ambient_risk"].ge(threshold).fillna(False)
        scored_excluded = int(flagged.sum())
        pool = pool[~flagged].drop(columns=["ensembl_id"])
    return pool.reset_index(drop=True), {
        "known_ambient_excluded": len(known_ambient),
        "high_ambient_risk_excluded": scored_excluded,
        "ambient_risk_threshold": threshold if np.isfinite(threshold) else None,
        "pattern_exclusion": EXCLUDE_PATTERN.pattern,
    }


def choose_match(
    target: pd.Series,
    pool: pd.DataFrame,
    used: set[str],
    rng: np.random.Generator,
    n_bins: int,
) -> pd.Series:
    for radius in range(n_bins):
        allowed_bins = {target["expression_bin"] - radius, target["expression_bin"] + radius}
        candidates = pool[
            pool["expression_bin"].isin(allowed_bins) & ~pool["gene_symbol"].isin(used)
        ]
        if len(candidates):
            distances = (candidates["mean_log1p_cp10k"] - target["mean_log1p_cp10k"]).abs()
            closest = candidates.loc[distances.nsmallest(min(20, len(candidates))).index]
            return closest.iloc[int(rng.integers(0, len(closest)))]
    raise RuntimeError(f"No unused expression match for {target['gene_symbol']}")


def gene_records(symbols: list[str], gene_map: dict[str, str]) -> list[dict]:
    return [{"gene": gene, "ensembl_id": gene_map[gene]} for gene in symbols]


def display_path(path: Path) -> str:
    return os.path.relpath(path, HERE)


def build_manifest(
    expression_path: Path = EXPRESSION,
    ambient_table: Path = AMBIENT_TABLE,
    ambient_manifest: Path = AMBIENT_MANIFEST,
    target_panel: Path = TARGET_PANEL,
    seed: int = SEED,
    n_bins: int = N_BINS,
    n_null: int = N_NULL,
) -> dict:
    rng = np.random.default_rng(seed)
    expression = pooled_expression(expression_path, n_bins)
    expression_by_symbol = expression.set_index("gene_symbol")
    gene_map = load_gene_map(target_panel)
    target_genes = {gene for genes in PROGRAMS.values() for gene in genes}
    missing_expression = target_genes - set(expression_by_symbol.index)
    if missing_expression:
        raise ValueError(f"No expression row for T4 genes: {sorted(missing_expression)}")
    pool, exclusions = null_pool(expression, target_genes, ambient_table, ambient_manifest)

    runs = []
    for program, symbols in PROGRAMS.items():
        runs.append(
            {
                "id": f"program_{program}",
                "phase": "program",
                "program": program,
                "genes": gene_records(symbols, gene_map),
            }
        )

    exhaustion_order = sorted(
        PROGRAMS["exhaustion"],
        key=lambda gene: (-expression_by_symbol.loc[gene, "mean_log1p_cp10k"], gene),
    )
    for size in range(1, len(exhaustion_order) + 1):
        symbols = exhaustion_order[:size]
        runs.append(
            {
                "id": f"exhaustion_nested_{size:02d}",
                "phase": "nested",
                "program": "exhaustion",
                "nested_size": size,
                "genes": gene_records(symbols, gene_map),
            }
        )

    null_audit = []
    pool_by_symbol = pool.set_index("gene_symbol")
    for program, target_symbols in PROGRAMS.items():
        for iteration in range(1, n_null + 1):
            used: set[str] = set()
            matches = []
            for target_symbol in target_symbols:
                target = expression_by_symbol.loc[target_symbol]
                match = choose_match(target, pool, used, rng, n_bins)
                used.add(str(match["gene_symbol"]))
                matches.append({"gene": str(match["gene_symbol"]), "ensembl_id": str(match["gene"])})
                null_audit.append(
                    {
                        "program": program,
                        "iteration": iteration,
                        "target_gene": target_symbol,
                        "target_expression_bin": int(target["expression_bin"]),
                        "target_mean_log1p_cp10k": float(target["mean_log1p_cp10k"]),
                        "matched_gene": str(match["gene_symbol"]),
                        "matched_expression_bin": int(match["expression_bin"]),
                        "matched_mean_log1p_cp10k": float(match["mean_log1p_cp10k"]),
                    }
                )
            runs.append(
                {
                    "id": f"null_{program}_{iteration:02d}",
                    "phase": "null",
                    "program": program,
                    "null_iteration": iteration,
                    "genes": matches,
                }
            )

    return {
        "schema_version": 1,
        "design": {
            "perturb_type": "overexpress",
            "sources": list(STATES),
            "random_seed": seed,
            "expression_bins": n_bins,
            "matched_null_sets_per_program": n_null,
            "empirical_one_sided_p_minimum": 1 / (n_null + 1),
            "expression_profile": "cell-count-weighted held-out-test mean_log1p_cp10k across three states",
            "null_reuse_across_sources": True,
            "exhaustion_nested_order": exhaustion_order,
            "null_pool_size": len(pool),
            "null_pool_exclusions": exclusions,
        },
        "inputs": {
            "expression": {"path": display_path(expression_path), "sha256": sha256(expression_path)},
            "target_panel": {"path": display_path(target_panel), "sha256": sha256(target_panel)},
            "ambient_table": {"path": display_path(ambient_table), "sha256": sha256(ambient_table)},
            "ambient_manifest": {"path": display_path(ambient_manifest), "sha256": sha256(ambient_manifest)},
        },
        "run_counts": {
            "program_definitions": len(PROGRAMS),
            "nested_definitions": len(exhaustion_order),
            "null_definitions": len(PROGRAMS) * n_null,
            "gpu_runs_all_sources": (len(PROGRAMS) + len(exhaustion_order) + len(PROGRAMS) * n_null) * len(STATES),
        },
        "runs": runs,
        "null_match_audit": null_audit,
        "null_pool_preview": pool_by_symbol.reset_index()[
            ["gene_symbol", "gene", "expression_bin", "mean_log1p_cp10k", "detect_rate"]
        ].head(20).to_dict(orient="records"),
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--n-null", type=int, default=N_NULL)
    parser.add_argument("--seed", type=int, default=SEED)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = build_manifest(seed=args.seed, n_null=args.n_null)
    args.output.write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest["run_counts"], indent=2))
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
