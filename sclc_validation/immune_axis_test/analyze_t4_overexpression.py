#!/usr/bin/env python3
"""Analyze T4 raw shifts by donor/subtype, matched null, and nested size."""
from __future__ import annotations

import argparse
import json
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

import run_t4_overexpression as runner

HERE = Path(__file__).resolve().parent
OUTPUT_DIR = HERE / "results"
CD4CD8_GROUP = {
    "CD4-positive helper T cell": "CD4",
    "CD8-positive, alpha-beta memory T cell": "CD8",
    "effector CD8-positive, alpha-beta T cell": "CD8",
    "regulatory T cell": "CD4 (Treg)",
}


def load_raw_shifts(path: Path) -> dict[str, np.ndarray]:
    # Internal experiment pickles only; never load an untrusted path.
    with path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected raw shift dict in {path}")
    shifts = {}
    for state_name, inner in payload.items():
        if state_name not in runner.STATE_SLUGS:
            continue
        if not isinstance(inner, dict) or len(inner) != 1:
            raise ValueError(f"Expected one group-perturbation key for {state_name} in {path}")
        values = np.asarray(next(iter(inner.values())), dtype=float)
        if values.ndim != 1 or not np.isfinite(values).all():
            raise ValueError(f"Invalid shifts for {state_name} in {path}: shape={values.shape}")
        shifts[runner.STATE_SLUGS[state_name]] = values
    missing = set(runner.STATE_NAMES) - set(shifts)
    if missing:
        raise ValueError(f"Raw shift file lacks states {sorted(missing)}: {path}")
    lengths = {len(values) for values in shifts.values()}
    if len(lengths) != 1:
        raise ValueError(f"State shift lengths differ in {path}: {sorted(lengths)}")
    return shifts


def ordered_metadata(source_path: Path, load_from_disk) -> pd.DataFrame:
    dataset = load_from_disk(str(source_path))
    required = {"individual", "celltype", "length"}
    missing = required - set(dataset.column_names)
    if missing:
        raise ValueError(f"Source dataset is missing metadata columns: {sorted(missing)}")
    frame = pd.DataFrame({column: dataset[column] for column in required})
    frame = frame.sort_values("length", ascending=False, kind="stable").reset_index(drop=True)
    frame["cd4cd8"] = frame["celltype"].map(CD4CD8_GROUP).fillna("other")
    return frame


def scalar_summary(values: np.ndarray, donors: pd.Series) -> dict:
    by_donor = pd.DataFrame({"donor": donors.to_numpy(), "shift": values}).groupby("donor")["shift"].mean()
    mean = float(np.mean(values))
    donor_mean = float(by_donor.mean())
    overall_sign = np.sign(mean)
    same_sign = (np.sign(by_donor) == overall_sign) | (by_donor == 0)
    return {
        "n_cells": len(values),
        "n_donors": len(by_donor),
        "mean_shift": mean,
        "median_shift": float(np.median(values)),
        "sd_shift": float(np.std(values, ddof=1)) if len(values) > 1 else np.nan,
        "cell_se": float(np.std(values, ddof=1) / np.sqrt(len(values))) if len(values) > 1 else np.nan,
        "mean_donor_shift": donor_mean,
        "donor_se": float(by_donor.std(ddof=1) / np.sqrt(len(by_donor))) if len(by_donor) > 1 else np.nan,
        "donor_shift_min": float(by_donor.min()),
        "donor_shift_max": float(by_donor.max()),
        "frac_donors_same_sign": float(same_sign.mean()),
    }


def summarize_run(
    run: dict,
    source: str,
    shifts: dict[str, np.ndarray],
    metadata: pd.DataFrame,
) -> tuple[list[dict], list[dict]]:
    n_values = len(next(iter(shifts.values())))
    if len(metadata) != n_values:
        raise ValueError(f"Metadata/shift mismatch for {run['id']}/{source}: {len(metadata)} != {n_values}")
    summary_rows, donor_rows = [], []
    scopes = [("all", pd.Series(True, index=metadata.index))]
    scopes += [(label, metadata["cd4cd8"].eq(label)) for label in sorted(metadata["cd4cd8"].unique())]
    for target, values in shifts.items():
        for scope, mask in scopes:
            selected = values[mask.to_numpy()]
            if not len(selected):
                continue
            row = {
                "item_id": run["id"],
                "phase": run["phase"],
                "program": run["program"],
                "nested_size": run.get("nested_size"),
                "null_iteration": run.get("null_iteration"),
                "source_state": source,
                "target_state": target,
                "scope": scope,
            }
            row.update(scalar_summary(selected, metadata.loc[mask, "individual"]))
            summary_rows.append(row)

            donor_frame = pd.DataFrame(
                {
                    "donor": metadata.loc[mask, "individual"].to_numpy(),
                    "shift": selected,
                }
            )
            for donor, donor_values in donor_frame.groupby("donor")["shift"]:
                donor_rows.append(
                    {
                        "item_id": run["id"],
                        "phase": run["phase"],
                        "program": run["program"],
                        "nested_size": run.get("nested_size"),
                        "null_iteration": run.get("null_iteration"),
                        "source_state": source,
                        "target_state": target,
                        "scope": scope,
                        "donor": donor,
                        "n_cells": len(donor_values),
                        "mean_shift": float(donor_values.mean()),
                    }
                )
    return summary_rows, donor_rows


def matched_null_tests(summary: pd.DataFrame, expected_null: int) -> pd.DataFrame:
    overall = summary[summary["scope"] == "all"]
    primary = overall[overall["phase"] == "program"]
    null = overall[overall["phase"] == "null"]
    rows = []
    for observed in primary.itertuples():
        distribution = null[
            (null["program"] == observed.program)
            & (null["source_state"] == observed.source_state)
            & (null["target_state"] == observed.target_state)
        ]["mean_shift"].to_numpy()
        if not len(distribution):
            continue
        center = float(np.median(distribution))
        if observed.mean_shift >= center:
            directional_count = int(np.sum(distribution >= observed.mean_shift))
            direction = "greater"
        else:
            directional_count = int(np.sum(distribution <= observed.mean_shift))
            direction = "less"
        two_sided_count = int(np.sum(np.abs(distribution - center) >= abs(observed.mean_shift - center)))
        null_sd = float(np.std(distribution, ddof=1)) if len(distribution) > 1 else np.nan
        rows.append(
            {
                "program": observed.program,
                "source_state": observed.source_state,
                "target_state": observed.target_state,
                "observed_mean_shift": observed.mean_shift,
                "n_null": len(distribution),
                "expected_n_null": expected_null,
                "null_complete": len(distribution) == expected_null,
                "null_mean": float(np.mean(distribution)),
                "null_median": center,
                "null_sd": null_sd,
                "null_z": (observed.mean_shift - np.mean(distribution)) / null_sd if null_sd > 0 else np.nan,
                "directional_alternative": direction,
                "empirical_p_directional": (directional_count + 1) / (len(distribution) + 1),
                "empirical_p_two_sided": (two_sided_count + 1) / (len(distribution) + 1),
            }
        )
    return pd.DataFrame(rows)


def nested_monotonicity(summary: pd.DataFrame) -> pd.DataFrame:
    nested = summary[(summary["scope"] == "all") & (summary["phase"] == "nested")]
    rows = []
    for (source, target), group in nested.groupby(["source_state", "target_state"]):
        group = group.sort_values("nested_size")
        values = group["mean_shift"].to_numpy()
        differences = np.diff(values)
        rows.append(
            {
                "source_state": source,
                "target_state": target,
                "n_sizes": len(group),
                "spearman_size_vs_shift": float(group["nested_size"].corr(group["mean_shift"], method="spearman")),
                "n_non_decreasing_steps": int(np.sum(differences >= 0)),
                "n_non_increasing_steps": int(np.sum(differences <= 0)),
                "fully_non_decreasing": bool(np.all(differences >= 0)),
                "fully_non_increasing": bool(np.all(differences <= 0)),
                "shift_at_size_1": float(values[0]),
                "shift_at_max_size": float(values[-1]),
                "pre_registered_primary": source == "sclc" and target == "luad",
            }
        )
    return pd.DataFrame(rows)


def analyze(
    manifest: dict,
    run_dir: Path,
    phases: set[str],
    allow_partial: bool,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[str]]:
    from datasets import load_from_disk

    expected = [
        (run, source)
        for run in manifest["runs"]
        if run["phase"] in phases
        for source in runner.STATE_NAMES
    ]
    missing = [
        f"{run['id']}/{source}"
        for run, source in expected
        if not runner.completion_marker(run_dir, run, source).exists()
    ]
    if missing and not allow_partial:
        raise FileNotFoundError(f"T4 is incomplete: {len(missing)} missing markers; first={missing[:12]}")

    metadata_cache: dict[str, pd.DataFrame] = {}
    summary_rows, donor_rows = [], []
    for run, source in expected:
        marker = runner.completion_marker(run_dir, run, source)
        if not marker.exists():
            continue
        completion = json.loads(marker.read_text())
        raw_path = run_dir / completion["raw_file"]
        shifts = load_raw_shifts(raw_path)
        if source not in metadata_cache:
            source_path = run_dir / "data/sources" / f"{source}.dataset"
            metadata_cache[source] = ordered_metadata(source_path, load_from_disk)
        summaries, donors = summarize_run(run, source, shifts, metadata_cache[source])
        summary_rows.extend(summaries)
        donor_rows.extend(donors)

    summary = pd.DataFrame(summary_rows)
    donors = pd.DataFrame(donor_rows)
    null_tests = matched_null_tests(summary, manifest["design"]["matched_null_sets_per_program"]) if len(summary) else pd.DataFrame()
    monotonicity = nested_monotonicity(summary) if len(summary) else pd.DataFrame()
    return summary, donors, null_tests, monotonicity, missing


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=runner.MANIFEST)
    parser.add_argument("--run-dir", type=Path, default=runner.RUN_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--phase", choices=("all",) + runner.PHASES, default="all")
    parser.add_argument("--allow-partial", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = runner.load_manifest(args.manifest)
    phases = set(runner.PHASES if args.phase == "all" else [args.phase])
    summary, donors, null_tests, monotonicity, missing = analyze(
        manifest, args.run_dir, phases, args.allow_partial
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    summary.to_csv(args.output_dir / "t4_shift_summary.csv", index=False)
    donors.to_csv(args.output_dir / "t4_donor_summary.csv", index=False)
    null_tests.to_csv(args.output_dir / "t4_matched_null_test.csv", index=False)
    monotonicity.to_csv(args.output_dir / "t4_nested_monotonicity.csv", index=False)
    result = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "phases": sorted(phases),
        "allow_partial": args.allow_partial,
        "n_summary_rows": len(summary),
        "n_donor_rows": len(donors),
        "n_missing_runs": len(missing),
        "first_missing_runs": missing[:20],
        "limitations": [
            "cell-level standard errors are descriptive because cells within donors are not independent",
            "donor summaries are the biological-replicate view",
            "Normal has one held-out donor, so its donor uncertainty is not estimable",
            "20 null sets give a minimum empirical one-sided p-value of 1/21",
        ],
    }
    (args.output_dir / "t4_analysis_manifest.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
