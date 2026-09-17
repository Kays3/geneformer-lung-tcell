#!/usr/bin/env python3
"""T6 — donor-balanced robustness for the immune axis.

WHAT QUESTION THIS ANSWERS
--------------------------
T2 (test-only population) and T5a (complete population) disagree on the SCLC-vs-LUAD
ordering of the curated 7-gene exhaustion program:

    T2  test-only, pooled:  Normal 0.180 < SCLC 0.186 < LUAD 0.242   -> LUAD > SCLC
    T5a complete,  pooled:  Normal 0.123 < LUAD 0.260 < SCLC 0.286   -> SCLC > LUAD

RESULTS_T5.md attributes this to donor composition: SCLC's held-out test split is 74.9 %
one donor (`PleuralEffusion`), which sits near the bottom of the full 19-donor SCLC
distribution. This module tests that explanation directly by changing ONE thing - the
unit that carries weight - and holding the population fixed.

THE REPLICATION UNIT IS THE DONOR, NOT THE CELL.
Cells within a donor are not independent observations of "a state": they share a patient,
a sample, a dissociation and a sequencing run. Every pooled number above is CELL-weighted,
so a donor contributing 1,816 cells counts 1,816 times as much as a donor contributing
212. Donor-level (pseudobulk) scores give each donor weight 1. Both are computed here,
from the same rows, so the difference between them is attributable to weighting alone.

WHAT THIS DOES NOT DO
---------------------
It does not promote the axis claim. T4's matched null (directional empirical p = .4286)
and the failed strict titration leave the axis conclusion qualified, and nothing here
changes that. This is reconciliation and robustness only: no retuning, no ranking rerun,
no new claim.

INPUTS - committed tables only, no raw per-cell data, CPU only
--------------------------------------------------------------
  results/baseline_expression_per_donor.csv
      T2's own authoritative input (test-only population, 8 donors: 4 LUAD, 3 SCLC,
      1 Normal). Reproduces t2_program_summary.csv exactly. Carries the CD4/CD8 strata,
      which exist for the test-only population ONLY.
  results/t5a_donor_composition_check.csv
      T5a's per-donor exhaustion score for the complete population (45 donors: 22 LUAD,
      19 SCLC, 4 Normal).

Both inputs use CP10k normalization from the current full-`X` row sum. T2 and T5 must
therefore agree when evaluated on the same held-out cells and seven genes; the regression
test records that invariant. This script still compares weighting schemes *within* each
population rather than treating distinct patient populations as interchangeable.
"""
from __future__ import annotations

import itertools
import json
from math import comb
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"

EXHAUSTION = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]
PERMUTATIONS = 100_000
SEED = 20260917


# ----------------------------------------------------------------- donor scores

def test_only_donor_scores() -> pd.DataFrame:
    """Per-donor exhaustion score for the test-only population, from T2's own table.

    A donor's score is the mean over the seven genes of that gene's cell-weighted mean
    across the donor's CD4/CD8 strata - i.e. the donor's pseudobulk value. Weighting by
    n_cells WITHIN a donor is correct (it reconstructs that donor's cell mean); weighting
    ACROSS donors is the thing under test and is done separately.
    """
    frame = pd.read_csv(RESULTS / "baseline_expression_per_donor.csv")
    frame = frame[frame["program"] == "exhaustion"]
    missing = set(EXHAUSTION) - set(frame["gene"])
    if missing:
        raise SystemExit(f"exhaustion genes absent from baseline table: {sorted(missing)}")

    per_gene = (
        frame.groupby(["state", "donor", "gene"])
        .apply(lambda d: pd.Series({
            "value": np.average(d["mean_log1p_cp10k"], weights=d["n_cells"]),
            "n_cells": d["n_cells"].sum(),
        }), include_groups=False)
        .reset_index()
    )
    donors = per_gene.groupby(["state", "donor"]).agg(
        score=("value", "mean"), n_cells=("n_cells", "mean")
    ).reset_index()
    donors["population"] = "test_only"
    return donors


def complete_donor_scores() -> pd.DataFrame:
    """Per-donor exhaustion score for the complete population, from T5a's own table."""
    frame = pd.read_csv(RESULTS / "t5a_donor_composition_check.csv")
    donors = frame.rename(columns={"mean_log1p_cp10k": "score"})[
        ["state", "donor", "score", "n_cells"]
    ].copy()
    donors["population"] = "complete"
    return donors


# ------------------------------------------------------------------ weighting

def weighting_table(donors: pd.DataFrame) -> pd.DataFrame:
    """Cell-weighted vs donor-level state means, from identical rows."""
    rows = []
    for (population, state), group in donors.groupby(["population", "state"]):
        rows.append({
            "population": population,
            "state": state,
            "n_donors": len(group),
            "n_cells": int(group["n_cells"].sum()),
            "cell_weighted_mean": float(np.average(group["score"], weights=group["n_cells"])),
            "donor_level_mean": float(group["score"].mean()),
            "donor_level_sd": float(group["score"].std(ddof=1)) if len(group) > 1 else float("nan"),
            "max_donor_cell_share": float(group["n_cells"].max() / group["n_cells"].sum()),
        })
    return pd.DataFrame(rows).sort_values(["population", "state"]).reset_index(drop=True)


# ---------------------------------------------------------------- permutation

def permutation_floor(n_a: int, n_b: int) -> dict:
    """Exact floor on an attainable two-sided p, BEFORE any p-value is interpreted.

    Under label permutation of donors, there are C(n_a+n_b, n_a) distinct assignments.
    The observed assignment and its mirror are always at least as extreme as itself, so
    no two-sided p below 2/C(n_a+n_b, n_a) exists. A criterion stricter than this floor
    is a check that cannot pass, whatever the data say. For a balanced split n_a=n_b=k
    this is the 2/C(2k,k) form.
    """
    total = comb(n_a + n_b, n_a)
    return {
        "n_a": n_a, "n_b": n_b, "n_assignments": total,
        "min_two_sided_p": 2 / total,
        "min_one_sided_p": 1 / total,
        "two_sided_alpha_05_attainable": bool(2 / total <= 0.05),
    }


def permutation_test(donors: pd.DataFrame, state_a: str, state_b: str,
                     label: str, rng: np.random.Generator) -> dict:
    """Donor-label permutation test on the difference of donor-level means.

    NULL IN FORCE: the donor-level exhaustion score is exchangeable between the two state
    labels - i.e. a donor's score carries no information about whether that donor is
    SCLC or LUAD. The alternative is two-sided (no direction is assumed).

    Exact enumeration when the assignment count is small enough; otherwise Monte Carlo,
    whose resolution is reported alongside so a p-value is never read finer than the
    method can express.
    """
    a = donors.loc[donors["state"] == state_a, "score"].to_numpy(float)
    b = donors.loc[donors["state"] == state_b, "score"].to_numpy(float)
    observed = a.mean() - b.mean()
    pooled = np.concatenate([a, b])
    n_a = len(a)
    floor = permutation_floor(n_a, len(b))

    if floor["n_assignments"] <= 200_000:
        stats = np.array([
            pooled[list(idx)].mean() - pooled[list(set(range(len(pooled))) - set(idx))].mean()
            for idx in itertools.combinations(range(len(pooled)), n_a)
        ])
        p_two = float(np.mean(np.abs(stats) >= abs(observed) - 1e-12))
        method, resolution = "exact_enumeration", 1 / floor["n_assignments"]
    else:
        count = 0
        for _ in range(PERMUTATIONS):
            shuffled = rng.permutation(pooled)
            stat = shuffled[:n_a].mean() - shuffled[n_a:].mean()
            if abs(stat) >= abs(observed) - 1e-12:
                count += 1
        p_two = (1 + count) / (PERMUTATIONS + 1)
        method, resolution = f"monte_carlo_B={PERMUTATIONS}", 1 / (PERMUTATIONS + 1)

    return {
        "comparison": label, "state_a": state_a, "state_b": state_b,
        "n_donors_a": n_a, "n_donors_b": len(b),
        "observed_donor_level_difference": float(observed),
        "method": method, "p_two_sided": p_two, "p_resolution": resolution,
        "min_attainable_two_sided_p": floor["min_two_sided_p"],
        "n_assignments": floor["n_assignments"],
        "two_sided_alpha_05_attainable": floor["two_sided_alpha_05_attainable"],
        "null": "donor-level score exchangeable between the two state labels",
    }


# ----------------------------------------------------------------------- LODO

def leave_one_donor_out(donors: pd.DataFrame, state_a: str, state_b: str,
                        population: str) -> pd.DataFrame:
    """Drop each donor in turn; recompute both weightings. Sign stability is the point."""
    subset = donors[donors["state"].isin([state_a, state_b])]
    full_a, full_b = subset[subset.state == state_a], subset[subset.state == state_b]
    full_donor = full_a["score"].mean() - full_b["score"].mean()
    full_cell = (np.average(full_a["score"], weights=full_a["n_cells"])
                 - np.average(full_b["score"], weights=full_b["n_cells"]))
    rows = []
    for donor in subset["donor"].unique():
        kept = subset[subset["donor"] != donor]
        a, b = kept[kept.state == state_a], kept[kept.state == state_b]
        rows.append({
            "population": population,
            "dropped_donor": donor,
            "dropped_from_state": subset.loc[subset.donor == donor, "state"].iloc[0],
            "dropped_n_cells": int(subset.loc[subset.donor == donor, "n_cells"].iloc[0]),
            "n_donors_remaining": len(kept),
            "donor_level_difference": float(a["score"].mean() - b["score"].mean()),
            "cell_weighted_difference": float(
                np.average(a["score"], weights=a["n_cells"])
                - np.average(b["score"], weights=b["n_cells"])),
            "donor_level_sign_matches_full": bool(
                np.sign(a["score"].mean() - b["score"].mean()) == np.sign(full_donor)),
            "cell_weighted_sign_matches_full": bool(
                np.sign(np.average(a["score"], weights=a["n_cells"])
                        - np.average(b["score"], weights=b["n_cells"])) == np.sign(full_cell)),
        })
    return pd.DataFrame(rows).sort_values(["dropped_from_state", "dropped_donor"])


# ------------------------------------------------------------------- CD4/CD8

def cd4cd8_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """Composition shares per donor, and donor-level scores WITHIN each stratum.

    Available for the test-only population only - the complete population's per-cell
    CD4/CD8 assignment is not in any committed table. Stratifying answers the
    composition-confounding question directly: if SCLC > LUAD holds inside every
    stratum, the difference is not produced by the states' differing CD4/CD8 mix.
    """
    frame = pd.read_csv(RESULTS / "baseline_expression_per_donor.csv")
    frame = frame[frame["program"] == "exhaustion"]

    composition = frame.groupby(["state", "donor", "cd4cd8"])["n_cells"].first().reset_index()
    composition["share_of_donor_cells"] = (
        composition["n_cells"] / composition.groupby(["state", "donor"])["n_cells"].transform("sum"))

    stratified = (frame.groupby(["state", "donor", "cd4cd8"])["mean_log1p_cp10k"]
                  .mean().reset_index().rename(columns={"mean_log1p_cp10k": "score"}))
    return composition, stratified


def stratified_summary(stratified: pd.DataFrame, state_a: str, state_b: str) -> pd.DataFrame:
    rows = []
    for stratum, group in stratified.groupby("cd4cd8"):
        a = group.loc[group.state == state_a, "score"]
        b = group.loc[group.state == state_b, "score"]
        rows.append({
            "cd4cd8": stratum,
            f"donor_level_{state_a}": float(a.mean()),
            f"donor_level_{state_b}": float(b.mean()),
            "difference": float(a.mean() - b.mean()),
            "direction_matches_unstratified": bool(a.mean() > b.mean()),
            "n_donors_a": len(a), "n_donors_b": len(b),
        })
    return pd.DataFrame(rows)


# -------------------------------------------------------------------- driver

def main() -> int:
    rng = np.random.default_rng(SEED)
    donors = pd.concat([test_only_donor_scores(), complete_donor_scores()], ignore_index=True)
    donors.to_csv(RESULTS / "t6_donor_level_scores.csv", index=False)

    print("PERMUTATION FLOORS — computed BEFORE any p-value below is read")
    floors = {}
    for population, (n_a, n_b) in {
        "test_only": (3, 4), "complete": (19, 22)
    }.items():
        floor = permutation_floor(n_a, n_b)
        floors[population] = floor
        verdict = "attainable" if floor["two_sided_alpha_05_attainable"] else "NOT ATTAINABLE"
        print(f"  {population:10s} SCLC {n_a} v LUAD {n_b}: C({n_a+n_b},{n_a})="
              f"{floor['n_assignments']:,}  min two-sided p = {floor['min_two_sided_p']:.4g}"
              f"  -> alpha=0.05 {verdict}")

    weighting = weighting_table(donors)
    weighting.to_csv(RESULTS / "t6_weighting_reconciliation.csv", index=False)
    print("\nWEIGHTING RECONCILIATION (identical rows, only the weight differs)")
    print(weighting.to_string(index=False))

    tests = [
        permutation_test(donors[donors.population == "test_only"], "sclc", "luad",
                         "test_only_sclc_vs_luad", rng),
        permutation_test(donors[donors.population == "complete"], "sclc", "luad",
                         "complete_sclc_vs_luad", rng),
    ]
    pd.DataFrame(tests).to_csv(RESULTS / "t6_permutation_tests.csv", index=False)
    print("\nPERMUTATION TESTS (null: donor-level score exchangeable between state labels)")
    for test in tests:
        print(f"  {test['comparison']:26s} diff={test['observed_donor_level_difference']:+.4f}"
              f"  p={test['p_two_sided']:.4g} ({test['method']}, floor {test['min_attainable_two_sided_p']:.4g})")

    lodo = pd.concat([
        leave_one_donor_out(donors[donors.population == "test_only"], "sclc", "luad", "test_only"),
        leave_one_donor_out(donors[donors.population == "complete"], "sclc", "luad", "complete"),
    ], ignore_index=True)
    lodo.to_csv(RESULTS / "t6_lodo.csv", index=False)
    print("\nLEAVE-ONE-DONOR-OUT, donor-level SCLC - LUAD difference")
    for population, group in lodo.groupby("population"):
        print(f"  {population:10s} n={len(group):2d} drops  range "
              f"[{group.donor_level_difference.min():+.4f}, {group.donor_level_difference.max():+.4f}]"
              f"  sign stable in {int(group.donor_level_sign_matches_full.sum())}/{len(group)}")
        cw = group["cell_weighted_difference"]
        print(f"             cell-weighted range [{cw.min():+.4f}, {cw.max():+.4f}]"
              f"  sign stable in {int(group.cell_weighted_sign_matches_full.sum())}/{len(group)}"
              f" (vs the FULL-data cell-weighted sign)")

    composition, stratified = cd4cd8_tables()
    composition.to_csv(RESULTS / "t6_cd4cd8_composition.csv", index=False)
    stratified.to_csv(RESULTS / "t6_cd4cd8_stratified_scores.csv", index=False)
    summary = stratified_summary(stratified, "sclc", "luad")
    summary.to_csv(RESULTS / "t6_cd4cd8_stratified_summary.csv", index=False)
    print("\nCD4/CD8 STRATIFIED (test-only population only — the complete population's\n"
          "strata are not in any committed table)")
    print(summary.to_string(index=False))

    manifest = {
        "analysis": "T6 donor-balanced robustness",
        "replication_unit": "donor",
        "permutation_floors": floors,
        "seed": SEED,
        "monte_carlo_permutations": PERMUTATIONS,
        "inputs": ["results/baseline_expression_per_donor.csv",
                   "results/t5a_donor_composition_check.csv"],
        "normalization_note": ("T2 and T5 use the current full-X row sum for CP10k "
                               "normalization. Matched test cells and genes must agree; "
                               "population comparisons remain limited by donor composition."),
        "claim_status": ("reconciliation and robustness only; the axis conclusion remains "
                         "qualified by T4's matched null and the failed strict titration"),
    }
    (RESULTS / "t6_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"\nWrote t6_* tables and manifest to {RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
