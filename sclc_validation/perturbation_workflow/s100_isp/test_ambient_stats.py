#!/usr/bin/env python3
"""Tests for ambient_stats.py and matched_controls.py's stub, all against
synthetic data -- no GPU, no real dataset, no geneformer import.

Run: python3 test_ambient_stats.py
"""
from __future__ import annotations

import itertools
import sys
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, ".")
import ambient_stats as astats
import matched_controls as mc
import retained_rows as rr


def test_compute_E():
    rows = pd.DataFrame([
        {"gene": "G1", "ensembl_id": "E1", "role": "r", "stratum": "s1",
         "source_state": "sclc", "goal_state": "luad", "alt_state": "normal",
         "perturbation_type": "delete", "donor_balanced_shift": 0.4},
        {"gene": "G1", "ensembl_id": "E1", "role": "r", "stratum": "s1",
         "source_state": "sclc", "goal_state": "luad", "alt_state": "normal",
         "perturbation_type": "overexpress", "donor_balanced_shift": -0.2},
        {"gene": "G2", "ensembl_id": "E2", "role": "r", "stratum": "s1",
         "source_state": "sclc", "goal_state": "luad", "alt_state": "normal",
         "perturbation_type": "delete", "donor_balanced_shift": None},
        {"gene": "G2", "ensembl_id": "E2", "role": "r", "stratum": "s1",
         "source_state": "sclc", "goal_state": "luad", "alt_state": "normal",
         "perturbation_type": "overexpress", "donor_balanced_shift": 0.1},
    ])
    E = astats.compute_E(rows)
    g1 = E[E.gene == "G1"].iloc[0]
    assert abs(g1["E"] - 0.3) < 1e-9, g1["E"]  # (|0.4| + |-0.2|)/2
    g2 = E[E.gene == "G2"].iloc[0]
    assert pd.isna(g2["E"]), "E must be NaN when one arm is missing, not silently the other arm's abs value"
    print("compute_E: correct value when both arms present, NaN when one is missing -- OK")


def test_midrank_percentile():
    # value tied with two controls, out of a 5-control reference set
    controls = [1.0, 2.0, 2.0, 3.0, 4.0]
    q = astats.midrank_percentile(2.0, controls)
    # combined = [1,2,2,3,4,2] (6 points); ranks (average, ascending) of the
    # three 2.0's (positions 2,3,6 in sorted order) = average of ranks 2,3,4 = 3
    assert abs(q - 100.0 * (3 - 0.5) / 6) < 1e-9, q
    assert astats.midrank_percentile(float("nan"), controls) is not None and np.isnan(astats.midrank_percentile(float("nan"), controls))
    assert np.isnan(astats.midrank_percentile(1.0, []))
    print("midrank_percentile: tie handling and empty/NaN inputs -- OK")


def test_spearman_exact_permutation_matches_naive_bruteforce():
    # Small n (5 -> 120 permutations) so a naive per-permutation
    # scipy.stats.spearmanr loop is fast enough to be the independent check.
    rng = np.random.default_rng(20260922)
    Q = rng.normal(size=5)
    risk = rng.normal(size=5)

    rho_fast, p_fast, perm_rhos_fast = astats.spearman_exact_permutation(Q, risk)

    naive_rhos = []
    observed_naive, _ = spearmanr(Q, risk)
    for perm in itertools.permutations(range(5)):
        r, _ = spearmanr(Q, risk[list(perm)])
        naive_rhos.append(r)
    naive_rhos = np.array(naive_rhos)
    p_naive = float(np.mean(np.abs(naive_rhos) >= np.abs(observed_naive) - 1e-9))

    assert abs(rho_fast - observed_naive) < 1e-9, (rho_fast, observed_naive)
    assert np.allclose(np.sort(perm_rhos_fast), np.sort(naive_rhos), atol=1e-9)
    assert abs(p_fast - p_naive) < 1e-9, (p_fast, p_naive)
    print(f"spearman_exact_permutation (n=5, 120 perms): matches naive brute force exactly "
          f"(rho={rho_fast:.4f}, p={p_fast:.4f}) -- OK")


def test_spearman_exact_permutation_n10_runs_and_is_symmetric():
    # The real design size (n=10 non-anchor genes, 3,628,800 permutations).
    # Correctness already proven at n=5 above; this checks it actually
    # completes in reasonable time/memory and gives a sane p-value, plus an
    # exact self-consistency check: permuting an IDENTICAL vector against
    # itself must give rho=1.0 and p equal to 1/10! (only the identity
    # permutation ties the observed extreme).
    identical = np.arange(10, dtype=float)
    started = time.time()
    rho, p, perm_rhos = astats.spearman_exact_permutation(identical, identical.copy())
    elapsed = time.time() - started
    assert abs(rho - 1.0) < 1e-9, rho
    assert len(perm_rhos) == 3_628_800
    # exactly the identity permutation reaches rho=+1.0 among 10! permutations
    n_at_max = int(np.sum(perm_rhos >= 1.0 - 1e-9))
    assert n_at_max == 1, n_at_max
    # two-sided: the full reversal permutation also hits rho=-1.0 (|rho|=1
    # too, for a strictly monotonic 10-point vector), so exactly 2 of the
    # 10! permutations are "as extreme", not 1 -- the two-sidedness the
    # design explicitly requires, not a bug in this check.
    assert abs(p - 2.0 / 3_628_800) < 1e-12, p
    print(f"spearman_exact_permutation (n=10, 3,628,800 perms): ran in {elapsed:.1f}s, "
          f"self-correlation rho=1.0 with two-sided p=2/10! exactly -- OK")


def test_compute_Q_and_stratum_gate():
    E_by_gene = {"S100X": 0.5}
    stratum_by_gene = {"S100X": "strat1"}
    enough_controls = {"strat1": {f"c{i}": float(i) for i in range(20)}}
    too_few_controls = {"strat1": {f"c{i}": float(i) for i in range(19)}}

    Q_ok = astats.compute_Q_for_contrast(E_by_gene, stratum_by_gene, enough_controls)
    assert not np.isnan(Q_ok["S100X"])

    Q_blocked = astats.compute_Q_for_contrast(E_by_gene, stratum_by_gene, too_few_controls)
    assert np.isnan(Q_blocked["S100X"]), "a stratum with <20 common controls must yield NaN Q, never a computed value from a widened set"
    print("compute_Q_for_contrast: 20 controls -> real Q, 19 controls -> not-estimable NaN -- OK")


def test_primary_test_gate():
    non_anchor = [f"g{i}" for i in range(10)]
    rng = np.random.default_rng(1)
    risk = {g: v for g, v in zip(non_anchor, rng.permutation(10))}
    # perfectly monotonic Q vs risk -> rho should be 1.0, gate passes
    Q_perfect = {g: risk[g] for g in non_anchor}
    result = astats.primary_test(Q_perfect, risk, non_anchor)
    assert abs(result["rho"] - 1.0) < 1e-9
    assert result["gate_pass"] is True

    Q_missing = dict(Q_perfect)
    Q_missing[non_anchor[0]] = np.nan
    result_missing = astats.primary_test(Q_missing, risk, non_anchor)
    assert result_missing["gate_pass"] is False and np.isnan(result_missing["rho"])
    print("primary_test: perfect monotonic Q passes the gate; a missing gene's Q fails it, not silently excluded -- OK")


def test_loo_and_bootstrap_gates():
    non_anchor = [f"g{i}" for i in range(10)]
    rng = np.random.default_rng(2)
    risk = {g: float(v) for g, v in zip(non_anchor, range(10))}
    E_by_gene = {g: float(v) for g, v in zip(non_anchor, range(10))}  # E perfectly tracks risk rank
    stratum_by_gene = {g: "strat1" for g in non_anchor}
    # 20 controls whose E values are unrelated to the S100 genes' own E/rank
    # (controls only set the PERCENTILE scale; a stable, evenly-spread
    # control distribution should leave rho basically unaffected by
    # dropping/resampling any one of them here).
    control_E = {"strat1": {f"c{i}": float(i) for i in range(20)}}

    loo_df = astats.leave_one_control_out(E_by_gene, stratum_by_gene, control_E, risk, non_anchor)
    assert len(loo_df) == 20, len(loo_df)
    print(f"leave_one_control_out: 20 rows for 1 stratum x 20 controls, "
          f"loo_gate_passes={astats.loo_gate_passes(loo_df)} -- OK")

    rhos = astats.bootstrap_stability(E_by_gene, stratum_by_gene, control_E, risk, non_anchor,
                                       seed=20260922, n_boot=500)
    assert len(rhos) == 500
    rhos_repeat = astats.bootstrap_stability(E_by_gene, stratum_by_gene, control_E, risk, non_anchor,
                                              seed=20260922, n_boot=500)
    assert np.array_equal(rhos, rhos_repeat, equal_nan=True), "same seed must reproduce identical bootstrap draws"
    print(f"bootstrap_stability: 500 draws, seed-reproducible, "
          f"bootstrap_gate_passes={astats.bootstrap_gate_passes(rhos)} -- OK")


def test_numerical_floor():
    assert astats.above_numerical_floor(0.001) is True
    assert astats.above_numerical_floor(0.0001) is False
    assert astats.above_numerical_floor(astats.NUMERICAL_FLOOR) is False  # strictly greater-than, not >=
    assert astats.above_numerical_floor(None) is None
    assert astats.above_numerical_floor(float("nan")) is None
    print(f"above_numerical_floor (floor={astats.NUMERICAL_FLOOR}): strict >, None for missing -- OK")


def test_matched_controls_stub_names_both_gaps():
    try:
        mc.build_matched_control_table(mc.STRATA, "unused.csv", "unused.csv", per_source_state=True)
        raise AssertionError("build_matched_control_table must raise -- it is blocked, not implemented")
    except NotImplementedError as exc:
        msg = str(exc)
        assert "median token rank" in msg, "stub must name open question 1"
        assert "per source-state or globally" in msg or "per-source-state or globally" in msg, "stub must name open question 2"
        print("build_matched_control_table: raises NotImplementedError naming both open questions -- OK")


def test_synthetic_control_table_is_structurally_valid_and_labeled_fake():
    table = mc.synthetic_control_table()
    assert set(table.keys()) == set(mc.STRATA.keys())
    for stratum, controls in table.items():
        assert len(controls) == mc.MIN_COMMON_CONTROLS, (stratum, len(controls))
        assert all(c.startswith("ENSG9SYNTH") for c in controls), "fake ids must be unmistakably fake"
    print("synthetic_control_table: 20 fabricated, unmistakably-fake ids per stratum, all 6 strata -- OK")


def test_end_to_end_with_synthetic_controls():
    """Everything above the stub, run together once against
    synthetic_control_table(), per Michael's instruction: 'the interface
    is only real once something has run through it.'"""
    panel = [
        {"gene": g, "ensembl_id": f"ENSG_{g}", "stratum": "clean_high", "role": "clean_control"}
        for g in ("S100A4", "S100A6", "S100A10", "S100A11")
    ] + [
        {"gene": g, "ensembl_id": f"ENSG_{g}", "stratum": "low_p_a16", "role": "flagged"}
        for g in ("S100P", "S100A16")
    ] + [
        {"gene": "S100A2", "ensembl_id": "ENSG_S100A2", "stratum": "low_a2_b", "role": "flagged_non_anchor"},
        {"gene": "S100B", "ensembl_id": "ENSG_S100B", "stratum": "low_a2_b", "role": "intermediate"},
        {"gene": "S100A13", "ensembl_id": "ENSG_S100A13", "stratum": "singleton_a13", "role": "intermediate"},
        {"gene": "S100A8", "ensembl_id": "ENSG_S100A8", "stratum": "singleton_a8", "role": "ambient_anchor"},
        {"gene": "S100A9", "ensembl_id": "ENSG_S100A9", "stratum": "low_a9_pbp", "role": "ambient_anchor"},
        {"gene": "S100PBP", "ensembl_id": "ENSG_S100PBP", "stratum": "low_a9_pbp", "role": "intermediate"},
    ]
    non_anchor_genes = [g["gene"] for g in panel if g["role"] != "ambient_anchor"]
    assert len(non_anchor_genes) == 10
    stratum_by_gene = {g["gene"]: g["stratum"] for g in panel}

    rng = np.random.default_rng(20260922)
    E_by_gene = {g["gene"]: float(rng.uniform(0, 1)) for g in panel}
    risk_by_gene = {g: float(rng.uniform(0, 1)) for g in non_anchor_genes}

    control_table = mc.synthetic_control_table()
    control_E_by_stratum = {
        stratum: {cid: float(rng.uniform(0, 1)) for cid in cids}
        for stratum, cids in control_table.items()
    }

    Q_by_gene = astats.compute_Q_for_contrast(E_by_gene, stratum_by_gene, control_E_by_stratum)
    assert all(not np.isnan(Q_by_gene[g["gene"]]) for g in panel), "every gene has 20 controls in this fixture, so every Q must be estimable"
    summary = astats.stratum_summary(Q_by_gene, stratum_by_gene)
    assert len(summary) == 6

    result = astats.primary_test(Q_by_gene, risk_by_gene, non_anchor_genes)
    assert not np.isnan(result["rho"])

    loo_df = astats.leave_one_control_out(E_by_gene, stratum_by_gene, control_E_by_stratum, risk_by_gene, non_anchor_genes)
    assert len(loo_df) == 6 * 20

    boot_rhos = astats.bootstrap_stability(E_by_gene, stratum_by_gene, control_E_by_stratum, risk_by_gene,
                                            non_anchor_genes, n_boot=200)
    assert len(boot_rhos) == 200

    print(f"end-to-end (synthetic controls): Q computed for all 12 genes, 6-stratum rollup, "
          f"primary rho={result['rho']:.3f} p={result['p_exact']:.4g}, "
          f"120 LOO rows, 200 bootstrap draws -- all ran through the same interface build_matched_control_table() will satisfy -- OK")


def main() -> None:
    test_compute_E()
    test_midrank_percentile()
    test_spearman_exact_permutation_matches_naive_bruteforce()
    test_spearman_exact_permutation_n10_runs_and_is_symmetric()
    test_compute_Q_and_stratum_gate()
    test_primary_test_gate()
    test_loo_and_bootstrap_gates()
    test_numerical_floor()
    test_matched_controls_stub_names_both_gaps()
    test_synthetic_control_table_is_structurally_valid_and_labeled_fake()
    test_end_to_end_with_synthetic_controls()
    print("\nALL AMBIENT-STATS CHECKS PASSED")


if __name__ == "__main__":
    main()
