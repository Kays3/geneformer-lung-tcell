"""Unit tests for ISP-BAL-1 (synthetic data only).

Run: python -m pytest protocol_criteria_e3/tests -q

Each test builds the situation a criterion exists for and checks the verdict, so
a failing test means a criterion no longer catches what it was written to catch.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "criteria"))

import isp_criteria as c  # noqa: E402


# ------------------------------------------------------------------ A: design

def test_a1_donor_disjoint_passes_and_catches_a_leak():
    ok = pd.DataFrame({"donor": ["d1", "d1", "d2"], "partition": ["train", "train", "test"]})
    assert c.check_donor_disjoint(ok)["pass"]
    leak = pd.DataFrame({"donor": ["d1", "d1", "d2"], "partition": ["train", "test", "test"]})
    res = c.check_donor_disjoint(leak)
    assert not res["pass"] and res["leaking_donors"] == ["d1"]


def test_a2_study_confounded_classes_fail():
    # The July defect: no study contributes both a normal and a tumour T cell.
    confounded = pd.DataFrame({"study": ["s1", "s1", "s2", "s3"], "label": ["normal", "normal", "LUAD", "LUAD"]})
    res = c.check_class_confounding(confounded)
    assert not res["pass"] and res["unbridged_pairs"] == [("LUAD", "normal")]
    bridged = pd.DataFrame({"study": ["s1", "s1", "s2"], "label": ["normal", "LUAD", "LUAD"]})
    assert c.check_class_confounding(bridged)["pass"]


def test_a3_confounder_balance_categorical_and_numeric():
    rng = np.random.default_rng(1)
    balanced = pd.DataFrame({"partition": ["a"] * 40 + ["b"] * 40, "sex": ["F", "M"] * 40,
                             "age": rng.normal(60, 5, 80)})
    assert c.check_confounder_balance(balanced, "sex")["pass"]
    assert c.check_confounder_balance(balanced, "age")["pass"]
    skewed = pd.DataFrame({"partition": ["a"] * 40 + ["b"] * 40, "sex": ["F"] * 40 + ["M"] * 40,
                           "age": np.r_[rng.normal(40, 3, 40), rng.normal(70, 3, 40)]})
    assert not c.check_confounder_balance(skewed, "sex")["pass"]
    assert not c.check_confounder_balance(skewed, "age")["pass"]


def test_a4_equal_donor_weight():
    assert c.check_equal_donor_weight(pd.Series({"d1": 100, "d2": 100}), 100)["pass"]
    assert not c.check_equal_donor_weight(pd.Series({"d1": 100, "d2": 2220}), 100)["pass"]


# -------------------------------------------------------------- B: classifier

def _predictions(n_donors, accuracy, seed=0, classes=("normal", "tumour")):
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(n_donors):
        for cls in classes:  # donor-paired: each donor has every class
            for _ in range(50):
                pred = cls if rng.random() < accuracy else rng.choice([x for x in classes if x != cls])
                rows.append({"donor": f"d{d}", "y_true": cls, "y_pred": pred})
    return pd.DataFrame(rows)


def test_b_gate_passes_a_good_classifier():
    res = c.classifier_gate(_predictions(19, 0.9), n_boot=200)
    assert res["B1_pass"] and res["B2_pass"] and res["pass"]
    assert res["donors_above_chance"] == 19 and res["sign_test_p"] < 1e-4


def test_b_gate_fails_guessing_and_permutation_check_uses_it():
    res = c.classifier_gate(_predictions(19, 0.5), n_boot=200)
    assert not res["pass"]
    assert c.permutation_calibration(res)["pass"]  # B4: a permuted-label run must fail


def test_b3_noise_band_is_reported_without_changing_the_decision():
    res = c.classifier_gate(_predictions(19, 0.62, seed=3), n_boot=200)
    assert res["B3_near_threshold"]


# ------------------------------------------------- C: perturbation integrity

def test_c1_noop_gate():
    assert c.noop_gate([0.0, 0.0, 0.0])["pass"]
    assert not c.noop_gate([0.0, 1e-9])["pass"]


def test_c2_precision_agreement():
    rng = np.random.default_rng(2)
    a = rng.normal(0, 1, 300)
    assert c.precision_agreement(a, a + rng.normal(0, 1e-3, 300))["pass"]
    # Two model sizes ranking genes differently (as 104M against 316M did) fail.
    assert not c.precision_agreement(a, 0.5 * a + rng.normal(0, 1, 300))["pass"]


def test_c3_cross_host_equivalence():
    a = np.linspace(-0.01, 0.01, 100)
    assert c.cross_host_equivalence(a, a + 1e-5)["pass"]
    assert not c.cross_host_equivalence(a, a + 5e-3)["pass"]


def test_c4_config_consistency():
    cfg = {"emb_layer": 0, "emb_mode": "cls", "filter_data": {"cell_type": ["T"]}}
    assert c.config_consistency(cfg, dict(cfg))["pass"]
    assert not c.config_consistency(cfg, {**cfg, "emb_layer": -1})["pass"]


# -------------------------------------------------- multiple testing helpers

def test_holm_and_bh_match_known_values():
    p = [0.01, 0.04, 0.03, 0.005]
    assert np.allclose(c.holm(p), [0.03, 0.06, 0.06, 0.02])
    assert np.allclose(c.bh(p), [0.02, 0.04, 0.04, 0.02])


def test_attainable_p_and_holm_family_limit():
    # 9 DCM donors in the protocol's atlas: min p = 2/512, so Holm caps families at 12.
    assert c.attainable_min_p(9) == pytest.approx(2 / 512)
    assert c.max_holm_family(9) == 12
    assert c.max_holm_family(19) > 10_000


def test_family_correction_switches_to_bh_when_holm_cannot_pass():
    _, method = c.family_correction([0.004] * 12, n_donors=9)
    assert method == "HOLM"
    adj, method = c.family_correction([0.004] * 13, n_donors=9)
    assert method == "BH_POWER_ADAPTED" and np.all(adj <= 0.05)


# --------------------------------------------------- D: gene-level evidence

def _null_cells(n_genes, n_donors, cells_per_donor, donor_sd, seed):
    """Per-cell shifts of random genes. donor_sd > 0 gives each (gene, donor) its own
    offset, as real donors do; cells inside a donor are then not independent."""
    rng = np.random.default_rng(seed)
    out = {}
    for g in range(n_genes):
        offsets = rng.normal(0, donor_sd, n_donors)
        out[f"null{g}"] = np.concatenate([rng.normal(o, 0.01, cells_per_donor) for o in offsets])
    return out


def test_d2_screen_is_calibrated_when_cells_are_independent():
    res = c.screen_calibration(_null_cells(60, 10, 100, donor_sd=0.0, seed=4))
    assert res["pass"], res
    # A rate a little above 0.05 by chance is not evidence of miscalibration.
    assert res["fpr_ci95"][0] <= 0.05


def test_d2_screen_is_flagged_under_pseudo_replication():
    # The failure P4 is exposed to: thousands of cells from a few donors, scored as
    # independent. Random genes then look significant far beyond 5%.
    res = c.screen_calibration(_null_cells(60, 10, 100, donor_sd=0.01, seed=5))
    assert not res["pass"] and res["fpr"] > 0.5


def test_d3_d4_d5_donor_level_test():
    donors = [f"d{i}" for i in range(10)]
    gene = pd.Series(0.02, index=donors)
    controls = pd.DataFrame(np.random.default_rng(6).normal(0, 0.001, (10, 20)), index=donors)
    adj, status = c.control_adjust(gene, controls)
    assert status == "OK"
    res = c.donor_test(adj)
    assert res["p"] == pytest.approx(2 / 2 ** 10) and res["median"] > 0
    _, status = c.control_adjust(gene, controls.iloc[:, :19])
    assert status == "NOT_ESTIMABLE_CONTROLS"


def test_d6_dose_concordance_and_alt_state():
    assert c.dose_concordance(True, 0.1, True, -0.1) == "COHERENT"
    assert c.dose_concordance(True, 0.1, True, 0.1) == "INCOHERENT"  # the S100A8/A9 pattern
    assert c.dose_concordance(True, 0.1, False, 0.0) == "UNRESOLVED"
    assert c.alt_state_ok(0.02, -0.01) is True
    assert c.alt_state_ok(0.02, 0.03) is False
    assert c.alt_state_ok(0.02, None) is None


def test_d7_beats_null_and_null_pair_correlation():
    rng = np.random.default_rng(7)
    null = rng.normal(0, 0.01, 100)
    assert c.beats_null(0.05, null)["pass"]
    assert not c.beats_null(0.005, null)["pass"]
    dele = rng.normal(0, 1, 100)
    res = c.null_pair_correlation(dele, -0.6 * dele + rng.normal(0, 0.5, 100), n_perm=500)
    assert res["rho"] < -0.5 and res["p_lower"] < 0.01


def test_d8_fold_stability():
    assert c.fold_stable([0.1, 0.2, -0.05, 0.1, 0.1], 0.1)["pass"]
    assert not c.fold_stable([0.1, -0.2, -0.05, 0.1, 0.1], 0.1)["pass"]


# ------------------------------------------------------- E: evidence levels

def test_e1_levels():
    base = dict(estimable=True, donor_sig=True, concordance="COHERENT", beats_null=True,
                fold_stable=True, ambient_clear=True, alt_ok=None)
    assert c.evidence_level(base) == "SUPPORTED"
    assert c.evidence_level({**base, "replicated": True}) == "REPLICATED"
    assert c.evidence_level({**base, "validated": True}) == "VALIDATED"
    assert c.evidence_level({**base, "concordance": "INCOHERENT"}) == "OPEN"
    assert c.evidence_level({**base, "ambient_clear": False}) == "OPEN"
    assert c.evidence_level({**base, "alt_ok": False}) == "OPEN"
    screen_only = dict(estimable=True, screen_sig=True, screen_calibrated=True)
    assert c.evidence_level(screen_only) == "SCREEN"
    assert c.evidence_level({**screen_only, "screen_calibrated": False}) == "OPEN"
    assert c.evidence_level({"estimable": False}) == "NOT_ESTIMABLE"


def _synthetic_run(seed=8, n_donors=12, cells=30):
    """One run: a real gene (consistent across donors, opposite under overexpression),
    a pseudo-replicated gene (a moderate shift in 3 of 12 donors, the same sign under
    both operations), 20 matched controls and 40 null genes, five folds."""
    rng = np.random.default_rng(seed)
    rows = []
    donors = [f"d{i}" for i in range(n_donors)]

    def add(gene, role, op, means):
        for d, m in zip(donors, means):
            for s in rng.normal(m, 0.002, cells):
                rows.append(dict(gene=gene, role=role, stratum="s1", op=op, donor=d, shift=s,
                                 fold=int(d[1:]) % 5))

    for op, sign in (("delete", 1), ("overexpress", -1)):
        add("REAL", "test", op, [sign * 0.02] * n_donors)
        add("FEW_DONORS", "test", op, [0.05] * 3 + [0.0] * (n_donors - 3))
        for k in range(20):
            add(f"ctl{k}", "control", op, rng.normal(0, 0.002, n_donors))
        for k in range(40):
            add(f"null{k}", "null", op, rng.normal(0, 0.003, n_donors))
    return pd.DataFrame(rows)


def test_evaluate_genes_separates_a_real_effect_from_pseudo_replication():
    res = c.evaluate_genes(_synthetic_run(), ["REAL", "FEW_DONORS"]).set_index("gene")
    assert res.loc["REAL", "level"] == "SUPPORTED"
    assert res.loc["REAL", "concordance"] == "COHERENT"
    # Shifted in 3 of 12 donors: the cell-pooled screen calls it, the donor test does not.
    assert bool(res.loc["FEW_DONORS", "screen_sig"])
    assert res.loc["FEW_DONORS", "level"] != "SUPPORTED"
    assert not bool(res.loc["REAL", "screen_calibrated"])  # donor offsets make the screen overcall here
