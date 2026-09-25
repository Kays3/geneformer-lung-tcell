"""Phase 7 analysis against SYNTHETIC Phase 6 output with planted answers (tests/synth.py).
Each registered status must be recovered, and each broken input must be refused, not absorbed."""
import json
import os
import sys

import numpy as np
import pytest

HERE = os.path.dirname(__file__)
sys.path[:0] = [os.path.join(HERE, "..", "scripts"), HERE]
import analyse as an  # noqa: E402
import stats_core as sc  # noqa: E402
import synth  # noqa: E402

RULES = {"min_controls_per_donor": 10, "holm_m": "panel"}


def genes_all():
    return sorted([g for g, v in synth.PLANT.items() if v[1] is not None and v[3] == "S0"] + synth.CONTROLS)


@pytest.fixture(scope="module")
def tree(tmp_path_factory):
    root = str(tmp_path_factory.mktemp("phase6"))
    synth.build(root)
    return root


@pytest.fixture(scope="module")
def result(tree):
    d = an.Design(**synth.design())
    calls = an.load_all(tree, genes_all(), d.donors)
    primary = an.analyse(calls, d, RULES)
    return d, calls, primary


def test_every_planted_status_is_recovered(result):
    _, _, primary = result
    got = {r["gene"]: r["status"] for r in primary["rows"]}
    want = {g: v[4] for g, v in synth.PLANT.items()}
    want["GB_TRAC"] = "NOT_RUN"
    assert got == want


def test_s3_downgrade_is_labelled(result):
    row = next(r for r in result[2]["rows"] if r["gene"] == "GB_FOLD")
    assert row["s3"].startswith("fold-heterogeneous (downgraded from T_CELL_SIGNAL_TOWARD)")
    assert sum(v < 0 for v in row["fold_medians_delete"].values()) == 2


def test_panel_a_non_replications_carry_the_registered_qualifier(result):
    for r in result[2]["rows"]:
        if r["panel"] == "A" and r["status"] in sc.NON_REPLICATION:
            assert r["claim_qualifier"] == sc.PANEL_A_QUALIFIER
        assert r["design_warrant"] == sc.WARRANT_PAIRED


def test_not_run_rows_state_their_reason(result):
    rows = {r["gene"]: r for r in result[2]["rows"]}
    assert rows["GA_LOWD"]["reason"] == "estimable donors = 4 < d_min = 10"
    assert rows["GA_NOCTL"]["status"] == "NOT_ESTIMABLE_CONTROLS"
    assert "token dictionary" in rows["GB_TRAC"]["reason"]


def test_control_adjustment_removes_a_shared_offset(tmp_path):
    # Add +0.5 to EVERY gene (panel and controls) in every donor: a must not change.
    root = str(tmp_path)
    synth.build(root)
    d = an.Design(**synth.design())
    calls = an.load_all(root, genes_all(), d.donors)
    base = an.analyse(calls, d, RULES)["arms"][("GA_REP", "delete")]
    for key, per_d in calls.items():
        for c in per_d.values():
            if c["status"] == "done":
                c["shifts"] = {s: v + 0.5 for s, v in c["shifts"].items()}
    shifted = an.analyse(calls, d, RULES)["arms"][("GA_REP", "delete")]
    assert all(abs(base[k] - shifted[k]) < 1e-9 for k in base)


def test_donor_values_come_from_the_goal_state_not_another(result):
    d, calls, primary = result
    donor = d.donors[0]
    s = np.mean(calls[("delete", "GA_REP")][donor]["shifts"]["normal_adjacent"])
    ctrl = np.median([np.mean(calls[("delete", c)][donor]["shifts"]["normal_adjacent"]) for c in synth.CONTROLS])
    assert abs(primary["arms"][("GA_REP", "delete")][donor] - (s - ctrl)) < 1e-12


# ---------------------------------------------------------------- must FAIL on broken input
def test_missing_marker_refuses_partial_arm(tmp_path):
    root = str(tmp_path); synth.build(root)
    os.remove(os.path.join(root, "delete", "GA_REP", f"{synth.DONORS[3]}.complete.json"))
    with pytest.raises(an.IncompleteRun, match="partial arm"):
        an.load_all(root, genes_all(), synth.DONORS)


def test_truncated_marker_is_refused(tmp_path):
    root = str(tmp_path); synth.build(root)
    p = os.path.join(root, "overexpress", "CTRL05", f"{synth.DONORS[0]}.complete.json")
    open(p, "w").write('{"gene": "CTRL05", "op"')        # a kill mid-write
    with pytest.raises(an.IncompleteRun, match="unparseable marker"):
        an.load_all(root, genes_all(), synth.DONORS)


def test_cell_count_mismatch_is_refused(tmp_path):
    root = str(tmp_path); synth.build(root)
    p = os.path.join(root, "delete", "GB_TOW", f"{synth.DONORS[1]}.complete.json")
    rec = json.load(open(p)); rec["n_token_cells"] = 29; json.dump(rec, open(p, "w"))
    with pytest.raises(an.IncompleteRun, match="per-cell count"):
        an.load_all(root, genes_all(), synth.DONORS)


def test_estimable_count_disagreeing_with_pre_gpu_table_is_refused(tree):
    dd = synth.design(); dd["expected_estimable"]["GA_REP"] = 39
    d = an.Design(**dd)
    calls = an.load_all(tree, genes_all(), d.donors)
    with pytest.raises(ValueError, match="pre-GPU table"):
        an.analyse(calls, d, RULES)


def test_validator_refuses_a_row_without_the_qualifier():
    row = sc.make_row("A", "X", "OPEN"); row.pop("claim_qualifier")
    with pytest.raises(ValueError):
        sc.validate_rows([row])


# ---------------------------------------------------------------- the two open rules behave as specified
def test_min_controls_rule_drops_donors_below_it(tmp_path):
    root = str(tmp_path)
    synth.build(root, ctrl_missing={synth.DONORS[0]: 11, synth.DONORS[1]: 10})   # 9 and 10 controls left
    d = an.Design(**synth.design())
    calls = an.load_all(root, genes_all(), d.donors)
    a = an.analyse(calls, d, RULES)["arms"][("GA_REP", "delete")]
    assert synth.DONORS[0] not in a and synth.DONORS[1] in a


def test_holm_panel_padding_is_stricter_than_tested_only():
    ps = {"g1": sc.Fraction(1, 100)}
    fam = ["g1"] + [f"u{i}" for i in range(9)]
    assert an.holm_family(ps, fam, {"holm_m": "tested"})["g1"] == sc.Fraction(1, 100)
    assert an.holm_family(ps, fam, {"holm_m": "panel"})["g1"] == sc.Fraction(10, 100)


def test_secondary_floor_blocks_small_panels(result):
    sec = result[2]["secondaries"]
    assert sec["A_all"]["status"].startswith("not run: n_eligible = 4 < 12")
    assert sec["B_all"]["status"].startswith("not run: n_eligible = 5 < 12")


def test_sensitivities_never_change_primary_status(result):
    d, calls, primary = result
    before = {r["gene"]: r["status"] for r in primary["rows"]}
    sens = an.sensitivities(calls, d, primary, RULES)
    after = {r["gene"]: r["status"] for r in primary["rows"]}
    assert before == after
    assert "GA_REP|delete" in sens["S2_leader_merad_vs_rest"] and "panel_B|delete" in sens["A3f_accuracy_vs_effect"]
    assert sens["A3f_accuracy_vs_effect"]["GA_REP|delete"]["conditioned_on_significance"] is True
