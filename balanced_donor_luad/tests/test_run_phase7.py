"""The Phase 7 launch gate: all gates satisfied passes; each kind of breach blocks (Amendment 3g.5)."""
import json
import os
import shutil
import sys

import pytest

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
import run_phase7 as rp  # noqa: E402

MANIFEST = os.path.join(HERE, "..", "registration", "required_gates.json")


def _write(base, rel, obj):
    f = os.path.join(base, rel); os.makedirs(os.path.dirname(f), exist_ok=True)
    if isinstance(obj, str):
        open(f, "w").write(obj)
    else:
        json.dump(obj, open(f, "w"))


def satisfied_base(tmp_path):
    b = str(tmp_path)
    _write(b, "phase4_results/classifier_gate.json", {"PASS": True})
    _write(b, "ambient/ambient_summary.json", {"valid": True})
    _write(b, "phase6_prep/pre_isp_gates_late.json", {"noop_gate": {"PASS": True}, "PASS": True})
    _write(b, "phase6_prep/probe_3e3_probe_3e3_result.json", {"delete": {"registered_PASS": True}, "overexpress": {"registered_PASS": True}})
    _write(b, "phase7_prep/probe_end_result.json", {"delete": {"registered_PASS": False}, "overexpress": {"registered_PASS": True}})
    _write(b, "phase6_prep/token_edit_verification.json", {"PASS": True})
    _write(b, "registration/PHASE5_ISP_REGISTRATION.md", "...\n## Amendment 3g — rules\n")
    _write(b, "registration/phase7_rules.json", {"control_min_cells": 10, "min_controls_per_donor": 10, "holm_m": "panel"})
    return b


def test_all_satisfied_passes_even_with_end_gate_fail(tmp_path):
    # an end-gate FAIL must not block: it sets host_drift instead (s.3.5)
    assert rp.check_gates(satisfied_base(tmp_path), MANIFEST) == []


@pytest.mark.parametrize("rel", ["phase4_results/classifier_gate.json", "phase6_prep/pre_isp_gates_late.json",
                                 "phase7_prep/probe_end_result.json", "registration/phase7_rules.json"])
def test_missing_file_blocks(tmp_path, rel):
    b = satisfied_base(tmp_path); os.remove(os.path.join(b, rel))
    assert any("missing" in f for f in rp.check_gates(b, MANIFEST))


def test_fail_value_blocks(tmp_path):
    b = satisfied_base(tmp_path)
    _write(b, "phase6_prep/pre_isp_gates_late.json", {"noop_gate": {"PASS": False}, "PASS": True})
    assert any(f.startswith("noop_gate:") for f in rp.check_gates(b, MANIFEST))


def test_changed_rules_block(tmp_path):
    b = satisfied_base(tmp_path)
    _write(b, "registration/phase7_rules.json", {"control_min_cells": 1, "min_controls_per_donor": 10, "holm_m": "panel"})
    assert any(f.startswith("rules_frozen:") for f in rp.check_gates(b, MANIFEST))


def test_absent_field_and_unparseable_file_block(tmp_path):
    b = satisfied_base(tmp_path)
    _write(b, "ambient/ambient_summary.json", {"auc": 0.9})
    _write(b, "phase4_results/classifier_gate.json", "{not json")
    fails = rp.check_gates(b, MANIFEST)
    assert any(f.startswith("ambient_produced:") for f in fails) and any(f.startswith("classifier_gate:") for f in fails)


def test_missing_amendment_text_blocks(tmp_path):
    b = satisfied_base(tmp_path)
    _write(b, "registration/PHASE5_ISP_REGISTRATION.md", "no amendment here\n")
    assert any(f.startswith("amendment_3g_registered:") for f in rp.check_gates(b, MANIFEST))
