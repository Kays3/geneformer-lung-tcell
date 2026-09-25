"""build_design.py against the COMMITTED pre-GPU inputs (no Phase 6 output exists in the repo), plus
mutated copies that must be refused."""
import csv
import json
import os
import shutil
import sys

import pytest

HERE = os.path.dirname(__file__)
BASE = os.environ.get("BDL_BASE", os.path.join(HERE, ".."))
sys.path[:0] = [os.path.join(HERE, "..", "scripts")]
import analyse as an  # noqa: E402
import build_design as bd  # noqa: E402


@pytest.fixture(scope="module")
def real():
    return bd.build(BASE)[0]


@pytest.fixture
def copy(tmp_path):
    for rel in bd.INPUTS.values():
        dst = tmp_path / rel
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy(os.path.join(BASE, rel), dst)
    return tmp_path


def edit_csv(path, fn):
    rows = list(csv.DictReader(open(path, newline="")))
    fields = list(rows[0])
    rows = fn(rows)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)


def edit_json(path, fn):
    d = json.load(open(path))
    fn(d)
    json.dump(d, open(path, "w"))


# ---------------------------------------------------------------- the real design
def test_real_design_counts(real):
    assert len(real["panel_a"]) == 15 and len(real["panel_b"]) == 36 and len(real["panel_b_not_run"]) == 3
    assert {g["symbol"] for g in real["panel_b_not_run"]} == {"TRAC", "TRBC1", "TRBC2"}
    assert len(real["donors"]) == 43 and set(real["donor_fold"].values()) == {0, 1, 2, 3, 4}
    assert real["not_run_genes"] == ["ENSG00000133639"]                       # BTG1
    assert real["strata"]["A01"]["status"] == "not_estimable_control_stratum"


def test_real_design_loads_into_analyse(real):
    d = an.Design(**real)
    assert len(d.expected_estimable) == 51 and set(d.ambient_flagged) == {g["ensembl_id"] for g in d.panel_a}


def test_real_design_is_deterministic(tmp_path):
    for n in ("a", "b"):
        assert bd.main(["--base", BASE, "--out", str(tmp_path / f"{n}.json")]) == 0
    assert (tmp_path / "a.json").read_bytes() == (tmp_path / "b.json").read_bytes()
    prov = json.load(open(tmp_path / "a_inputs.json"))
    assert set(prov) == set(bd.INPUTS) and all(len(v["sha256"]) == 64 for v in prov.values())


# ---------------------------------------------------------------- must be REFUSED
def test_refuses_eligible_flag_disagreeing_with_count(copy):
    p = copy / bd.INPUTS["eligibility"]
    edit_csv(p, lambda rs: [dict(r, estimable_donors="9") if r["symbol"] == "POLR2J3" else r for r in rs])
    with pytest.raises(bd.DesignError, match="eligible flag"):
        bd.build(str(copy))


def test_refuses_donor_missing_from_classifier_gate(copy):
    p = copy / bd.INPUTS["classifier_gate"]
    edit_json(p, lambda d: d["per_donor_balanced_accuracy"].pop(sorted(d["per_donor_balanced_accuracy"])[0]))
    with pytest.raises(bd.DesignError, match="classifier gate donors"):
        bd.build(str(copy))


def test_refuses_ambient_table_missing_a_panel_a_gene(copy):
    p = copy / bd.INPUTS["ambient"]
    edit_csv(p, lambda rs: rs[1:])
    with pytest.raises(bd.DesignError, match="ambient table"):
        bd.build(str(copy))


def test_refuses_stratum_membership_mismatch(copy):
    p = copy / bd.INPUTS["strata"]
    edit_json(p, lambda d: d["B00"]["members"].append("ENSG00000188389"))
    with pytest.raises(bd.DesignError, match="strata members"):
        bd.build(str(copy))


def test_refuses_panel_gene_used_as_control(copy):
    p = copy / bd.INPUTS["strata"]
    edit_json(p, lambda d: d["B00"]["controls"].__setitem__(0, "ENSG00000168255"))
    with pytest.raises(bd.DesignError, match="panel gene is used as a control"):
        bd.build(str(copy))


def test_refuses_host_map_that_differs_from_run_set(copy):
    p = copy / bd.INPUTS["host_map"]
    edit_csv(p, lambda rs: rs[:-1])
    with pytest.raises(bd.DesignError, match="host map"):
        bd.build(str(copy))


def test_refuses_missing_input(copy):
    os.remove(copy / bd.INPUTS["folds"])
    with pytest.raises(bd.DesignError, match="missing input folds"):
        bd.build(str(copy))


def test_main_exits_nonzero_on_refusal(copy, tmp_path):
    os.remove(copy / bd.INPUTS["folds"])
    assert bd.main(["--base", str(copy), "--out", str(tmp_path / "x.json")]) == 1
    assert not (tmp_path / "x.json").exists()
