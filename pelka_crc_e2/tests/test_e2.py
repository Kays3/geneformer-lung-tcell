"""Unit tests for the E2-specific scripts (synthetic data only). Run: python -m pytest pelka_crc_e2/tests -q"""
import os
import sys
from fractions import Fraction

import numpy as np
import pandas as pd
import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "scripts"))
sys.path.insert(0, os.path.join(HERE, "..", "..", "balanced_donor_luad", "scripts"))

import build_cohort_pelka as bc  # noqa: E402
import compare_panel_b as cp  # noqa: E402
import make_run_order as ro  # noqa: E402
import select_null_e2 as sn  # noqa: E402


# ---------------------------------------------------------------- cohort
def test_strip_version_handles_suffix_and_par_y():
    assert bc.strip_version(["ENSG00000243485.5_4", "ENSG00000182378.13_5_PAR_Y", "ENSG00000000003.14"]) == \
        ["ENSG00000243485", "ENSG00000182378", "ENSG00000000003"]


def test_strip_version_refuses_non_gene_ids():
    with pytest.raises(SystemExit):
        bc.strip_version(["CD4-SK3"])


def test_collapse_sums_duplicates():
    uniq, inv, n = bc.collapse_genes(["ENSG1", "ENSG2", "ENSG1"])
    assert list(uniq) == ["ENSG1", "ENSG2"] and list(inv) == [0, 1, 0] and n == 1


def _cells(counts):
    rows = []
    for (pid, spec, proc, ct), n in counts.items():
        rows += [{"cell_id": f"{pid}_{spec}_{proc}_{ct}_{i}", "PID": pid, "SPECIMEN_TYPE": spec,
                  "PROCESSING_TYPE": proc, "clMidwayPr": ct} for i in range(n)]
    return pd.DataFrame(rows)


def test_eligible_cells_refuses_when_donor_set_differs(monkeypatch):
    monkeypatch.setattr(bc, "REGISTERED_DONORS", ("P1",))
    x = _cells({("P1", "T", "unsorted", "TCD4"): 100, ("P1", "N", "unsorted", "TCD8"): 100,
                ("P2", "T", "unsorted", "TCD4"): 100, ("P2", "N", "unsorted", "TCD4"): 100})
    with pytest.raises(SystemExit):
        bc.eligible_cells(x)


def test_eligible_cells_filters_processing_and_type(monkeypatch):
    monkeypatch.setattr(bc, "REGISTERED_DONORS", ("P1",))
    x = _cells({("P1", "T", "unsorted", "TCD4"): 100, ("P1", "N", "unsorted", "TCD8"): 100,
                ("P1", "T", "CD45pMACS", "TCD4"): 500, ("P1", "N", "unsorted", "Myeloid"): 50,
                ("P2", "T", "CD45pMACS", "TCD4"): 300, ("P2", "N", "unsorted", "TCD4"): 300})
    e = bc.eligible_cells(x)
    assert len(e) == 200 and set(e.origin) == {"tumor_primary", "normal_adjacent"}
    assert set(e.PROCESSING_TYPE) == {"unsorted"}


def test_draw_matches_balanced_donor_draw():
    import build_cohort as lung
    e = pd.DataFrame({"cell_id": [f"c{i:04d}" for i in range(350)], "donor_id": "D", "origin": "tumor_primary"})
    a = bc.draw(e.copy()).sort_values("cell_id").draw_rank.to_numpy()
    b = lung.draw(e.copy()).sort_values("cell_id").draw_rank.to_numpy()
    assert np.array_equal(a, b)
    d = bc.draw(e.copy())
    assert d.in_analysis100.sum() == 100 and d.in_pool300.sum() == 300


# ---------------------------------------------------------------- null draw
def test_walk_stops_at_target_and_counts_estimable():
    tok = {"g1": 1, "g2": 2, "g3": 3}
    sets = {f"d{k}": [{1, 3}] * 10 + [set()] * 90 for k in range(10)}   # g1, g3 in 10 cells of all 10 donors
    rows, chosen = sn.walk(["g2", "g1", "g3"], tok, sets, target=1)
    assert chosen == ["g1"] and [r["gene"] for r in rows] == ["g2", "g1"]
    assert rows[0]["estimable"] is False and rows[1]["n_estimable_donors"] == 10


def test_walk_needs_ten_donors():
    tok = {"g1": 1}
    sets = {f"d{k}": ([{1}] * 10 if k < 9 else []) + [set()] * 90 for k in range(12)}
    with pytest.raises(SystemExit):
        sn.walk(["g1"], tok, sets, target=1)


def test_frame_excludes_special_and_listed():
    tok = {"<pad>": 0, "<cls>": 1, "A": 2, "B": 3}
    assert sn.frame(tok, {"B"}) == ["A"]


# ---------------------------------------------------------------- run order
def test_run_order_skips_controlless_strata_and_dedupes():
    strata = {"B00": {"status": "eligible", "controls": [f"c{i:02d}" for i in range(20)], "members": ["p1"]},
              "B01": {"status": "eligible", "controls": [f"c{i:02d}" for i in range(5, 25)], "members": ["p2"]},
              "B02": {"status": "NOT_ESTIMABLE_CONTROLS", "controls": None, "members": ["p3"]}}
    elig = [{"ensembl_id": "p1", "eligible": "True", "stratum": "B00"},
            {"ensembl_id": "p2", "eligible": "True", "stratum": "B01"},
            {"ensembl_id": "p3", "eligible": "True", "stratum": "B02"},
            {"ensembl_id": "p4", "eligible": "False", "stratum": ""}]
    null = [f"n{i:03d}" for i in range(100)]
    genes, panel, controls, spot = ro.order(null, strata, elig)
    assert panel == ["p1", "p2"] and len(controls) == 25 and genes[:100] == null
    assert len(genes) == len(set(genes)) == 127
    assert spot[:5] == ["n019", "n039", "n059", "n079", "n099"]


def test_run_order_refuses_overlap():
    strata = {"B00": {"status": "eligible", "controls": ["n000"] + [f"c{i}" for i in range(19)], "members": ["p1"]}}
    with pytest.raises(SystemExit):
        ro.order(["n000"], strata, [{"ensembl_id": "p1", "eligible": "True", "stratum": "B00"}])


# ---------------------------------------------------------------- H2c
def _luad(n_ref=6):
    rows = []
    for i in range(n_ref):
        rows.append({"panel": "B", "gene": f"g{i}", "symbol": f"S{i}", "status": "T_CELL_SIGNAL_TOWARD",
                     "del_median": 0.01, "ovx_median": -0.01})
    rows.append({"panel": "B", "gene": "gx", "symbol": "SX", "status": "OPEN", "del_median": 0.0, "ovx_median": 0.0})
    return rows


def test_binom_upper_exact():
    assert cp.binom_upper(5, 5) == Fraction(1, 32)
    assert cp.binom_upper(0, 5) == 1


def test_h2c_holds_when_all_agree():
    e2 = [{"panel": "B", "gene": f"g{i}", "status": "OPEN", "del_median": 0.002, "ovx_median": 0.001} for i in range(6)]
    r = cp.h2c(_luad(), e2)
    assert r["n_tested"] == 6 and r["del_agree"] == 6 and r["reading"] == "pattern_holds"
    assert r["p_exact"] == "1/64" and r["ovx_agree_descriptive"] == 0


def test_h2c_not_replicated_and_zero_counts_as_disagree():
    e2 = [{"panel": "B", "gene": f"g{i}", "status": "OPEN", "del_median": (0.0 if i < 2 else 0.002)} for i in range(6)]
    r = cp.h2c(_luad(), e2)
    assert r["del_agree"] == 4 and r["reading"] == "pattern_not_replicated"


def test_h2c_not_testable_below_five():
    e2 = [{"panel": "B", "gene": f"g{i}", "status": "OPEN", "del_median": 0.002} for i in range(4)]
    e2 += [{"panel": "B", "gene": "g4", "status": "NOT_RUN"}, {"panel": "B", "gene": "g5", "status": "NOT_RUN"}]
    r = cp.h2c(_luad(), e2)
    assert r["n_tested"] == 4 and r["reading"] == "not_testable" and r["p_one_sided"] is None


def test_mmr_split():
    pd_ = {"g0|delete": {"a": 1.0, "b": 3.0, "c": -1.0}}
    out = cp.mmr_split(pd_, {"a": "MMRd", "b": "MMRd", "c": "MMRp"}, ["g0"])
    assert out["g0"]["MMRd"] == {"n": 2, "median": 2.0} and out["g0"]["MMRp"]["median"] == -1.0
