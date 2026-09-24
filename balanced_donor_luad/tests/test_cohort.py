"""Tests for the cohort draw and the QC assertions.

Each QC test breaks the cohort in one specific way and requires the matching
check to FAIL; a suite that only ever sees valid cohorts cannot tell a working
check from one that always passes. Donors are given UNEQUAL eligible counts so
the equal-cap result is produced by the draw, not by the fixture.
"""
import os
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from build_cohort import draw  # noqa: E402
from qc_cohort import checks  # noqa: E402

ORIGINS = ("tumor_primary", "normal_adjacent")


def eligible(n_donors=12, seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for d in range(n_donors):
        for o in ORIGINS:
            n = int(rng.integers(150, 2000))  # unequal on purpose
            for i in range(n):
                rows.append({"cell_id": f"D{d}_{o}_{i}", "donor_id": f"D{d}", "origin": o,
                             "study": "S", "assay": "10x 3' v2", "cell_type_major": "T cell CD4"})
    return pd.DataFrame(rows)


@pytest.fixture(scope="module")
def cohort():
    c = draw(eligible())
    return c[c.in_pool300].reset_index(drop=True)


def failed(c):
    return {name for name, _, _, ok in checks(c) if not ok}


def test_valid_cohort_passes_everything(cohort):
    assert failed(cohort) == set()


def test_draw_is_independent_of_row_order():
    e = eligible(n_donors=3, seed=1)
    a = draw(e).set_index("cell_id").draw_rank
    b = draw(e.sample(frac=1, random_state=7)).set_index("cell_id").draw_rank
    assert a.sort_index().equals(b.sort_index())


def test_draw_depends_on_seed():
    e = eligible(n_donors=2, seed=2)
    a = draw(e, seed=1).set_index("cell_id").draw_rank.sort_index()
    b = draw(e, seed=2).set_index("cell_id").draw_rank.sort_index()
    assert not a.equals(b)


def test_one_missing_analysis_cell_fires_a2_a3_a4(cohort):
    c = cohort.copy()
    i = c.index[(c.donor_id == "D0") & (c.origin == ORIGINS[0]) & c.in_analysis100][0]
    c.loc[i, "in_analysis100"] = False
    f = failed(c)
    assert f"A2 [{ORIGINS[0]}] analysis cells per donor (min, max)" in f
    assert f"A4 [{ORIGINS[0]}] analysis max/min donor cells" in f
    assert f"A3 [{ORIGINS[0]}] largest donor share of analysis cells" in f


def test_unpaired_donor_fires_f3(cohort):
    c = cohort[~((cohort.donor_id == "D1") & (cohort.origin == ORIGINS[1]))]
    assert "F3 every analysis donor present on both sides (paired)" in failed(c)


def test_too_few_donors_fires_a1():
    c = draw(eligible(n_donors=11))
    assert f"A1 [{ORIGINS[0]}] qualifying donors" in failed(c[c.in_pool300])


def test_mixed_assay_fires_f4(cohort):
    c = cohort.copy()
    c.loc[(c.donor_id == "D2") & (c.origin == ORIGINS[0]), "assay"] = "10x 5' v1"
    assert "F4 assay set identical on both sides per donor" in failed(c)


def test_pool_over_cap_fires(cohort):
    # Pick a donor already AT the cap; one below it would absorb the extra
    # cell and the check would rightly pass, making this test uninformative.
    c = cohort.copy()
    sizes = c[c.origin == ORIGINS[0]].groupby("donor_id").size()
    donor = sizes[sizes == 300].index[0]
    extra = c[(c.donor_id == donor) & (c.origin == ORIGINS[0]) & ~c.in_analysis100].iloc[:1].copy()
    extra["cell_id"] = "extra"
    c = pd.concat([c, extra])
    assert f"A4 [{ORIGINS[0]}] training-pool cells per donor capped" in failed(c)


def test_duplicate_cell_fires(cohort):
    c = pd.concat([cohort, cohort.iloc[:1]])
    assert "no duplicate cell_id" in failed(c)
