"""3h.6 selection rule: categories in order, distinct pairs, ID tie-break, empty never substituted."""
import os
import sys

HERE = os.path.dirname(__file__)
sys.path[:0] = [os.path.join(HERE, "..", "scripts")]
import select_identity_pairs as sp  # noqa: E402


def c(g, d, k, m, host="h1"):
    return {"gene": g, "donor": d, "host": host, "k": k, "mixed_tie_positives": m}


def test_categories_distinct_and_tie_broken_by_id():
    cands = [c("G2", "D1", 50, 9), c("G1", "D1", 50, 9), c("G3", "D1", 60, 0), c("G4", "D1", 40, 0),
             c("G5", "D2", 11, 3), c("G6", "D2", 95, 4)]
    got = {p["category"]: (p["gene"], p["donor"]) for p in sp.select(cands)}
    assert got == {"tie_rich": ("G1", "D1"), "tie_free": ("G3", "D1"), "low_k": ("G5", "D2"), "high_k": ("G6", "D2")}


def test_pair_taken_by_an_earlier_category_is_not_reused():
    cands = [c("G1", "D1", 11, 9), c("G2", "D1", 12, 1), c("G3", "D1", 30, 0)]
    got = {p["category"]: p.get("gene") for p in sp.select(cands)}
    assert got["tie_rich"] == "G1" and got["low_k"] == "G2"


def test_empty_category_reported_not_substituted():
    got = {p["category"]: p for p in sp.select([c("G1", "D1", 50, 2)])}
    assert got["tie_free"].get("empty") and got["high_k"].get("empty") and got["low_k"].get("empty")


def test_mixed_tie_positive_count():
    assert sp.mixed_tie_positives([5, 5, 7, 7], [{1}, set(), {1}, {1}], 1) == 1
