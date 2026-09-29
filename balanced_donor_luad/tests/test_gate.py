"""Tests for the classifier gate: exact sign test vs enumeration, and the gate
FAILS on each input it exists to catch (low pooled BA, too few donors above
chance, significance in the wrong direction)."""
import itertools
import os
import sys
from fractions import Fraction

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import classifier_gate as cg  # noqa: E402


def brute_sign_p(n_pos, n_neg):
    n = n_pos + n_neg
    obs = abs(n_pos - n_neg)
    hits = sum(1 for s in itertools.product([0, 1], repeat=n) if abs(2 * sum(s) - n) >= obs)
    return Fraction(hits, 2 ** n)


def test_sign_test_matches_enumeration():
    for n_pos in range(0, 9):
        for n_neg in range(0, 9 - n_pos):
            assert cg.sign_test_two_sided(n_pos, n_neg) == brute_sign_p(n_pos, n_neg)


def cells(donor_ba, n=100):
    """Build per-cell predictions giving each donor the requested per-class accuracy."""
    rows = []
    for d, acc in donor_ba.items():
        k = int(round(acc * n))
        for cls in ("tumor_primary", "normal_adjacent"):
            other = "normal_adjacent" if cls == "tumor_primary" else "tumor_primary"
            rows += [{"donor": d, "y_true": cls, "y_pred": cls, "fold": 0}] * k
            rows += [{"donor": d, "y_true": cls, "y_pred": other, "fold": 0}] * (n - k)
    return pd.DataFrame(rows)


def test_passes_when_clearly_above_chance():
    r = cg.evaluate(cells({f"d{i}": 0.8 for i in range(12)}))
    assert r["PASS"] and abs(r["pooled_balanced_accuracy"] - 0.8) < 1e-9


def test_fails_on_low_pooled_ba_even_if_every_donor_above_half():
    r = cg.evaluate(cells({f"d{i}": 0.55 for i in range(20)}))
    assert r["donors_above_0.5"] == 20 and not r["PASS"]


def test_fails_when_too_few_donors_for_sign_test():
    r = cg.evaluate(cells({f"d{i}": 0.9 for i in range(5)}))     # min p = 2/32 > 0.05
    assert not r["PASS"]


def test_direction_guard_alone_blocks_wrong_direction(monkeypatch):
    # With 12 donors below chance the pooled BA can never reach 0.60, so the
    # pooled bar would mask the direction guard. Set the bar to 0 to isolate it:
    # without `n_pos > n_neg` this case would PASS.
    monkeypatch.setattr(cg, "GATE_BA", 0.0)
    r = cg.evaluate(cells({f"d{i}": 0.2 for i in range(12)}))
    assert r["sign_test_p_float"] <= 0.05 and r["donors_below_0.5"] == 12
    assert not r["PASS"]


def test_3e_margins_flag_thin_pooled_margin_but_not_decision():
    r = cg.evaluate(cells({f"d{i}": 0.62 for i in range(12)}))
    assert r["PASS"] and r["pooled_within_nondet_of_threshold"] and not r["clean_result"]


def test_3e_clean_when_far_from_both_thresholds():
    r = cg.evaluate(cells({f"d{i}": 0.8 for i in range(12)}))
    assert r["PASS"] and r["n_flippable_donors"] == 0 and r["clean_result"]


def test_3e_flippable_donors_can_flip_sign_test():
    # 7 firm above, 5 flippable just above 0.5: best case 12-0 passes, worst 7-5 does not.
    donors = {**{f"a{i}": 0.8 for i in range(7)}, **{f"f{i}": 0.53 for i in range(5)}}
    r = cg.evaluate(cells(donors))
    assert r["n_flippable_donors"] == 5 and r["sign_test_within_nondet"] and not r["clean_result"]
