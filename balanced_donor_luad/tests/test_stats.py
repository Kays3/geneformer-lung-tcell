"""Tests for stats_core. Each check is compared against an independent
computation (brute-force enumeration, a published table value, or a hand
example), and each rule is shown to FAIL on the input it exists to catch."""
import itertools
import os
import sys
from fractions import Fraction

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import stats_core as sc  # noqa: E402


def brute_force_p(x):
    """Two-sided exact p by enumerating all 2^n sign flips of the midranks."""
    x = np.asarray(x, dtype=float)
    x = x[x != 0]
    n = len(x)
    if n == 0:
        return Fraction(1)
    r2 = sc._doubled_midranks(np.abs(x))
    total2 = int(r2.sum())
    obs = abs(2 * int(r2[x > 0].sum()) - total2)
    hits = 0
    for signs in itertools.product([0, 1], repeat=n):
        w2 = sum(int(r) for r, s in zip(r2, signs) if s)
        hits += abs(2 * w2 - total2) >= obs
    return Fraction(hits, 2 ** n)


@pytest.mark.parametrize("seed", range(40))
def test_wilcoxon_exact_matches_brute_force_including_ties_and_zeros(seed):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(1, 11))
    x = rng.choice([-3, -2, -1, 0, 1, 2, 3, 2.5, -2.5], size=n)  # forces ties and zeros
    assert sc.wilcoxon_exact(x)["p"] == brute_force_p(x)


def test_all_positive_untied_reaches_min_attainable_p():
    for n in range(1, 15):
        x = np.arange(1, n + 1, dtype=float)
        assert sc.wilcoxon_exact(x)["p"] == sc.min_attainable_p(n) == Fraction(2, 2 ** n)


def test_wilcoxon_known_value_against_table():
    # n=10, W+ = 8 is the classical two-sided 0.05 critical boundary: P(W<=8) = 25/1024
    x = [1, 2, 3, -4, -5, -6, -7, -8, -9, -10]   # W+ = 6  (below 8)
    assert sc.wilcoxon_exact(x)["p"] == 2 * Fraction(sum(1 for s in itertools.product([0, 1], repeat=10)
                                                          if sum(r for r, b in zip(range(1, 11), s) if b) <= 6), 1024)


def test_registered_d_min_by_scanning():
    assert sc.smallest_d_passing(Fraction(5, 100) / 15) == 10   # Panel A
    assert sc.smallest_d_passing(Fraction(5, 100) / 36) == 11   # Panel B
    assert sc.smallest_d_passing(Fraction(5, 100) / 51) == 11   # combined (not used)
    # the neighbours really fail / pass
    assert sc.min_attainable_p(9) > Fraction(5, 100) / 15 >= sc.min_attainable_p(10)


def test_walsh_ci_n10_uses_classical_critical_value():
    x = np.arange(1, 11, dtype=float)
    med, lo, hi, cov = sc.walsh_ci(x)
    walsh = sorted((a + b) / 2 for a, b in itertools.combinations_with_replacement(x, 2))
    assert (lo, hi) == (walsh[8], walsh[-9])          # c = 8 -> 9th from each end
    assert cov == 1 - 2 * Fraction(25, 1024)


def test_walsh_ci_too_small_n_has_no_interval():
    assert sc.walsh_ci([1.0, 2.0, 3.0])[1] is None     # n=3 cannot reach 95%


def test_holm_hand_example_and_monotone():
    p = [Fraction(1, 100), Fraction(4, 100), Fraction(3, 100), Fraction(5, 1000)]
    assert sc.holm(p) == [Fraction(3, 100), Fraction(6, 100), Fraction(6, 100), Fraction(2, 100)]
    # an unadjusted p that is not multiplied would be caught here
    assert sc.holm([Fraction(1, 100), Fraction(1, 100)]) == [Fraction(2, 100), Fraction(2, 100)]


@pytest.mark.parametrize("args,expected", [
    ((True, 0.1, True, -0.1), "COHERENT"),
    ((True, 0.1, True, 0.1), "INCOHERENT"),      # the S100A8/A9 signature
    ((True, -0.1, True, -0.1), "INCOHERENT"),
    ((True, 0.1, False, 0.0), "UNRESOLVED"),
    ((False, 0.0, True, -0.1), "UNRESOLVED"),
    ((False, 0.1, False, -0.1), "NONE"),
])
def test_concordance_table(args, expected):
    assert sc.concordance(*args) == expected


@pytest.mark.parametrize("args,expected", [
    ((True, True, True, 0.1, "COHERENT", False), "REPLICATED"),
    ((True, True, True, 0.1, "COHERENT", True), "REPLICATED_AMBIENT"),
    ((True, True, True, -0.1, "COHERENT", False), "REVERSED"),
    ((True, True, True, 0.1, "INCOHERENT", False), "DOSE_INCOHERENT"),
    ((True, True, True, 0.1, "UNRESOLVED", False), "DELETION_ONLY"),
    ((True, True, False, 0.1, "UNRESOLVED", False), "OPEN"),
    ((False, True, False, 0.0, "NONE", False), "NOT_RUN"),
    ((True, False, True, 0.1, "COHERENT", False), "NOT_ESTIMABLE_CONTROLS"),
])
def test_panel_a_status_table(args, expected):
    assert sc.panel_a_status(*args) == expected


def test_deletion_only_is_never_replicated():
    assert sc.panel_a_status(True, True, True, 0.5, "UNRESOLVED", False) != "REPLICATED"


def test_panel_b_statuses():
    assert sc.panel_b_status(True, True, True, 0.1, "COHERENT") == "T_CELL_SIGNAL_TOWARD"
    assert sc.panel_b_status(True, True, True, -0.1, "COHERENT") == "T_CELL_SIGNAL_AWAY"
    assert sc.panel_b_status(True, True, True, 0.1, "INCOHERENT") == "DOSE_INCOHERENT"


def test_rows_carry_registered_fields_and_validator_refuses_stripped_qualifier():
    good = [sc.make_row("A", "g1", s) for s in sc.NON_REPLICATION] + [sc.make_row("A", "g2", "REPLICATED"),
                                                                        sc.make_row("B", "g3", "OPEN")]
    assert sc.validate_rows(good)
    stripped = dict(good[0]); stripped.pop("claim_qualifier")
    with pytest.raises(ValueError, match="claim_qualifier"):
        sc.validate_rows([stripped])
    reworded = dict(good[0]); reworded["claim_qualifier"] = "did not replicate"
    with pytest.raises(ValueError, match="claim_qualifier"):
        sc.validate_rows([reworded])


def test_lusc_rows_carry_between_donor_warrant_and_missing_warrant_fails():
    r = sc.make_row("A", "g", "REPLICATED", involves_lusc=True)
    assert r["design_warrant"] == sc.WARRANT_LUSC
    bad = dict(r); bad.pop("design_warrant")
    with pytest.raises(ValueError, match="design_warrant"):
        sc.validate_rows([bad])
