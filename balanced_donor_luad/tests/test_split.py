"""Tests for make_split: invariants hold on a realistic study mix, the
assignment is deterministic, and check() FAILS on each broken split."""
import copy
import os
import sys

import pandas as pd
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
import make_split as ms  # noqa: E402


@pytest.fixture
def donor_study():
    # the real study mix: 23 / 10 / 5 / 3 / 2
    rows = [(f"{s}_{i}", s) for s, n in [("Leader", 23), ("Kim", 10), ("He", 5), ("Lamb", 3), ("Laugh", 2)]
            for i in range(n)]
    return pd.Series({d: s for d, s in rows})


def test_invariants_and_balance(donor_study):
    split = ms.build(donor_study)
    assert ms.check(split)
    assert sorted(len(f["test"]) for f in split["folds"]) == [8, 8, 9, 9, 9]
    for f in split["folds"]:
        assert len(f["eval"]) == 4
        # the 23-donor study is spread: 4 or 5 per fold, never piled into one
        assert sum(d.startswith("Leader") for d in f["test"]) in (4, 5)


def test_deterministic_and_order_independent(donor_study):
    a = ms.build(donor_study)
    b = ms.build(donor_study.sample(frac=1, random_state=3))
    assert a["donor_test_fold"] == b["donor_test_fold"]
    assert [f["eval"] for f in a["folds"]] == [f["eval"] for f in b["folds"]]


def test_check_fails_on_donor_tested_twice(donor_study):
    s = ms.build(donor_study)
    bad = copy.deepcopy(s)
    bad["folds"][1]["test"].append(bad["folds"][0]["test"][0])
    with pytest.raises(ValueError):
        ms.check(bad)


def test_check_fails_on_eval_test_overlap(donor_study):
    s = ms.build(donor_study)
    bad = copy.deepcopy(s)
    bad["folds"][0]["eval"][0] = bad["folds"][0]["test"][0]
    with pytest.raises(ValueError, match="overlap"):
        ms.check(bad)


def test_check_fails_on_unbalanced_folds(donor_study):
    s = ms.build(donor_study)
    bad = copy.deepcopy(s)
    moved = bad["folds"][0]["test"][:3]
    bad["folds"][0]["test"] = bad["folds"][0]["test"][3:]
    bad["folds"][1]["test"] += moved
    for f in bad["folds"]:
        f["train"] = sorted(set(f["train"]) - set(moved)) if f["fold"] == 1 else sorted(set(f["train"]) | (set(moved) if f["fold"] == 0 else set()))
    with pytest.raises(ValueError):
        ms.check(bad)
