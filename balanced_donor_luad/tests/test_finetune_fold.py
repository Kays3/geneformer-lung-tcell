"""finetune_fold's metric extraction: accepts validate()'s one-element per-split
lists (the 2026-09-24 Phase 4 crash) and FAILS rather than silently picking a
value when there is more than one split."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))
from finetune_fold import single_split_value  # noqa: E402


def test_one_element_list_as_validate_returns():
    assert single_split_value({"macro_f1": [0.71]}, "macro_f1") == 0.71


def test_scalar_and_missing():
    assert single_split_value({"acc": 0.5}, "acc") == 0.5
    assert single_split_value({}, "acc") is None


def test_fails_on_more_than_one_split():
    with pytest.raises(ValueError):
        single_split_value({"acc": [0.6, 0.7]}, "acc")
