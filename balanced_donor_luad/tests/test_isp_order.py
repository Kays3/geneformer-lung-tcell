"""isp_order replays Geneformer's start-state filter + length sort; ties must keep dataset order."""
import os
import sys

from datasets import Dataset

HERE = os.path.dirname(__file__)
sys.path[:0] = [os.path.join(HERE, "..", "scripts")]
import isp_order as io  # noqa: E402

T = 7


def ds(cells):
    """cells: list of (length, has_token, origin)."""
    rows = {"input_ids": [], "length": [], "origin": []}
    for i, (L, pos, org) in enumerate(cells):
        ids = [T] + [100 + j for j in range(L - 1)] if pos else [100 + j for j in range(L)]
        rows["input_ids"].append(ids); rows["length"].append(L); rows["origin"].append(org)
    d = Dataset.from_dict(rows)
    return d.add_column("_orig_idx", list(range(len(d))))


def test_order_is_length_descending_with_ties_in_dataset_order():
    d = ds([(5, 0, "tumor_primary"), (9, 1, "tumor_primary"), (5, 1, "tumor_primary"), (7, 0, "tumor_primary"),
            (5, 0, "tumor_primary")])
    assert io.replay_order(d) == [1, 3, 0, 2, 4]            # the three length-5 cells keep 0, 2, 4


def test_start_state_filter_drops_other_origins():
    d = ds([(5, 1, "normal_adjacent"), (6, 1, "tumor_primary")])
    assert io.replay_order(d) == [1]


def test_positive_positions_across_a_mixed_tie():
    d = ds([(5, 0, "tumor_primary"), (9, 1, "tumor_primary"), (5, 1, "tumor_primary"), (7, 0, "tumor_primary")])
    order = io.replay_order(d)                              # [1, 3, 0, 2]
    sets = [set(x) for x in d["input_ids"]]
    assert io.positive_positions(order, sets, T) == [0, 3]
    assert io.has_mixed_tie(d["length"], sets, T) is True


def test_no_mixed_tie_when_tied_cells_agree():
    d = ds([(5, 1, "tumor_primary"), (5, 1, "tumor_primary"), (7, 0, "tumor_primary")])
    assert io.has_mixed_tie(d["length"], [set(x) for x in d["input_ids"]], T) is False
