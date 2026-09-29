"""Which positions of an OVERexpress raw pickle belong to token-positive cells (Amendment 3h).

Geneformer's InSilicoPerturber (vendored f45a6c7) filters to token-positive cells only for delete;
overexpress with a gene list perturbs every start-state cell. Both then pass through
perturber_utils.downsample_and_sort (max_ncells=None -> Dataset.sort("length", reverse=True)),
so per-cell output order is the length-descending order of the input, ties kept in dataset order
(pyarrow sort_indices is stable). We replay that exact call on the same input rather than
re-implementing it, and then select the token-positive positions.

No shift value is read here: only input_ids, length and the start-state column.
"""
from datasets import load_from_disk

STATE_KEY, START_STATE = "origin", "tumor_primary"


def replay_order(ds):
    """Original row indices in the order Geneformer processed them."""
    ds = ds.filter(lambda r: r[STATE_KEY] in [START_STATE], num_proc=None)      # filter_data_by_start_state
    ds = ds.select(range(len(ds)))                                               # downsample_and_sort, no cap
    return list(ds.sort("length", reverse=True)["_orig_idx"])


def donor_orders(manifest, donors):
    """{donor: (order, token_sets, lengths)} from each donor's ISP input."""
    out = {}
    for d in donors:
        ds = load_from_disk(manifest[d]["isp_input"])
        ds = ds.add_column("_orig_idx", list(range(len(ds))))
        out[d] = (replay_order(ds), [set(x) for x in ds["input_ids"]], list(ds["length"]))
    return out


def has_mixed_tie(lengths, token_sets, token):
    """True if some length is shared by a token-positive and a token-negative cell. Such ties are
    resolved only by the sort's stability, which the GPU identity check tests (3h)."""
    seen = {}
    for L, s in zip(lengths, token_sets):
        seen.setdefault(L, set()).add(token in s)
    return any(len(v) == 2 for v in seen.values())


def positive_positions(order, token_sets, token):
    """Positions in the processed order whose cell contains `token`."""
    return [p for p, i in enumerate(order) if token in token_sets[i]]


def delete_positions_check(order, token_sets, token):
    """Delete output has one entry per token-positive cell: its length is the must-match control."""
    return sum(token in token_sets[i] for i in order)


