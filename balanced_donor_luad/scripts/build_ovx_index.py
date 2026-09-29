"""Amendment 3h: token-positive positions for every overexpress call, plus the two COUNT checks.

  (i)  delete: the replayed token-positive count equals each delete pickle's per-cell length
  (ii) overexpress: the reconstructed positive count equals the marker's n_token_cells
These are count checks only; identity is tested separately on the GPU (identity_check.py).
Reads markers, pickle LENGTHS (never values), ISP inputs and the token dictionary.
Writes ovx_index.json ({"n_total", "positions", "ties"}) and refuses (exit 1) on any mismatch.
"""
import argparse
import collections
import glob
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyse as an  # noqa: E402
import isp_order as io  # noqa: E402


def pickle_len(root, op, gene, donor):
    fs = glob.glob(os.path.join(root, op, gene, f"in_silico_{op}_{donor.replace('/', '_')}_cell_embs_dict_*_raw.pickle"))
    if len(fs) != 1:
        return None
    lens = {len(next(iter(v.values()))) for v in pickle.load(open(fs[0], "rb")).values()}
    return lens.pop() if len(lens) == 1 else -1


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6-root", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--manifest", required=True, help="goals/donor_manifest.json (isp_input per donor)")
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    design = an.Design(**json.load(open(a.design)))
    genes = an.genes_to_load(design)
    tok = pickle.load(open(a.token_dict, "rb"))
    manifest = json.load(open(a.manifest))
    orders = io.donor_orders(manifest, design.donors)

    positions, fails, n_checked, n_tie = {}, [], collections.Counter(), 0
    for d, (order, sets, lengths) in orders.items():
        for g in genes:
            t = tok[g]
            pos = io.positive_positions(order, sets, t)
            k = len(pos)
            for op in an.OPS:
                rec = json.load(open(os.path.join(a.phase6_root, op, g, f"{d.replace('/', '_')}.complete.json")))
                if rec["status"] == "no_token_cells":
                    if k != 0:
                        fails.append(f"{op}/{g}/{d}: marker no_token_cells but replay finds {k}")
                    continue
                L = pickle_len(a.phase6_root, op, g, d)
                if op == "overexpress":
                    n_checked["ii"] += 1
                    if k != rec["n_token_cells"]:
                        fails.append(f"(ii) {g}/{d}: reconstructed {k} != marker {rec['n_token_cells']}")
                    if L != len(order):
                        fails.append(f"{op}/{g}/{d}: pickle length {L} != donor cells {len(order)}")
                else:
                    n_checked["i"] += 1
                    if L != k:
                        fails.append(f"(i) {g}/{d}: delete pickle length {L} != replayed positive count {k}")
            if k:
                positions[f"{g}|{d}"] = pos
                n_tie += io.has_mixed_tie(lengths, sets, t)
    out = {"n_total": {d: len(o) for d, (o, _, _) in orders.items()}, "positions": positions,
           "checks": {"i_delete_lengths": n_checked["i"], "ii_overexpress_counts": n_checked["ii"],
                      "failures": len(fails), "pairs_with_positive": len(positions),
                      "pairs_with_mixed_length_tie": n_tie}}
    if fails:
        print(f"REFUSED: {len(fails)} count-check failures; first: {fails[:5]}", file=sys.stderr)
        return 1
    json.dump(out, open(a.out, "w"), sort_keys=True)
    print(json.dumps(out["checks"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
