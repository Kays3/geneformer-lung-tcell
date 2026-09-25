"""Amendment 3h.6 pair selection for the GPU identity check, applied exactly as registered.

Candidates: (control gene, donor), 10 <= k <= 99. Host = the host that ran the gene in Phase 6.
"Mixed-tie positives" = token-positive cells whose length is shared with a token-negative cell.
Per host, four distinct pairs in this order (ties broken by ascending (Ensembl ID, donor)):
  tie_rich   most mixed-tie positives
  tie_free   zero mixed-tie positives, largest k
  low_k      10 <= k <= 12, >= 1 mixed-tie positive, most mixed-tie positives
  high_k     90 <= k <= 99, >= 1 mixed-tie positive, most mixed-tie positives
An empty category is reported empty, never substituted. Reads tokens and lengths only.
"""
import argparse
import csv
import json
import os
import pickle
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyse as an  # noqa: E402
import isp_order as io  # noqa: E402

RULES = [
    ("tie_rich", lambda c: True, lambda c: (-c["mixed_tie_positives"], c["gene"], c["donor"])),
    ("tie_free", lambda c: c["mixed_tie_positives"] == 0, lambda c: (-c["k"], c["gene"], c["donor"])),
    ("low_k", lambda c: 10 <= c["k"] <= 12 and c["mixed_tie_positives"] >= 1,
     lambda c: (-c["mixed_tie_positives"], c["gene"], c["donor"])),
    ("high_k", lambda c: 90 <= c["k"] <= 99 and c["mixed_tie_positives"] >= 1,
     lambda c: (-c["mixed_tie_positives"], c["gene"], c["donor"])),
]


def mixed_tie_positives(lengths, sets, token):
    neg_lengths = {L for L, s in zip(lengths, sets) if token not in s}
    return sum(1 for L, s in zip(lengths, sets) if token in s and L in neg_lengths)


def candidates(design, host_of, tok, orders):
    controls = sorted({c for s in design.strata.values() for c in (s.get("controls") or [])})
    out = []
    for d, (order, sets, lengths) in orders.items():
        for g in controls:
            k = sum(tok[g] in s for s in sets)
            if 10 <= k <= 99:
                out.append({"gene": g, "donor": d, "host": host_of[g], "k": k,
                            "mixed_tie_positives": mixed_tie_positives(lengths, sets, tok[g])})
    return out


def select(cands):
    chosen = []
    for host in sorted({c["host"] for c in cands}):
        taken = set()
        for name, keep, key in RULES:
            pool = sorted((c for c in cands if c["host"] == host and keep(c)
                           and (c["gene"], c["donor"]) not in taken), key=key)
            if pool:
                taken.add((pool[0]["gene"], pool[0]["donor"]))
                chosen.append({"category": name, **pool[0]})
            else:
                chosen.append({"category": name, "host": host, "empty": True})
    return chosen


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--design", required=True)
    ap.add_argument("--host-map", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    design = an.Design(**json.load(open(a.design)))
    host_of = {r["ensembl_id"]: r["host"] for r in csv.DictReader(open(a.host_map))}
    tok = pickle.load(open(a.token_dict, "rb"))
    orders = io.donor_orders(json.load(open(a.manifest)), design.donors)
    cands = candidates(design, host_of, tok, orders)
    pairs = select(cands)
    json.dump({"rule": "Amendment 3h.6", "n_candidates": len(cands), "pairs": pairs}, open(a.out, "w"), indent=1)
    for p in pairs:
        print(json.dumps(p))
    return 0


if __name__ == "__main__":
    sys.exit(main())
