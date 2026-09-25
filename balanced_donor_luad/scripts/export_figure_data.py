"""Figure inputs that need the Phase 6 tree: the control-gene null cloud for the dose-concordance map.

For each control c in each eligible stratum: a_op(c,d) = s_op(c,d) - median over the OTHER controls of the
stratum (same 3g.1 rules), median over donors. Descriptive only: controls are not outcomes, no test is run.
Also records, per panel gene, the donor count and the per-fold delete medians already in the rows.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyse as an  # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6-root", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--ovx-index", required=True)
    ap.add_argument("--rules", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    design = an.Design(**json.load(open(a.design)))
    rules = json.load(open(a.rules))
    ix = json.load(open(a.ovx_index))
    positions = {tuple(k.split("|")): v for k, v in ix["positions"].items()}
    calls = an.load_all(a.phase6_root, an.genes_to_load(design), design.donors, positions, ix["n_total"])
    cloud = {}
    for sname, st in design.strata.items():
        if st["status"] != "eligible":
            continue
        for c in st["controls"]:
            rec = {"stratum": sname}
            for op in an.OPS:
                vals = {x: an.donor_means(calls, x, op, an.GOAL, rules["control_min_cells"]) for x in st["controls"]}
                adj = []
                for d, (s, _) in vals[c].items():
                    others = [vals[x][d][0] for x in st["controls"] if x != c and d in vals[x]]
                    if len(others) >= rules["min_controls_per_donor"]:
                        adj.append(s - float(np.median(others)))
                rec[op] = {"median": float(np.median(adj)) if adj else None, "n": len(adj)}
            cloud[f"{sname}|{c}"] = rec
    json.dump({"about": __doc__.strip(), "control_cloud": cloud}, open(a.out, "w"), indent=1)
    print(len(cloud), "control entries")
    return 0


if __name__ == "__main__":
    sys.exit(main())
