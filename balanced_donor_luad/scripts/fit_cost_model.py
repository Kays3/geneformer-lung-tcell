"""Amendment 5 (ISP-STD-1 E.5): fixed + per-cell GPU-arm cost model, fit from Phase 6's own
recorded per-call data. Least squares: seconds ~ a + b * n_token_cells, separately per operation.
Reads run_isp.py's own "seconds" and "n_token_cells" fields from *.complete.json markers -- never
wall clock, never an inferred rate. CPU only.
"""
import argparse
import glob
import json

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6-roots", nargs="+", required=True, help="one or more hosts' phase6 output roots")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    xs = {"delete": [], "overexpress": []}
    ys = {"delete": [], "overexpress": []}
    for root in a.phase6_roots:
        for f in glob.glob(f"{root}/*/*/*.complete.json"):
            r = json.load(open(f))
            op = r.get("op")
            if op in xs and "seconds" in r and "n_token_cells" in r:
                xs[op].append(r["n_token_cells"])
                ys[op].append(r["seconds"])

    fit = {}
    for op in ("delete", "overexpress"):
        x = np.array(xs[op], dtype=float)
        y = np.array(ys[op], dtype=float)
        A = np.vstack([x, np.ones_like(x)]).T
        b, aa = np.linalg.lstsq(A, y, rcond=None)[0]
        pred = aa + b * x
        resid = y - pred
        fit[op] = {
            "n_calls": int(len(x)), "intercept_a": float(aa), "slope_b_per_cell": float(b),
            "mean_n_cells": float(x.mean()), "max_n_cells": float(x.max()),
            "rmse": float(np.sqrt((resid ** 2).mean())), "mean_resid": float(resid.mean()),
            "p95_abs_resid": float(np.percentile(np.abs(resid), 95)), "max_resid": float(resid.max()),
            "sum_actual_gpu_h": float(y.sum() / 3600), "sum_predicted_gpu_h": float(pred.sum() / 3600),
        }
    json.dump({"about": __doc__.strip(), "sources": a.phase6_roots, "fit": fit}, open(a.out, "w"), indent=1)
    print(json.dumps(fit, indent=1))


if __name__ == "__main__":
    main()
