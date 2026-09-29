"""Amendment 4 primary analysis: genome-wide null distribution of delete/overexpress concordance.

CPU only. Reads the completed phase8 markers (same shape as Phase 6's), computes each drawn gene's raw
donor-level median shift per operation (no stratum control-adjustment -- Amendment 4 s.4.4), then the
one-sided permutation test on Spearman rho (H1: rho < 0), leave-one-gene-out, and a gene-level bootstrap.

Never reads Shift_to_goal_end; reads the raw per-cell pickles Amendment 3g/3h established, via the same
donor_means() helper as analyse.py, imported unchanged.
"""
import argparse
import glob
import json
import os
import sys

import numpy as np
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyse as an  # noqa: E402

OPS = an.OPS
CONTROL_MIN_CELLS = 10
D_MIN = 10
N_PERM = 100_000
N_BOOT = 10_000
ALPHA = 0.05
SEED = 20260929


def gene_medians(phase6_root, gene, donors, positions, n_total, control_min_cells):
    calls = an.load_all(phase6_root, [gene], donors, positions, n_total)
    out = {}
    for op in OPS:
        dm = an.donor_means(calls, gene, op, an.GOAL, control_min_cells)
        estimable = {d: s for d, (s, n) in dm.items() if n >= control_min_cells}
        if len(estimable) < D_MIN:
            out[op] = None
        else:
            out[op] = {"median": float(np.median(list(estimable.values()))), "n_donors": len(estimable)}
    return out


def perm_test(x, y, rng, n_perm=N_PERM):
    obs = spearmanr(x, y).correlation
    y = np.asarray(y)
    count = sum(1 for _ in range(n_perm) if spearmanr(x, rng.permutation(y)).correlation <= obs)
    return obs, (1 + count) / (n_perm + 1)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase6-root", required=True, help="combined phase8_null output tree, both ops")
    ap.add_argument("--null-genes", required=True)
    ap.add_argument("--design", required=True, help="phase7_results/design.json, for donors + positions inputs")
    ap.add_argument("--ovx-index", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ng = json.load(open(a.null_genes))
    design = an.Design(**json.load(open(a.design)))
    ix = json.load(open(a.ovx_index))
    positions = {tuple(k.split("|")): v for k, v in ix["positions"].items()}

    rows = []
    for g in ng["draw_order"]:
        # only genes with completed markers for both ops, all donors, are eligible; anything else is
        # stopped_not_analysed (the budget ceiling stopped the run before this gene, or it never started).
        marker_glob = os.path.join(a.phase6_root, "*", g, "*.complete.json")
        if not glob.glob(marker_glob):
            rows.append({"gene": g, "status": "stopped_not_analysed"})
            continue
        m = gene_medians(a.phase6_root, g, design.donors, positions, ix["n_total"], CONTROL_MIN_CELLS)
        if m["delete"] is None or m["overexpress"] is None:
            rows.append({"gene": g, "status": "not_estimable", "detail": m})
            continue
        rows.append({"gene": g, "status": "estimable", "delete_median": m["delete"]["median"],
                     "overexpress_median": m["overexpress"]["median"],
                     "delete_n_donors": m["delete"]["n_donors"], "overexpress_n_donors": m["overexpress"]["n_donors"]})

    estimable = [r for r in rows if r["status"] == "estimable"]
    n = len(estimable)
    result = {"about": __doc__.strip(), "n_drawn": len(ng["draw_order"]), "n_estimable": n,
              "n_stopped_not_analysed": sum(r["status"] == "stopped_not_analysed" for r in rows),
              "n_not_estimable": sum(r["status"] == "not_estimable" for r in rows), "rows": rows}

    if n < 10:
        result["primary_status"] = "not_estimable"
    else:
        x = np.array([r["delete_median"] for r in estimable])
        y = np.array([r["overexpress_median"] for r in estimable])
        rng = np.random.default_rng(SEED)
        rho, p = perm_test(x, y, rng)
        # leave-one-out
        loo_signs = []
        for i in range(n):
            keep = [j for j in range(n) if j != i]
            loo_signs.append(np.sign(spearmanr(x[keep], y[keep]).correlation))
        loo_stable = all(s == np.sign(rho) for s in loo_signs)
        # bootstrap over genes
        boot_hits = 0
        for _ in range(N_BOOT):
            idx = rng.integers(0, n, n)
            rb = spearmanr(x[idx], y[idx]).correlation
            if np.sign(rb) == np.sign(rho):
                _, pb = perm_test(x[idx], y[idx], rng, n_perm=2000)
                if pb <= ALPHA:
                    boot_hits += 1
        boot_frac = boot_hits / N_BOOT
        primary_pass = rho < 0 and p <= ALPHA
        if primary_pass and loo_stable and boot_frac >= 0.95:
            status = "negative"
        elif rho > 0 and p <= ALPHA:
            status = "opposite_direction"
        elif primary_pass:
            status = "control_draw_sensitive_open"
        else:
            status = "primary_not_met"
        result.update({"primary_status": status, "rho": float(rho), "p_permutation": float(p),
                       "n_perm": N_PERM, "loo_stable": bool(loo_stable), "bootstrap_fraction_stable": boot_frac,
                       "n_bootstrap": N_BOOT, "alpha": ALPHA, "sidedness": "one-sided, H1: rho<0"})

    json.dump(result, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items() if k != "rows"}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
