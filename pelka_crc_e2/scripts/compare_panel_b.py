"""E2 H2c (registration s.6.3): does the LUAD Panel B pattern hold in colorectal tumour T cells?

Inputs: the LUAD outcome rows (balanced_donor_luad/phase7_results/outcome_rows.json on main) and the E2
outcome rows and per-donor adjusted values written by analyse.py (unchanged) on the E2 ISP output.

Primary (registered): among the LUAD reference genes (Panel B status T_CELL_SIGNAL_TOWARD or
T_CELL_SIGNAL_AWAY; 13 genes) that are TESTED in E2 (have a control-adjusted deletion median), count those
whose E2 deletion median has the same sign as the LUAD deletion median (a zero E2 median counts as not
agreeing). Exact one-sided binomial test against 1/2, H1: agreement > 1/2. Needs n_tested >= 5 (the
smallest n at which p can reach 0.05 is 5: 1/32 = 0.031).
  pattern_holds          n >= 5 and p <= 0.05
  pattern_not_replicated n >= 5 and p > 0.05  (reported with a/n and p; never "no effect")
  not_testable           n < 5
Descriptive only (never change the reading): the same count for the overexpression median; the E2 status
of every Panel B gene next to its LUAD status; Spearman rho between LUAD and E2 deletion medians over genes
tested in both; per-gene E2 deletion medians within MMRd and within MMRp donors.
"""
import argparse
import json
from fractions import Fraction
from math import comb

import numpy as np

REFERENCE_STATUSES = {"T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY"}
MIN_N = 5
ALPHA = Fraction(5, 100)


def binom_upper(a, n):
    """P(X >= a), X ~ Binomial(n, 1/2), exact."""
    return Fraction(sum(comb(n, k) for k in range(a, n + 1)), 2 ** n)


def sign_agree(x, ref):
    return x is not None and ref is not None and x != 0 and np.sign(x) == np.sign(ref)


def h2c(luad_rows, e2_rows):
    luad = {r["gene"]: r for r in luad_rows if r["panel"] == "B"}
    e2 = {r["gene"]: r for r in e2_rows if r["panel"] == "B"}
    ref = sorted(g for g, r in luad.items() if r["status"] in REFERENCE_STATUSES)
    per_gene, tested = [], []
    for g in ref:
        r = e2.get(g, {})
        is_tested = r.get("del_median") is not None
        row = {"gene": g, "symbol": luad[g]["symbol"], "luad_status": luad[g]["status"],
               "luad_del_median": luad[g]["del_median"], "luad_ovx_median": luad[g]["ovx_median"],
               "e2_status": r.get("status", "ABSENT"), "e2_tested": is_tested,
               "e2_del_median": r.get("del_median"), "e2_ovx_median": r.get("ovx_median")}
        if is_tested:
            row["del_agree"] = bool(sign_agree(r["del_median"], luad[g]["del_median"]))
            row["ovx_agree"] = bool(sign_agree(r.get("ovx_median"), luad[g]["ovx_median"]))
            tested.append(row)
        per_gene.append(row)
    n = len(tested)
    a = sum(r["del_agree"] for r in tested)
    a_ovx = sum(r["ovx_agree"] for r in tested)
    if n < MIN_N:
        reading, p = "not_testable", None
    else:
        p = binom_upper(a, n)
        reading = "pattern_holds" if p <= ALPHA else "pattern_not_replicated"
    return {"reference_genes": len(ref), "n_tested": n, "del_agree": a, "ovx_agree_descriptive": a_ovx,
            "p_one_sided": None if p is None else float(p), "p_exact": None if p is None else str(p),
            "reading": reading, "per_gene": per_gene}


def status_table(luad_rows, e2_rows):
    luad = {r["gene"]: r for r in luad_rows if r["panel"] == "B"}
    e2 = {r["gene"]: r for r in e2_rows if r["panel"] == "B"}
    return [{"gene": g, "symbol": luad[g]["symbol"], "luad_status": luad[g]["status"],
             "e2_status": e2.get(g, {}).get("status", "ABSENT"),
             "same_status": luad[g]["status"] == e2.get(g, {}).get("status")} for g in sorted(luad)]


def median_rho(luad_rows, e2_rows):
    from scipy.stats import spearmanr
    luad = {r["gene"]: r.get("del_median") for r in luad_rows if r["panel"] == "B"}
    e2 = {r["gene"]: r.get("del_median") for r in e2_rows if r["panel"] == "B"}
    both = sorted(g for g in luad if luad[g] is not None and e2.get(g) is not None)
    if len(both) < 3:
        return {"n": len(both), "rho": None}
    return {"n": len(both), "rho": float(spearmanr([luad[g] for g in both], [e2[g] for g in both]).correlation)}


def mmr_split(per_donor, donor_mmr, genes):
    out = {}
    for g in genes:
        v = per_donor.get(f"{g}|delete") or {}
        for grp in ("MMRd", "MMRp"):
            xs = [x for d, x in v.items() if donor_mmr.get(d) == grp]
            out.setdefault(g, {})[grp] = {"n": len(xs), "median": float(np.median(xs)) if xs else None}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--luad-rows", required=True)
    ap.add_argument("--e2-rows", required=True)
    ap.add_argument("--e2-per-donor", required=True, help="analyse.py per_donor_adjusted.json")
    ap.add_argument("--cohort-definition", required=True, help="for donor MMR status")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    luad_rows, e2_rows = json.load(open(a.luad_rows)), json.load(open(a.e2_rows))
    res = h2c(luad_rows, e2_rows)
    tested_genes = [r["gene"] for r in e2_rows if r["panel"] == "B" and r.get("del_median") is not None]
    res["descriptive"] = {
        "status_table": status_table(luad_rows, e2_rows),
        "del_median_spearman_luad_vs_e2": median_rho(luad_rows, e2_rows),
        "mmr_split_del_median": mmr_split(json.load(open(a.e2_per_donor)),
                                          json.load(open(a.cohort_definition))["donor_mmr"], tested_genes),
        "note": "descriptive only; no reading changes on these",
    }
    json.dump(res, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in res.items() if k not in ("per_gene", "descriptive")}, indent=1))


if __name__ == "__main__":
    main()
