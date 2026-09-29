"""Amendment 5 C2: secondary percentile join (s.4.3/R4), committed and hashed before any phase8 output
is read, per Stanley. Descriptive only -- does not gate any status (ISP-STD-1 B.9).

For every closed-analysis panel gene with status T_CELL_SIGNAL_TOWARD or T_CELL_SIGNAL_AWAY
(outcome_rows.json), recomputes its RAW (non-control-adjusted) donor-level median shift -- the same
raw quantity the null study uses (Amendment 4 s.4.4), NOT the control-adjusted value RESULTS.md
reports, so the comparison to the null distribution is apples-to-apples (raw vs raw). Then reports
that gene's percentile rank within the null's raw delete/overexpress median distributions
(null_analysis.py's output). CPU only; reads the CLOSED analysis's already-computed phase6_all tree,
no new GPU.
"""
import argparse
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import analyse as an  # noqa: E402

CONTROL_MIN_CELLS = 10
D_MIN = 10


def raw_median(phase6_root, gene, donors, positions, n_total, op):
    calls = an.load_all(phase6_root, [gene], donors, positions, n_total)
    dm = an.donor_means(calls, gene, op, an.GOAL, CONTROL_MIN_CELLS)
    estimable = {d: s for d, (s, n) in dm.items() if n >= CONTROL_MIN_CELLS}
    if len(estimable) < D_MIN:
        return None
    return float(np.median(list(estimable.values())))


def percentile_rank(value, distribution):
    """Fraction of the null distribution at or below value (descriptive only)."""
    arr = np.asarray(distribution)
    return float((arr <= value).mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--outcome-rows", required=True, help="the CLOSED analysis's phase7_results/outcome_rows.json")
    ap.add_argument("--phase6-panel-root", required=True, help="the CLOSED analysis's phase6_all tree")
    ap.add_argument("--panel-ovx-index", required=True, help="the CLOSED analysis's ovx_index.json")
    ap.add_argument("--design", required=True)
    ap.add_argument("--null-result", required=True, help="null_analysis.py's output (for null_delete_medians "
                                                           "and null_overexpress_medians)")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    outcome_rows = json.load(open(a.outcome_rows))
    rows = outcome_rows if isinstance(outcome_rows, list) else outcome_rows["rows"]
    target_genes = [r for r in rows if r.get("status") in ("T_CELL_SIGNAL_TOWARD", "T_CELL_SIGNAL_AWAY")]

    design = an.Design(**json.load(open(a.design)))
    ix = json.load(open(a.panel_ovx_index))
    positions = {tuple(k.split("|")): v for k, v in ix["positions"].items()}

    null_result = json.load(open(a.null_result))
    null_delete = null_result.get("null_delete_medians")
    null_overexpress = null_result.get("null_overexpress_medians")
    if null_delete is None or null_overexpress is None:
        json.dump({"about": __doc__.strip(), "status": "SKIPPED",
                   "reason": "null_analysis.py's primary did not reach a computable null distribution "
                             "(n_estimable < 10); no distribution to place panel genes against"},
                  open(a.out, "w"), indent=1)
        print("SKIPPED: no null distribution available")
        return 0

    out_rows = []
    for r in target_genes:
        gene = r["gene"]
        d_med = raw_median(a.phase6_panel_root, gene, design.donors, positions, ix["n_total"], "delete")
        o_med = raw_median(a.phase6_panel_root, gene, design.donors, positions, ix["n_total"], "overexpress")
        out_rows.append({
            "gene": gene, "symbol": r.get("symbol"), "closed_analysis_status": r["status"],
            "raw_delete_median": d_med, "raw_overexpress_median": o_med,
            "delete_percentile_in_null": percentile_rank(d_med, null_delete) if d_med is not None else None,
            "overexpress_percentile_in_null": percentile_rank(o_med, null_overexpress) if o_med is not None else None,
        })

    json.dump({"about": __doc__.strip(), "n_null_genes": len(null_delete), "rows": out_rows},
               open(a.out, "w"), indent=1)
    print(json.dumps(out_rows, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
