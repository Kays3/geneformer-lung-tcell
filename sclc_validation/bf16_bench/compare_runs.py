#!/usr/bin/env python3
"""Our port of the colleague's ISP accuracy gate (analysis/14_compare_isp.py
in petadimensionlab/Geneformer, gate "G5": a lower-precision run is only
acceptable if it reproduces the *ranking* of the shift score -- Spearman
rho, top-N overlap, sign agreement -- not just a scalar loss), adapted to
our two workloads and table schemas instead of theirs.

Two subcommands:

  panel   -- per-gene targeted-panel comparison. Reads the per-comparison
             InSilicoPerturberStats CSVs directly from each run's
             <stats-root>/{delete,overexpress}/targeted_*_to_*.csv (the same
             files run_stats() in run_targeted_panel.py writes; we do not
             require the separate analyze_targeted_results.py merge step).
             Key = (comparison, Gene_name); value = Shift_to_goal_end.

  t4      -- set-level T4 program-phase comparison. Reads
             t4_shift_summary.csv (scope == "all" rows only) and
             t4_matched_null_test.csv from each run.

Every number in RESULTS_BF16.md must come from the JSON this script writes
-- no hand-typed numbers (see task boundary "Verification").
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

# Pre-registered acceptance gate (frozen before the bf16 arm runs; see
# hive/handoffs/bf16-isp-replication-plan-20260917.md "Pre-registered
# acceptance gate"). Changing these after the bf16 arm has run is a dated
# amendment, not a silent edit.
RHO_MIN = 0.95
TOPN = 20
TOPN_MIN_OVERLAP = 19


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) < 3:
        return float("nan")
    try:
        from scipy.stats import spearmanr

        return float(spearmanr(a, b).statistic)
    except Exception:
        ra = pd.Series(a).rank().to_numpy()
        rb = pd.Series(b).rank().to_numpy()
        return float(np.corrcoef(ra, rb)[0, 1])


def topn_overlap(keys: pd.Series, values: np.ndarray, other_values: np.ndarray, n: int) -> int:
    order_a = np.argsort(-np.abs(values))[:n]
    order_b = np.argsort(-np.abs(other_values))[:n]
    top_a = set(keys.iloc[order_a])
    top_b = set(keys.iloc[order_b])
    return len(top_a & top_b)


# ---------------------------------------------------------------- panel ----

def load_panel_stats(stats_root: Path, perturb_type: str) -> pd.DataFrame:
    rows = []
    for f in sorted((stats_root / perturb_type).glob("targeted_*_to_*.csv")):
        # filename: targeted_{delete,overexpress}_{source}_to_{target}.csv
        comparison = f.stem.split(f"targeted_{perturb_type}_", 1)[1]
        df = pd.read_csv(f)
        df["comparison"] = comparison
        rows.append(df)
    if not rows:
        raise SystemExit(f"no stats CSVs found under {stats_root / perturb_type}")
    return pd.concat(rows, ignore_index=True)


def compare_panel_arm(df_a: pd.DataFrame, df_b: pd.DataFrame) -> dict:
    key = ["comparison", "Gene_name"]
    m = df_a.merge(df_b, on=key, suffixes=("_a", "_b"))
    if m.empty:
        return {"error": "no shared (comparison, Gene_name) rows"}
    x = m["Shift_to_goal_end_a"].to_numpy(dtype=float)
    y = m["Shift_to_goal_end_b"].to_numpy(dtype=float)
    diff = np.abs(x - y)
    keys = m["comparison"] + "/" + m["Gene_name"]
    return {
        "n_rows": int(len(m)),
        "n_genes": int(m["Gene_name"].nunique()),
        "spearman": round(spearman(x, y), 5),
        "sign_agreement_all": round(float(np.mean(np.sign(x) == np.sign(y))), 4),
        "topn": TOPN,
        "topn_overlap": topn_overlap(keys, x, y, TOPN),
        "mean_abs_diff": round(float(diff.mean()), 6),
        "median_abs_diff": round(float(np.median(diff)), 6),
        "max_abs_diff": round(float(diff.max()), 6),
        "mean_abs_shift_a": round(float(np.abs(x).mean()), 6),
        "_diff": diff,       # dropped before json.dump; used for signal split
        "_abs_shift_a": np.abs(x),
        "_sign_match": np.sign(x) == np.sign(y),
    }


def panel_command(args: argparse.Namespace) -> dict:
    out = {"a_stats_root": str(args.a_stats), "b_stats_root": str(args.b_stats),
           "floor_stats_root": str(args.floor_stats), "gate": {"rho_min": RHO_MIN,
           "topn": TOPN, "topn_min_overlap": TOPN_MIN_OVERLAP}, "by_perturb_type": {}}
    overall_pass = True
    for ptype in args.perturb_types:
        df_a = load_panel_stats(args.a_stats, ptype)
        df_b = load_panel_stats(args.b_stats, ptype)
        df_floor = load_panel_stats(args.floor_stats, ptype)

        main = compare_panel_arm(df_a, df_b)
        floor = compare_panel_arm(df_a, df_floor)
        if "error" in main or "error" in floor:
            out["by_perturb_type"][ptype] = {"error": main.get("error") or floor.get("error")}
            overall_pass = False
            continue

        # "genes whose |shift| exceeds the fp32-vs-fp32 noise floor": use the
        # floor run's max |diff| as the noise threshold on the fp32-baseline
        # |shift| scale -- below it, a gene's fp32-vs-fp32 rerun already
        # can't be trusted to agree in sign, so bf16 isn't held to a
        # different standard there either.
        floor_threshold = floor["max_abs_diff"]
        signal_mask = main["_abs_shift_a"] > floor_threshold
        n_signal = int(signal_mask.sum())
        sign_agreement_signal = (
            round(float(main["_sign_match"][signal_mask].mean()), 4) if n_signal else None
        )
        bf16_excess_signal = (
            round(float(np.median(main["_diff"][signal_mask])), 6) if n_signal else None
        )
        bf16_excess_nearzero = (
            round(float(np.median(main["_diff"][~signal_mask])), 6) if (~signal_mask).any() else None
        )

        rho_pass = not np.isnan(main["spearman"]) and main["spearman"] >= RHO_MIN
        topn_pass = main["topn_overlap"] >= TOPN_MIN_OVERLAP
        sign_pass = n_signal == 0 or sign_agreement_signal == 1.0
        ptype_pass = rho_pass and topn_pass and sign_pass
        overall_pass = overall_pass and ptype_pass

        for d in (main, floor):
            d.pop("_diff", None); d.pop("_abs_shift_a", None); d.pop("_sign_match", None)

        out["by_perturb_type"][ptype] = {
            "bf16_vs_fp32": main,
            "fp32_vs_fp32_floor": floor,
            "n_signal_genes": n_signal,
            "n_genes_total": int(len(signal_mask)),
            "sign_agreement_signal_genes": sign_agreement_signal,
            "bf16_median_abs_diff_signal_genes": bf16_excess_signal,
            "bf16_median_abs_diff_near_zero_genes": bf16_excess_nearzero,
            "gate_checks": {"rho_pass": rho_pass, "topn_pass": topn_pass, "sign_pass": sign_pass},
            "pass": ptype_pass,
        }
    out["pass"] = overall_pass
    return out


# -------------------------------------------------------------------- t4 ---

def t4_command(args: argparse.Namespace) -> dict:
    shift_a = pd.read_csv(args.a_shift)
    shift_b = pd.read_csv(args.b_shift)
    shift_a = shift_a[shift_a["scope"] == "all"]
    shift_b = shift_b[shift_b["scope"] == "all"]
    key = ["item_id", "phase", "program", "nested_size", "null_iteration", "source_state", "target_state"]
    m = shift_a.merge(shift_b, on=key, suffixes=("_a", "_b"), how="inner")
    shift_result: dict = {"n_rows": int(len(m))}
    if m.empty:
        shift_result["error"] = "no shared rows between the two shift summaries (scope=='all')"
    else:
        x = m["mean_shift_a"].to_numpy(dtype=float)
        y = m["mean_shift_b"].to_numpy(dtype=float)
        shift_result.update({
            "spearman": round(spearman(x, y), 5),
            "pearson": round(float(np.corrcoef(x, y)[0, 1]), 5),
            "mean_abs_diff": round(float(np.abs(x - y).mean()), 6),
            "max_abs_diff": round(float(np.abs(x - y).max()), 6),
        })

    null_a = pd.read_csv(args.a_null)
    null_b = pd.read_csv(args.b_null)
    nkey = ["program", "source_state", "target_state"]
    nm = null_a.merge(null_b, on=nkey, suffixes=("_a", "_b"), how="inner")
    nm["sig_a"] = nm["empirical_p_directional_a"] < args.alpha
    nm["sig_b"] = nm["empirical_p_directional_b"] < args.alpha
    flips = nm[nm["sig_a"] != nm["sig_b"]]
    null_result = {
        "n_rows": int(len(nm)),
        "alpha": args.alpha,
        "n_flips": int(len(flips)),
        "flipped": flips[nkey + ["sig_a", "sig_b", "empirical_p_directional_a", "empirical_p_directional_b"]]
                   .to_dict(orient="records"),
    }

    out = {
        "a_shift": str(args.a_shift), "b_shift": str(args.b_shift),
        "a_null": str(args.a_null), "b_null": str(args.b_null),
        "shift_correlation": shift_result,
        "matched_null_qualitative_agreement": null_result,
        "pass": null_result["n_flips"] == 0 and "error" not in shift_result,
    }
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)

    p_panel = sub.add_parser("panel", help="per-gene targeted-panel gate")
    p_panel.add_argument("--a-stats", type=Path, required=True, help="fp32 baseline stats root")
    p_panel.add_argument("--b-stats", type=Path, required=True, help="bf16 (candidate) stats root")
    p_panel.add_argument("--floor-stats", type=Path, required=True, help="fp32 repeat (noise floor) stats root")
    p_panel.add_argument("--perturb-types", nargs="+", choices=("delete", "overexpress"),
                          default=["delete", "overexpress"],
                          help="Which perturb types were run for this arm set (default: both). "
                               "Pass '--perturb-types overexpress' for an overexpress-only sizing amendment.")
    p_panel.add_argument("--out", type=Path, required=True)

    p_t4 = sub.add_parser("t4", help="set-level T4 program-phase gate")
    p_t4.add_argument("--a-shift", type=Path, required=True)
    p_t4.add_argument("--b-shift", type=Path, required=True)
    p_t4.add_argument("--a-null", type=Path, required=True)
    p_t4.add_argument("--b-null", type=Path, required=True)
    p_t4.add_argument("--alpha", type=float, default=0.05)
    p_t4.add_argument("--out", type=Path, required=True)

    args = ap.parse_args()
    result = panel_command(args) if args.mode == "panel" else t4_command(args)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(result, indent=2, default=str) + "\n")
    print(json.dumps(result, indent=2, default=str))
    print(f"\nverdict -> {'PASS' if result['pass'] else 'FAIL'}  ({args.out})")


if __name__ == "__main__":
    main()
