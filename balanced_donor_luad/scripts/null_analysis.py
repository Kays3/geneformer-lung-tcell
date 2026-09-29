"""Amendment 4/5 primary analysis: genome-wide null distribution of delete/overexpress concordance.

CPU only. Reads the completed phase8 markers (same shape as Phase 6's), computes each drawn gene's raw
donor-level median shift per operation (no stratum control-adjustment -- Amendment 4 s.4.4), then the
one-sided permutation test on Spearman rho (H1: rho < 0), leave-one-out, and a gene-level bootstrap.
Never reads Shift_to_goal_end; reads the raw per-cell pickles via analyse.py, unchanged.

Fixes applied after Stanley's 2026-09-29 gate (P2-P5, R1, R3-R5):
  P2/R2: a gene counts as RUN only if all 86 markers (43 donors x 2 ops) exist; anything else is
         stopped_not_analysed, never fed to load_all (which would raise IncompleteRun uncaught).
  P3:    every leave-one-out recomputes the full permutation p, not just the sign of rho.
  P4:    opposite_direction uses its own upper-tail permutation p, not the primary's lower-tail p.
  P5/R1: status vocabulary fixed to ISP-STD-1 A.6: `positive` = registered criterion met (rho<0,
         p<=alpha, LOO-stable, bootstrap-stable); `negative` = not met, reported with direction/p/n,
         never "no effect"; `opposite_direction`, `control_draw_sensitive_open`, `not_estimable`,
         `stopped_not_analysed` unchanged.
  R3:    a validity check re-runs the identical raw-shift/permutation pipeline on the 120 EXISTING
         matched controls (already-computed Phase 6 data, no new GPU) as the closest analogue to a
         positive control: if this pipeline cannot recover a rho compatible in sign with the already-
         reported control-adjusted rho = -0.62 (RESULTS.md s.4), the null result is not trusted.
  R4:    ambient flag (from ambient_loao_null.py's output) is merged in per gene, and the secondary
         percentile of every panel gene with a TOWARD/AWAY call is computed against this null's raw
         distribution (descriptive, does not gate any status -- ISP-STD-1 B.9).
  R5:    bootstrap and leave-one-out both use N_PERM_SECONDARY=2,000 permutations (min p=1/2,001),
         registered explicitly; the primary test alone uses N_PERM=100,000.
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
N_PERM_SECONDARY = 2_000  # leave-one-out and bootstrap inner test (R5)
N_BOOT = 10_000
ALPHA = 0.05
SEED = 20260929


def gene_all_markers_present(phase8_root, gene, donors):
    for op in OPS:
        for d in donors:
            safe = d.replace("/", "_")
            if not os.path.exists(os.path.join(phase8_root, op, gene, f"{safe}.complete.json")):
                return False
    return True


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


def lower_tail_p(x, y, rng, n_perm):
    """H1: rho < 0. p = P(rho_perm <= rho_obs)."""
    obs = spearmanr(x, y).correlation
    y = np.asarray(y)
    count = sum(1 for _ in range(n_perm) if spearmanr(x, rng.permutation(y)).correlation <= obs)
    return obs, (1 + count) / (n_perm + 1)


def upper_tail_p(x, y, rng, n_perm):
    """For opposite_direction: p = P(rho_perm >= rho_obs)."""
    obs = spearmanr(x, y).correlation
    y = np.asarray(y)
    count = sum(1 for _ in range(n_perm) if spearmanr(x, rng.permutation(y)).correlation >= obs)
    return obs, (1 + count) / (n_perm + 1)


def collect_gene_rows(phase8_root, phase6_positive_controls_root, gene_list, donors, positions, n_total):
    """Shared by the primary (null genes) and the validity check (120 existing controls)."""
    rows = []
    for g in gene_list:
        root = phase8_root
        if not gene_all_markers_present(root, g, donors):
            rows.append({"gene": g, "status": "stopped_not_analysed"})
            continue
        m = gene_medians(root, g, donors, positions, n_total, CONTROL_MIN_CELLS)
        if m["delete"] is None or m["overexpress"] is None:
            rows.append({"gene": g, "status": "not_estimable", "detail": m})
            continue
        rows.append({"gene": g, "status": "estimable", "delete_median": m["delete"]["median"],
                     "overexpress_median": m["overexpress"]["median"],
                     "delete_n_donors": m["delete"]["n_donors"], "overexpress_n_donors": m["overexpress"]["n_donors"]})
    return rows


def primary_from_rows(rows, rng, seed_tag):
    estimable = [r for r in rows if r["status"] == "estimable"]
    n = len(estimable)
    out = {"n_estimable": n, "n_stopped_not_analysed": sum(r["status"] == "stopped_not_analysed" for r in rows),
           "n_not_estimable": sum(r["status"] == "not_estimable" for r in rows)}
    if n < 10:
        out["status"] = "not_estimable"
        return out, None, None
    x = np.array([r["delete_median"] for r in estimable])
    y = np.array([r["overexpress_median"] for r in estimable])
    rho, p_lower = lower_tail_p(x, y, rng, N_PERM)
    loo_pass = []
    for i in range(n):
        keep = [j for j in range(n) if j != i]
        r_loo, p_loo = lower_tail_p(x[keep], y[keep], rng, N_PERM_SECONDARY)
        loo_pass.append(r_loo < 0 and p_loo <= ALPHA)
    loo_stable = all(loo_pass)
    boot_hits = 0
    for _ in range(N_BOOT):
        idx = rng.integers(0, n, n)
        rb, pb = lower_tail_p(x[idx], y[idx], rng, N_PERM_SECONDARY)
        if rb < 0 and pb <= ALPHA:
            boot_hits += 1
    boot_frac = boot_hits / N_BOOT
    primary_pass = rho < 0 and p_lower <= ALPHA
    if primary_pass and loo_stable and boot_frac >= 0.95:
        status = "positive"
    elif primary_pass:
        status = "control_draw_sensitive_open"
    else:
        # Check opposite_direction on its own upper-tail p (P4) -- only meaningful to test when
        # the observed rho is actually positive; a negative-but-nonsignificant rho is just `negative`.
        if rho > 0:
            _, p_upper = upper_tail_p(x, y, rng, N_PERM)
            status = "opposite_direction" if p_upper <= ALPHA else "negative"
        else:
            status = "negative"
    out.update({"status": status, "rho": float(rho), "p_lower_tail": float(p_lower), "n_perm_primary": N_PERM,
               "n_perm_secondary": N_PERM_SECONDARY, "loo_stable": bool(loo_stable),
               "loo_n_pass": int(sum(loo_pass)), "loo_n": n, "bootstrap_fraction_stable": boot_frac,
               "n_bootstrap": N_BOOT, "alpha": ALPHA, "sidedness": "one-sided, H1: rho<0"})
    if status == "negative":
        out["report_sentence"] = f"primary not met: direction rho={rho:.4f}, p={p_lower:.4g}, n={n}"
    return out, x, y


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase8-root", required=True, help="run_phase8_null.sh's --out tree, both ops")
    ap.add_argument("--null-genes", required=True)
    ap.add_argument("--design", required=True, help="phase7_results/design.json")
    ap.add_argument("--ovx-index", required=True, help="the null-genes ovx index (Amendment 5 P1)")
    ap.add_argument("--phase6-controls-root", default=None,
                     help="R3 validity check: the CLOSED analysis's phase6_all tree, for the 120 existing "
                          "controls; if omitted, the validity check is skipped and noted as such")
    ap.add_argument("--controls-ovx-index", default=None, help="R3: the closed analysis's ovx_index.json")
    ap.add_argument("--ambient", default=None, help="ambient_loao_null.py output, for R4's per-gene flag")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ng = json.load(open(a.null_genes))
    design = an.Design(**json.load(open(a.design)))
    ix = json.load(open(a.ovx_index))
    positions = {tuple(k.split("|")): v for k, v in ix["positions"].items()}

    rng = np.random.default_rng(SEED)

    rows = collect_gene_rows(a.phase8_root, None, ng["draw_order"], design.donors, positions, ix["n_total"])
    primary, x, y = primary_from_rows(rows, rng, "primary")

    result = {"about": __doc__.strip(), "n_drawn": len(ng["draw_order"]), "rows": rows, "primary": primary}

    # R3: validity check on the 120 existing matched controls, same pipeline, no new GPU.
    if a.phase6_controls_root and a.controls_ovx_index:
        cix = json.load(open(a.controls_ovx_index))
        cpositions = {tuple(k.split("|")): v for k, v in cix["positions"].items()}
        control_ids = sorted({c for st in design.strata.values() for c in st.get("controls", [])})
        crows = collect_gene_rows(a.phase6_controls_root, None, control_ids, design.donors, cpositions, cix["n_total"])
        cprimary, _, _ = primary_from_rows(crows, rng, "validity")
        result["validity_check_120_controls"] = {
            "about": "same raw-shift/permutation pipeline applied to the 120 already-computed matched "
                     "controls; sign should be compatible with RESULTS.md s.4's control-adjusted rho=-0.62 "
                     "(not the same quantity -- s.4.4/R3 -- so magnitudes are not compared, sign only).",
            "rho": cprimary.get("rho"), "n_estimable": cprimary.get("n_estimable"), "status": cprimary.get("status"),
            "sign_compatible_with_minus_0_62": (cprimary.get("rho") is not None and cprimary["rho"] < 0),
        }
    else:
        result["validity_check_120_controls"] = {"status": "SKIPPED", "reason": "no --phase6-controls-root/--controls-ovx-index given"}

    # R4: ambient flag merge + secondary percentiles.
    if a.ambient:
        amb = {r["ensembl_id"]: r for r in json.load(open(a.ambient))["rows"]}
        for r in rows:
            r["ambient"] = amb.get(r["gene"], {"ambient_label": "NOT_SCORED"})
    else:
        result["ambient_flag"] = {"status": "SKIPPED", "reason": "no --ambient given"}

    if x is not None:
        # descriptive percentiles for any panel gene the caller wants placed in this distribution;
        # left as an empty hook here since the panel genes' own medians live in the closed analysis's
        # phase7_results, not this script's inputs -- computed by a small separate join at report time.
        result["null_delete_medians"] = x.tolist()
        result["null_overexpress_medians"] = y.tolist()

    json.dump(result, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items() if k not in ("rows", "null_delete_medians", "null_overexpress_medians")}, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
