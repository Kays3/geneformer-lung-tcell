"""Amendment 4/5 primary analysis: null distribution of delete/overexpress concordance over eligible
random genes detectable in LUAD tumor T cells (never "genome-wide" -- ISP-STD-1 B.9, Stanley C1).

CPU only. Reads the completed phase8 markers (same shape as Phase 6's), computes each gene's raw
donor-level median shift per operation (no stratum control-adjustment -- Amendment 4 s.4.4), then the
one-sided permutation test on Spearman rho (H1: rho < 0), leave-one-out, and a gene-level bootstrap.
Never reads Shift_to_goal_end; reads the raw per-cell pickles via analyse.py, unchanged.

Fixes applied after Stanley's 2026-09-29 re-gate (B2, B3, C1; P2-P5/R1/R3-R5 from the first gate):
  B2: --frozen (frozen_100_estimable.json) is required. Every run gene's marker n_token_cells is
      compared against the frozen pre-GPU npos_by_donor for that (gene, donor). Any mismatch HALTS
      the analysis (exit nonzero, a diagnostic file is written, no primary status is emitted) -- this
      is a pipeline-integrity failure, not one of the ISP-STD-1 A.6 outcome statuses.
  B3: --noop-results (noop_spotchecks_<host>/results.jsonl) is required. Any failed spot check, or a
      missing spot check for a gene that should have one (every 20th gene, s.5.5), forces the primary
      status to `no_op_failed`, overriding whatever the statistics alone would have said.
  C1: the validity check runs on the 318 unique control genes in design.json (union across 19 strata),
      not "120" -- corrected everywhere. If it does not recover a negative rho, the top-level result
      carries `validity_check_failed: true` and a `WARNING` string; the primary status is still
      reported, but the human-readable summary says the result is not trusted.
  P2/R2: a gene counts as RUN only if all 86 markers (43 donors x 2 ops) exist; anything else is
         stopped_not_analysed, never fed to load_all (which would raise IncompleteRun uncaught).
  P3:    every leave-one-out recomputes the full permutation p, not just the sign of rho.
  P4:    opposite_direction uses its own upper-tail permutation p, not the primary's lower-tail p.
  P5/R1: status vocabulary per ISP-STD-1 A.6: `positive` = registered criterion met; `negative` = not
         met, reported with direction/p/n, never "no effect"; `opposite_direction`,
         `control_draw_sensitive_open`, `not_estimable`, `no_op_failed`, `stopped_not_analysed`.
  R5:    bootstrap and leave-one-out both use N_PERM_SECONDARY=2,000 permutations (min p=1/2,001);
         the primary test alone uses N_PERM=100,000.
"""
import argparse
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
N_PERM_SECONDARY = 2_000
N_BOOT = 10_000
ALPHA = 0.05
SEED = 20260929
NOOP_EVERY = 20


def gene_marker_ntoken(phase8_root, op, gene, donor):
    safe = donor.replace("/", "_")
    p = os.path.join(phase8_root, op, gene, f"{safe}.complete.json")
    if not os.path.exists(p):
        return None
    return json.load(open(p)).get("n_token_cells")


def gene_all_markers_present(phase8_root, gene, donors):
    for op in OPS:
        for d in donors:
            if gene_marker_ntoken(phase8_root, op, gene, d) is None:
                return False
    return True


def check_gate_identity(phase8_root, gene, donors, frozen_npos):
    """B2: post-GPU marker n_token_cells must equal the pre-GPU frozen count for every donor, both ops
    (n_token_cells does not depend on delete/overexpress, but both markers record it independently, so
    both are checked). Returns a list of mismatch dicts, empty if none."""
    mismatches = []
    for op in OPS:
        for d in donors:
            marker_n = gene_marker_ntoken(phase8_root, op, gene, d)
            frozen_n = frozen_npos.get(d)
            if marker_n is not None and frozen_n is not None and marker_n != frozen_n:
                mismatches.append({"gene": gene, "op": op, "donor": d, "marker_n_token_cells": marker_n,
                                    "frozen_npos": frozen_n})
    return mismatches


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


def collect_gene_rows(phase8_root, gene_list, donors, positions, n_total, frozen_npos_by_gene=None):
    """frozen_npos_by_gene: {gene: {donor: n_pos}} for the gate-identity check (B2); None skips it
    (used for the R3 validity check on the existing controls, which has no pre-GPU freeze of its own --
    that data was already gated by Amendment 3h's own count checks at the time it was produced)."""
    rows = []
    all_mismatches = []
    for g in gene_list:
        if not gene_all_markers_present(phase8_root, g, donors):
            rows.append({"gene": g, "status": "stopped_not_analysed"})
            continue
        if frozen_npos_by_gene is not None:
            mism = check_gate_identity(phase8_root, g, donors, frozen_npos_by_gene.get(g, {}))
            if mism:
                all_mismatches.extend(mism)
                continue  # don't compute medians for a gene under investigation
        m = gene_medians(phase8_root, g, donors, positions, n_total, CONTROL_MIN_CELLS)
        if m["delete"] is None or m["overexpress"] is None:
            rows.append({"gene": g, "status": "not_estimable", "detail": m})
            continue
        rows.append({"gene": g, "status": "estimable", "delete_median": m["delete"]["median"],
                     "overexpress_median": m["overexpress"]["median"],
                     "delete_n_donors": m["delete"]["n_donors"], "overexpress_n_donors": m["overexpress"]["n_donors"]})
    return rows, all_mismatches


def primary_from_rows(rows, rng):
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
    elif rho > 0:
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


def apply_noop_gate(primary, noop_results_path, expected_genes):
    """B3: consume noop_spotcheck.py's results. Any fail, or a missing check for a gene that was due
    one (every NOOP_EVERY-th gene in run order), forces no_op_failed regardless of the statistics."""
    checked = {}
    if os.path.exists(noop_results_path):
        for line in open(noop_results_path):
            line = line.strip()
            if not line:
                continue
            r = json.loads(line)
            checked[r["gene"]] = r.get("pass", False)
    due = {expected_genes[i] for i in range(len(expected_genes)) if (i + 1) % NOOP_EVERY == 0}
    missing = sorted(due - set(checked))
    failed = sorted(g for g, ok in checked.items() if not ok)
    noop_report = {"due": sorted(due), "checked": checked, "missing": missing, "failed": failed}
    if missing or failed:
        primary = dict(primary, status="no_op_failed", no_op_detail=noop_report)
    else:
        primary = dict(primary, no_op_detail=noop_report)
    return primary


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase8-root", required=True, help="run_phase8_null.sh's --out tree, both ops")
    ap.add_argument("--null-genes", required=True, help="null_genes_100_estimable.json (Amendment 5 B1)")
    ap.add_argument("--frozen", required=True, help="frozen_100_estimable.json, for the gate-identity check (B2)")
    ap.add_argument("--noop-results", required=True, help="noop_spotchecks_<host>/results.jsonl (B3)")
    ap.add_argument("--design", required=True, help="phase7_results/design.json")
    ap.add_argument("--ovx-index", required=True, help="the null-genes ovx index (built post-GPU, Amendment 5 P1)")
    ap.add_argument("--phase6-controls-root", default=None,
                     help="R3 validity check: the CLOSED analysis's phase6_all tree, for the 318 existing "
                          "control genes; if omitted, the validity check is skipped and noted as such")
    ap.add_argument("--controls-ovx-index", default=None, help="R3: the closed analysis's ovx_index.json")
    ap.add_argument("--ambient", default=None, help="ambient_loao_null.py output, for R4's per-gene flag")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    ng = json.load(open(a.null_genes))
    assert ng["n"] == 100 and len(ng["draw_order"]) == 100, "Amendment 5 registers exactly N=100 estimable genes"
    design = an.Design(**json.load(open(a.design)))
    ix = json.load(open(a.ovx_index))
    positions = {tuple(k.split("|")): v for k, v in ix["positions"].items()}

    frozen = json.load(open(a.frozen))
    frozen_npos_by_gene = {r["gene"]: r["npos_by_donor"] for r in frozen["rows"] if r["estimable"]}

    rng = np.random.default_rng(SEED)

    rows, mismatches = collect_gene_rows(a.phase8_root, ng["draw_order"], design.donors, positions,
                                          ix["n_total"], frozen_npos_by_gene=frozen_npos_by_gene)
    if mismatches:
        # B2: a gate-identity mismatch HALTS the analysis -- this is not an ISP-STD-1 A.6 outcome
        # status, it is a pipeline-integrity failure that needs investigation before any status is
        # reported at all.
        json.dump({"about": "HALTED: gate-identity mismatch (Amendment 5 s.5.2/B2)",
                   "n_mismatches": len(mismatches), "mismatches": mismatches}, open(a.out, "w"), indent=1)
        print(f"HALTED: {len(mismatches)} gate-identity mismatches, see {a.out}", file=sys.stderr)
        return 2

    primary, x, y = primary_from_rows(rows, rng)
    if primary.get("status") != "not_estimable":
        noop_path = a.noop_results
        primary = apply_noop_gate(primary, noop_path, ng["draw_order"])

    result = {"about": __doc__.strip(), "n_registered": ng["n"], "rows": rows, "primary": primary}

    if a.phase6_controls_root and a.controls_ovx_index:
        cix = json.load(open(a.controls_ovx_index))
        cpositions = {tuple(k.split("|")): v for k, v in cix["positions"].items()}
        control_ids = sorted({c for st in design.strata.values() for c in st.get("controls", [])})
        crows, _ = collect_gene_rows(a.phase6_controls_root, control_ids, design.donors, cpositions,
                                      cix["n_total"], frozen_npos_by_gene=None)
        cprimary, _, _ = primary_from_rows(crows, rng)
        validity_failed = not (cprimary.get("rho") is not None and cprimary["rho"] < 0)
        result["validity_check_318_control_genes"] = {
            "about": "same raw-shift/permutation pipeline applied to the 318 unique control genes already "
                     "computed in design.json (union across 19 strata); sign should be compatible with "
                     "RESULTS.md s.4's control-adjusted rho=-0.62 over 360 control x stratum entries (not the "
                     "same quantity -- s.4.4/R3/Stanley C1 -- so magnitudes are not compared, sign only).",
            "rho": cprimary.get("rho"), "n_estimable": cprimary.get("n_estimable"), "status": cprimary.get("status"),
            "sign_compatible_with_minus_0_62": not validity_failed,
        }
        result["validity_check_failed"] = validity_failed
        if validity_failed:
            result["WARNING"] = ("validity check did not recover a negative rho on the 318 control genes -- "
                                  "the primary result below is reported but NOT TRUSTED until this is resolved.")
    else:
        result["validity_check_318_control_genes"] = {"status": "SKIPPED", "reason": "no --phase6-controls-root/--controls-ovx-index given"}
        result["validity_check_failed"] = None

    if a.ambient:
        amb = {r["ensembl_id"]: r for r in json.load(open(a.ambient))["rows"]}
        for r in rows:
            r["ambient"] = amb.get(r["gene"], {"ambient_label": "NOT_SCORED"})
    else:
        result["ambient_flag"] = {"status": "SKIPPED", "reason": "no --ambient given"}

    if x is not None:
        result["null_delete_medians"] = x.tolist()
        result["null_overexpress_medians"] = y.tolist()
        result["null_gene_order"] = [r["gene"] for r in rows if r["status"] == "estimable"]

    json.dump(result, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("rows", "null_delete_medians", "null_overexpress_medians", "null_gene_order")},
                     indent=1))
    return 0 if primary.get("status") != "no_op_failed" else 1


if __name__ == "__main__":
    sys.exit(main())
