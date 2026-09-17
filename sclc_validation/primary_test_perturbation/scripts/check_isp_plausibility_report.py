#!/usr/bin/env python3
"""Check every load-bearing number in isp_plausibility_20260917.md against the data.

The repository's rule, applied to prose: if a number is in a committed report, a
committed script must produce it. Expected values are DERIVED here from the source
tables and the generated candidate table -- never typed in -- and each check asserts
the report text CONTAINS the derived string.

Two deliberate design points:
  * The report is whitespace-collapsed before searching. It is hard-wrapped, so a
    phrase spanning a line break is invisible to a raw substring search and the
    check would pass or fail for the wrong reason.
  * A control check runs first. If it fails, the report could not be read and every
    "not found" below would be meaningless rather than informative.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
PTP = HERE.parent
SCLC = PTP.parent
REPORT = PTP / "reports" / "isp_plausibility_20260917.md"
CANDIDATES = PTP / "tables" / "isp_plausibility_top_candidates.csv"
SHIFT = PTP / "tables" / "allgene_delete_overexpress_shift.csv"
AMBIENT = PTP / "tables" / "ambient_risk_all_genes.csv"
MANIFEST = PTP / "tables" / "ambient_risk_manifest.json"
PAM = SCLC / "immune_axis_test" / "results" / "axis_sensitivity_summary.json"

failures: list[str] = []
checks = 0


def check(label: str, needle: str, haystack: str) -> None:
    global checks
    checks += 1
    if needle not in haystack:
        failures.append(f"{label}: report does not contain {needle!r}")


def main() -> int:
    global checks
    text = re.sub(r"\s+", " ", REPORT.read_text())

    # ---- control: prove the search can find something before trusting an absence
    control = "Biological plausibility of the complete ISP screen"
    if control not in text:
        print(f"CONTROL FAILED: {control!r} not found -- the report was not read; "
              "every check below would be meaningless.")
        return 2
    print(f"control ok: report read, {len(text)} chars collapsed")

    shift = pd.read_csv(SHIFT)
    cand = pd.read_csv(CANDIDATES)
    amb = pd.read_csv(AMBIENT)
    manifest = json.loads(MANIFEST.read_text())
    pam = json.loads(PAM.read_text())

    sel = cand[cand["selected_by"] == "top_n_by_effect"]
    focus = cand[cand["Gene_name"].isin(["S100A8", "S100A9"])]

    # ---- screen-level counts
    gate = (shift["delete_fdr"] < 0.05) & (shift["overexpress_fdr"] < 0.05) & (shift["delete_n"] >= 25)
    same_sign = (shift["delete_shift"] * shift["overexpress_shift"]) > 0
    check("screen rows", f"{len(shift):,}", text)
    check("gated rows", f"{int(gate.sum()):,}", text)
    check("gated same-sign", f"{int((gate & same_sign).sum()):,}", text)
    check("concordant rows", f"{int((shift['concordant'] == True).sum()):,}", text)
    check("selected rows", f"{len(cand)} (120 top-N", text)
    check("distinct genes in top-N", f"{sel['Gene_name'].nunique()}", text)
    check("median delete_n", f"{sel['delete_n'].median():.0f}", text)

    # ---- top-N composition
    check("ambient-flagged", f"**{int(sel['ambient_flag'].sum())} / {len(sel)}**", text)
    check("ambient anchors", f"**{int(sel['is_ambient_anchor'].sum())}**", text)
    check("tcell anchors zero", f"**{int(sel['is_tcell_anchor'].sum())}**", text)
    check("same-sign in top-N", f"**{int(sel['same_sign_delete_overexpress'].sum())} / {len(sel)}**", text)
    check("non-T lineage", f"**{int(sel['non_t_lineage'].sum())} / {len(sel)}**", text)
    check("pam excluded", f"| {int(sel['pam_excluded'].sum())} |", text)

    # ---- lineage counts
    for lin, n in sel["lineage_call"].value_counts().items():
        if lin in ("lung_epithelial", "myeloid", "neuroendocrine"):
            check(f"lineage {lin}", f"| {n} |", text)

    # ---- flag threshold and cohort size
    check("flag threshold", f"{manifest['flag_threshold']:.4f}"[:6], text)
    check("n_cells", f"{manifest['n_cells']:,}", text)
    check("anchors_tcell count", f"{len(manifest['anchors_tcell'])} canonical T-cell anchors", text)
    check("pam exclusion count", str(pam["n_excluded"]), text)

    # ---- anchor detection rates quoted in the prose
    a = amb.set_index("gene")
    for gene in ("CD3E", "IL7R", "CD8A"):
        check(f"{gene} detect_frac", f"{a.loc[gene, 'detect_frac']:.3f}", text)

    # ---- every S100A8/A9 row: shift, fdr, n
    for _, r in focus.iterrows():
        tag = f"{r['Gene_name']}/{r['comparison']}"
        sign = "+" if r["delete_shift"] > 0 else "−"
        check(f"{tag} delete_shift", f"{sign}{abs(r['delete_shift']):.5f}", text)
        # Python's %.2e already zero-pads the exponent (1.84e-03); stripping that
        # padding was a bug in this checker, not in the report -- the first run
        # "failed" seven true numbers because the check, not the prose, was wrong.
        check(f"{tag} delete_fdr", f"{r['delete_fdr']:.2e}", text)
        check(f"{tag} delete_n", f"| {r['delete_n']} |", text)

    # ---- S100 family table
    for gene in ("S100A16", "S100A8", "S100A9", "S100A2", "S100A13", "S100A10", "S100A11", "S100A6", "S100A4"):
        check(f"{gene} risk quoted", gene, text)

    # ---- the load-bearing S100 claims
    n_focus_concordant = int(focus["concordant"].sum())
    if n_focus_concordant != 0:
        failures.append(f"focus genes have {n_focus_concordant} concordant rows; report claims none")
    checks += 1
    n_same = int(focus["same_sign_delete_overexpress"].sum())
    check("focus same-sign count", f"In {n_same} of {len(focus)}", text)

    print(f"\n{checks} checks run, {len(failures)} failed")
    for f in failures:
        print("  FAIL:", f)
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
