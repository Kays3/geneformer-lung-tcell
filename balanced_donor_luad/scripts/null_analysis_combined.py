"""Amendment 6 B1: the combined N=200 confirmatory-extension analysis. Leaves null_analysis.py's own
N=100 primary path completely unchanged (imports it as a library, no edits) -- this is a new, separate
entry point per Stanley's run-package gate on run_phase8_null_extension.sh (commit 3807610).

CPU only. Reads genes 1-100 (the registered primary) from phase8_null/out (A5's tree) and genes
101-200 (the confirmatory extension) from phase8_null_ext/out (Amendment 6's tree), using the same
gate-identity check, eligibility rule, and one-sided Spearman permutation test as null_analysis.py,
unmodified.

  (a) genes/roots: null_genes_200_estimable.json's draw_order[:100] read from --phase8-root-a5,
      draw_order[100:200] read from --phase8-root-ext.
  (b) gate-identity: checked against frozen_200_estimable.json for all 200 genes (na.check_gate_identity,
      unchanged), any mismatch HALTS before any status is reported (same as A5's B2).
  (c) no-op merge: --noop-results-a5 and --noop-results-ext are concatenated into one file before
      na.apply_noop_gate (unchanged) computes due-positions over the combined 200-gene run order,
      excluding stopped_not_analysed, every 20th -- exactly Amendment 5 s.5.13's rule, extended.
  (d) ovx indexes: --ovx-index-a5 and --ovx-index-ext are both required; their `positions` dicts are
      merged (disjoint gene sets, so no key collision) and their `n_total` values must agree.
  (e) R-C: --a5-output-manifest (Stanley's sha256 manifest of A5's 17,192 output files -- every
      phase8_null/out marker+pickle and every noop_spotchecks_thinkstation1 file, recorded at A5's
      post-run clearance 2026-09-29T18:10Z, committed as phase8_null/a5_output_manifest.sha256) is
      itself pinned (EXPECTED_A5_OUTPUT_MANIFEST_SHA256 below) so the manifest file cannot be silently
      swapped, THEN verified in full against the real files BEFORE any file under --phase8-root-a5 is
      read. Any mismatch refuses to run (exit 3 or 5) -- A5's output must not have changed since Stanley
      cleared it. The manifest's paths are relative to phase8_null/ (one level above --phase8-root-a5,
      which is the .../out subtree), since it also covers noop_spotchecks_thinkstation1/ alongside out/.
  (f) s.6.7 pair rule: --a5-primary-result (A5's own null_analysis.py output, already produced and
      unaffected by this script) supplies the registered N=100 primary status; combine_headline()
      below implements Amendment 6 s.6.7's 6 points.
  (g) ambient/secondary-percentiles: unchanged mechanisms (ambient_loao_null.py, secondary_percentiles.py
      already take an arbitrary gene list / any null_analysis-shaped result file) -- run them again
      pointed at this script's --out for the combined 200-gene versions; no new code needed there.
"""
import argparse
import hashlib
import json
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import null_analysis as na  # noqa: E402

SEED = na.SEED
SPECIAL_STATUSES = {"stopped_not_analysed", "no_op_failed"}
# R-C: the A5 output manifest file itself is pinned, so it cannot be silently swapped for one that
# always matches. Stanley's manifest, hive/reports/a5_phase8_null_output_manifest_20260929.sha256,
# recorded 2026-09-29T18:10Z after his A5 post-run clearance (17,192 files: every phase8_null/out
# marker+pickle plus noop_spotchecks_thinkstation1).
EXPECTED_A5_OUTPUT_MANIFEST_SHA256 = "b872f8a2375b7f6da42caa9204fe3b011785f1bd13ca85040eae5e4e20860695"


def verify_manifest(manifest_path, root_for_relative):
    """R-C: manifest_path holds `<sha256>  <path>` lines (same format as goals/xfer_sha256.txt).
    Paths are resolved relative to root_for_relative. Returns a list of (path, reason) failures,
    empty if everything matches. Raises if the manifest itself cannot be read."""
    failures = []
    for line in open(manifest_path):
        line = line.rstrip("\n")
        if not line.strip():
            continue
        expected, _, rel = line.partition("  ")
        rel = rel.strip()
        path = os.path.join(root_for_relative, rel) if not os.path.isabs(rel) else rel
        if not os.path.exists(path):
            failures.append((path, "missing"))
            continue
        actual = hashlib.sha256(open(path, "rb").read()).hexdigest()
        if actual != expected:
            failures.append((path, f"hash mismatch: expected {expected}, got {actual}"))
    return failures


def collect_gene_rows_dual(root_a5, root_ext, a5_gene_set, gene_list, donors, positions, n_total,
                            frozen_npos_by_gene):
    """Same logic as null_analysis.collect_gene_rows, but picks the marker root per gene depending
    on whether it is one of A5's original 100 or one of the extension's 100 new genes. Delegates every
    single-item check to na's own unchanged functions."""
    rows = []
    all_mismatches = []
    for g in gene_list:
        root = root_a5 if g in a5_gene_set else root_ext
        if not na.gene_all_markers_present(root, g, donors):
            rows.append({"gene": g, "status": "stopped_not_analysed"})
            continue
        mism = na.check_gate_identity(root, g, donors, frozen_npos_by_gene.get(g, {}))
        if mism:
            all_mismatches.extend(mism)
            continue
        m = na.gene_medians(root, g, donors, positions, n_total, na.CONTROL_MIN_CELLS)
        if m["delete"] is None or m["overexpress"] is None:
            rows.append({"gene": g, "status": "not_estimable", "detail": m})
            continue
        rows.append({"gene": g, "status": "estimable", "delete_median": m["delete"]["median"],
                     "overexpress_median": m["overexpress"]["median"],
                     "delete_n_donors": m["delete"]["n_donors"],
                     "overexpress_n_donors": m["overexpress"]["n_donors"]})
    return rows, all_mismatches


def combine_headline(a5_status, n200_status):
    """Amendment 6 s.6.7, points 1-6, verbatim:
      1/2. the N=100 primary is the registered status; N=200 is a nested extension, never an
           independent replication.
      3. primary positive + N=200 positive -> "positive, confirmed at N=200".
      4. primary positive + N=200 not positive -> "positive (N=100), not confirmed at N=200".
      5. primary not positive -> headline stays "negative (primary)"; a positive N=200 never upgrades
         it, reported separately as "N=200 extension positive, not a registered primary result".
      6. stopped_not_analysed/no_op_failed on EITHER run is reported as that status for that run, and
         does not override or get overridden by the other run's status.
    Point 6 takes precedence for whichever run it applies to: if the PRIMARY run itself is
    stopped_not_analysed/no_op_failed, the headline is that status directly (points 3-5 do not apply,
    since there is no primary statistical outcome to report against). If only the EXTENSION run is one
    of those statuses, points 3-5's wording is used for the primary side, with the extension's status
    named plainly rather than folded into "not confirmed"."""
    if a5_status in SPECIAL_STATUSES:
        return f"primary {a5_status}"
    if a5_status == "positive":
        if n200_status == "positive":
            return "positive, confirmed at N=200"
        if n200_status in SPECIAL_STATUSES:
            return f"positive (N=100), N=200 extension {n200_status}"
        return "positive (N=100), not confirmed at N=200"
    headline = "negative (primary)"
    if n200_status == "positive":
        headline += "; N=200 extension positive, not a registered primary result"
    elif n200_status in SPECIAL_STATUSES:
        headline += f"; N=200 extension {n200_status}"
    return headline


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase8-root-a5", required=True, help="run_phase8_null.sh's --out tree (A5, genes 1-100)")
    ap.add_argument("--phase8-root-ext", required=True, help="run_phase8_null_extension.sh's --out tree (genes 101-200)")
    ap.add_argument("--null-genes-200", required=True, help="null_genes_200_estimable.json (Amendment 6)")
    ap.add_argument("--frozen-200", required=True, help="frozen_200_estimable.json (Amendment 6)")
    ap.add_argument("--noop-results-a5", required=True)
    ap.add_argument("--noop-results-ext", required=True)
    ap.add_argument("--a5-output-manifest", required=True,
                     help="R-C: Stanley's sha256 manifest of every A5 marker/pickle, recorded at A5's "
                          "post-run clearance; verified before anything under --phase8-root-a5 is read")
    ap.add_argument("--a5-primary-result", required=True,
                     help="A5's own null_analysis.py output JSON (unaffected by this script); supplies "
                          "the registered N=100 primary status for the s.6.7 pair rule")
    ap.add_argument("--design", required=True)
    ap.add_argument("--ovx-index-a5", required=True)
    ap.add_argument("--ovx-index-ext", required=True)
    ap.add_argument("--phase6-controls-root", default=None)
    ap.add_argument("--controls-ovx-index", default=None)
    ap.add_argument("--ambient", default=None, help="ambient_loao_null.py output scored for all 200 genes")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    # --- R-C: pin the manifest FILE itself first (it must not be silently swappable), then verify
    #     every file it lists, BEFORE anything under --phase8-root-a5 is read. Paths in the manifest
    #     are relative to phase8_null/ (one level above --phase8-root-a5, the .../out subtree), since
    #     Stanley's manifest also covers noop_spotchecks_thinkstation1/ alongside out/. ---
    manifest_hash = hashlib.sha256(open(a.a5_output_manifest, "rb").read()).hexdigest()
    if manifest_hash != EXPECTED_A5_OUTPUT_MANIFEST_SHA256:
        json.dump({"about": "HALTED: the A5 output manifest file itself does not match the pinned hash -- "
                             "refusing before checking any of its listed files",
                   "expected": EXPECTED_A5_OUTPUT_MANIFEST_SHA256, "actual": manifest_hash},
                  open(a.out, "w"), indent=1)
        print("HALTED: a5-output-manifest file hash mismatch", file=sys.stderr)
        return 5
    manifest_root = os.path.dirname(os.path.normpath(a.phase8_root_a5))
    failures = verify_manifest(a.a5_output_manifest, manifest_root)
    if failures:
        json.dump({"about": "HALTED: R-C manifest verification failed -- A5 output changed since Stanley's "
                             "post-run clearance manifest was recorded",
                   "n_failures": len(failures),
                   "failures": [{"path": p, "reason": r} for p, r in failures]}, open(a.out, "w"), indent=1)
        print(f"HALTED: {len(failures)} R-C manifest failures, see {a.out}", file=sys.stderr)
        return 3

    ng = json.load(open(a.null_genes_200))
    assert ng["n"] == 200 and len(ng["draw_order"]) == 200, "Amendment 6 registers exactly N=200 estimable genes"
    a5_gene_set = set(ng["draw_order"][:100])

    design = na.an.Design(**json.load(open(a.design)))
    ix_a5 = json.load(open(a.ovx_index_a5))
    ix_ext = json.load(open(a.ovx_index_ext))
    if ix_a5["n_total"] != ix_ext["n_total"]:
        json.dump({"about": "HALTED: ovx index n_total disagrees between A5 and the extension",
                   "n_total_a5": ix_a5["n_total"], "n_total_ext": ix_ext["n_total"]}, open(a.out, "w"), indent=1)
        print("HALTED: ovx index n_total mismatch", file=sys.stderr)
        return 4
    positions = {tuple(k.split("|")): v for k, v in ix_a5["positions"].items()}
    for k, v in ix_ext["positions"].items():
        positions[tuple(k.split("|"))] = v
    n_total = ix_a5["n_total"]

    frozen = json.load(open(a.frozen_200))
    frozen_npos_by_gene = {r["gene"]: r["npos_by_donor"] for r in frozen["rows"] if r["estimable"]}

    rng = np.random.default_rng(SEED)

    rows, mismatches = collect_gene_rows_dual(a.phase8_root_a5, a.phase8_root_ext, a5_gene_set,
                                               ng["draw_order"], design.donors, positions, n_total,
                                               frozen_npos_by_gene)
    if mismatches:
        json.dump({"about": "HALTED: gate-identity mismatch (Amendment 6, same rule as A5 s.5.2/B2)",
                   "n_mismatches": len(mismatches), "mismatches": mismatches}, open(a.out, "w"), indent=1)
        print(f"HALTED: {len(mismatches)} gate-identity mismatches, see {a.out}", file=sys.stderr)
        return 2

    primary_200, x, y = na.primary_from_rows(rows, rng)
    if primary_200.get("status") != "not_estimable":
        merged_noop_path = a.out + ".merged_noop.jsonl"
        with open(merged_noop_path, "w") as f:
            for p in (a.noop_results_a5, a.noop_results_ext):
                if os.path.exists(p):
                    f.write(open(p).read())
        primary_200 = na.apply_noop_gate(primary_200, merged_noop_path, rows)

    a5_primary = json.load(open(a.a5_primary_result))["primary"]
    a5_status = a5_primary["status"]
    headline = combine_headline(a5_status, primary_200["status"])

    result = {"about": __doc__.strip(), "n_registered": ng["n"],
              "pair_report": {"primary_n100_status": a5_status, "extension_n200_status": primary_200["status"],
                               "headline": headline},
              "rows": rows, "extension_n200": primary_200}

    if a.phase6_controls_root and a.controls_ovx_index:
        cix = json.load(open(a.controls_ovx_index))
        cpositions = {tuple(k.split("|")): v for k, v in cix["positions"].items()}
        control_ids = sorted({c for st in design.strata.values() for c in st.get("controls", [])})
        crows, _ = na.collect_gene_rows(a.phase6_controls_root, control_ids, design.donors, cpositions,
                                         cix["n_total"], frozen_npos_by_gene=None)
        cprimary, _, _ = na.primary_from_rows(crows, rng)
        validity_failed = not (cprimary.get("rho") is not None and cprimary["rho"] < 0)
        result["validity_check_318_control_genes"] = {
            "rho": cprimary.get("rho"), "n_estimable": cprimary.get("n_estimable"),
            "status": cprimary.get("status"), "sign_compatible_with_minus_0_62": not validity_failed,
        }
        result["validity_check_failed"] = validity_failed
    else:
        result["validity_check_318_control_genes"] = {"status": "SKIPPED"}
        result["validity_check_failed"] = None

    if a.ambient:
        amb = {r["ensembl_id"]: r for r in json.load(open(a.ambient))["rows"]}
        for r in rows:
            r["ambient"] = amb.get(r["gene"], {"ambient_label": "NOT_SCORED"})

    if x is not None:
        result["null_delete_medians"] = x.tolist()
        result["null_overexpress_medians"] = y.tolist()
        result["null_gene_order"] = [r["gene"] for r in rows if r["status"] == "estimable"]

    json.dump(result, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in result.items()
                      if k not in ("rows", "null_delete_medians", "null_overexpress_medians", "null_gene_order")},
                     indent=1))
    return 0 if primary_200.get("status") != "no_op_failed" else 1


if __name__ == "__main__":
    sys.exit(main())
