import hashlib
import json
import os
import shutil
import sys
import tempfile
import types

import numpy as np

SCRIPTS_DIR = sys.argv[1] if len(sys.argv) > 1 else "."
sys.path.insert(0, SCRIPTS_DIR)

# --- stub the `analyse` module exactly like the A5 smoke test did, so gene_medians is exercised
#     (routing logic real) without needing real per-cell pickle data ---
analyse_stub = types.ModuleType("analyse")
analyse_stub.OPS = ("delete", "overexpress")
analyse_stub.GOAL = "goal"


class Design:
    def __init__(self, **kw):
        self.__dict__.update(kw)


analyse_stub.Design = Design
analyse_stub.load_all = lambda *a, **k: {}


# donor_means keyed by (root marker) -- deterministic per gene+op+root, using a hash of gene/op/root
# so different (gene,op) pairs give different medians, and we control which genes are "not_estimable"
# via a small blacklist.
NOT_ESTIMABLE_GENES = {"GENE_NOTEST"}


def fake_donor_means(calls, gene, op, goal, control_min_cells):
    if gene in NOT_ESTIMABLE_GENES:
        return {}
    donors = [f"D{i}" for i in range(12)]
    h = int(hashlib.sha256(f"{gene}|{op}".encode()).hexdigest(), 16)
    base = ((h % 2000) - 1000) / 1000.0
    return {d: (base + 0.001 * i, 20) for i, d in enumerate(donors)}


analyse_stub.donor_means = fake_donor_means
sys.modules["analyse"] = analyse_stub

import null_analysis as na  # noqa: E402
import null_analysis_combined as nac  # noqa: E402

PASS = []
FAIL = []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append((name, detail))
    print(("PASS " if cond else "FAIL ") + name + (f" -- {detail}" if detail and not cond else ""))


tmp = tempfile.mkdtemp(prefix="nac_test_")
try:
    root_a5 = os.path.join(tmp, "phase8_null", "out")
    root_ext = os.path.join(tmp, "phase8_null_ext", "out")
    donors = [f"D{i}" for i in range(12)]

    def write_markers(root, gene, npos_by_donor, op_seconds=1.0):
        for op in ("delete", "overexpress"):
            d = os.path.join(root, op, gene)
            os.makedirs(d, exist_ok=True)
            for donor, npos in npos_by_donor.items():
                json.dump({"gene": gene, "op": op, "donor": donor, "status": "done",
                           "n_token_cells": npos, "seconds": op_seconds},
                          open(os.path.join(d, f"{donor}.complete.json"), "w"))

    npos = {d: 15 for d in donors}
    # A5 genes: 3 real genes (positions 1..3), one of them not-estimable
    a5_genes = ["GENE_A", "GENE_NOTEST", "GENE_B"]
    for g in a5_genes:
        write_markers(root_a5, g, npos)
    # extension genes: 3 more
    ext_genes = ["GENE_C", "GENE_D", "GENE_E"]
    for g in ext_genes:
        write_markers(root_ext, g, npos)

    all_genes = a5_genes + ext_genes
    frozen_npos_by_gene = {g: dict(npos) for g in all_genes}

    # --- test 1: routing -- collect_gene_rows_dual reads A5 genes from root_a5, ext genes from root_ext ---
    positions = {}  # gene_medians ignores positions/n_total via the stubbed load_all
    rows, mism = nac.collect_gene_rows_dual(root_a5, root_ext, set(a5_genes), all_genes, donors,
                                             positions, 100, frozen_npos_by_gene)
    check("routing: no mismatches on matched frozen counts", mism == [], mism)
    check("routing: 6 rows returned", len(rows) == 6, len(rows))
    check("routing: GENE_NOTEST -> not_estimable",
          next(r for r in rows if r["gene"] == "GENE_NOTEST")["status"] == "not_estimable")
    estimable_genes = {r["gene"] for r in rows if r["status"] == "estimable"}
    check("routing: 5 genes estimable (all but GENE_NOTEST)",
          estimable_genes == set(all_genes) - {"GENE_NOTEST"}, estimable_genes)

    # --- test 2: gate-identity halt when a marker's n_token_cells disagrees with frozen ---
    bad_frozen = dict(frozen_npos_by_gene)
    bad_frozen["GENE_C"] = {d: 999 for d in donors}  # deliberately wrong
    rows2, mism2 = nac.collect_gene_rows_dual(root_a5, root_ext, set(a5_genes), all_genes, donors,
                                               positions, 100, bad_frozen)
    check("gate-identity: mismatch detected for tampered GENE_C", len(mism2) > 0, len(mism2))
    check("gate-identity: mismatches reference GENE_C only",
          all(m["gene"] == "GENE_C" for m in mism2), mism2)

    # --- test 3: R-C manifest verification ---
    good_file = os.path.join(root_a5, "delete", "GENE_A", "D0.complete.json")
    good_hash = hashlib.sha256(open(good_file, "rb").read()).hexdigest()
    manifest_path = os.path.join(tmp, "a5_manifest.txt")
    rel = os.path.relpath(good_file, root_a5)
    with open(manifest_path, "w") as f:
        f.write(f"{good_hash}  {rel}\n")
    failures = nac.verify_manifest(manifest_path, root_a5)
    check("R-C: manifest matches unmodified file", failures == [], failures)

    with open(good_file, "a") as f:
        f.write("TAMPERED")
    failures2 = nac.verify_manifest(manifest_path, root_a5)
    check("R-C: manifest detects tampering", len(failures2) == 1 and failures2[0][0] == good_file, failures2)
    # restore
    with open(good_file, "r+") as f:
        content = f.read()
        f.seek(0)
        f.write(content.replace("TAMPERED", ""))
        f.truncate()

    missing_manifest = os.path.join(tmp, "missing_manifest.txt")
    with open(missing_manifest, "w") as f:
        f.write(f"{good_hash}  does_not_exist.json\n")
    failures3 = nac.verify_manifest(missing_manifest, root_a5)
    check("R-C: manifest detects a missing file", len(failures3) == 1 and failures3[0][1] == "missing", failures3)

    # --- test 4: combine_headline, all 6 points ---
    check("s.6.7.3 positive+positive", nac.combine_headline("positive", "positive") == "positive, confirmed at N=200")
    check("s.6.7.4 positive+negative", nac.combine_headline("positive", "negative") == "positive (N=100), not confirmed at N=200")
    check("s.6.7.4 positive+not_estimable", nac.combine_headline("positive", "not_estimable") == "positive (N=100), not confirmed at N=200")
    check("s.6.7.5 negative stays negative", nac.combine_headline("negative", "positive") == "negative (primary); N=200 extension positive, not a registered primary result")
    check("s.6.7.5 negative+negative, no extra note",
          nac.combine_headline("negative", "negative") == "negative (primary)")
    check("s.6.7.5 opposite_direction treated as not-positive",
          nac.combine_headline("opposite_direction", "positive") == "negative (primary); N=200 extension positive, not a registered primary result")
    check("s.6.7.6 primary stopped_not_analysed reported as itself, N=200 ignored",
          nac.combine_headline("stopped_not_analysed", "positive") == "primary stopped_not_analysed")
    check("s.6.7.6 primary no_op_failed reported as itself",
          nac.combine_headline("no_op_failed", "negative") == "primary no_op_failed")
    check("s.6.7.4+6 positive primary, extension stopped_not_analysed named explicitly",
          nac.combine_headline("positive", "stopped_not_analysed") == "positive (N=100), N=200 extension stopped_not_analysed")
    check("s.6.7.5+6 negative primary, extension no_op_failed named explicitly",
          nac.combine_headline("negative", "no_op_failed") == "negative (primary); N=200 extension no_op_failed")

    # --- test 5: no-op merge across two files, due positions over the combined 6-gene order ---
    noop_a5_path = os.path.join(tmp, "noop_a5.jsonl")
    noop_ext_path = os.path.join(tmp, "noop_ext.jsonl")
    with open(noop_a5_path, "w") as f:
        pass  # no due positions in a 3-gene A5 subset here (due is every 20th overall, none in 6 genes)
    with open(noop_ext_path, "w") as f:
        pass
    rows_est = [r for r in rows if r["status"] != "stopped_not_analysed"]
    merged_path = os.path.join(tmp, "merged.jsonl")
    with open(merged_path, "w") as f:
        for p in (noop_a5_path, noop_ext_path):
            f.write(open(p).read())
    primary = {"status": "positive"}
    gated = na.apply_noop_gate(primary, merged_path, rows)
    check("no-op merge: with NOOP_EVERY=20 and only 6 genes, due set is empty, no forced failure",
          gated["status"] == "positive" and gated["no_op_detail"]["due"] == [], gated)

finally:
    shutil.rmtree(tmp, ignore_errors=True)

print(f"\n{len(PASS)} passed, {len(FAIL)} failed")
if FAIL:
    sys.exit(1)
