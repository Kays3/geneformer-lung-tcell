import hashlib
import json
import os
import shutil
import sys
import tempfile
import types

import numpy as np

SCRIPTS_DIR = sys.argv[1]
sys.path.insert(0, SCRIPTS_DIR)

analyse_stub = types.ModuleType("analyse")
analyse_stub.OPS = ("delete", "overexpress")
analyse_stub.GOAL = "goal"


class Design:
    def __init__(self, **kw):
        self.__dict__.update(kw)


analyse_stub.Design = Design
analyse_stub.load_all = lambda *a, **k: {}


def fake_donor_means(calls, gene, op, goal, control_min_cells):
    donors = [f"D{i}" for i in range(12)]
    h = int(hashlib.sha256(f"{gene}|{op}".encode()).hexdigest(), 16)
    base = ((h % 2000) - 1000) / 1000.0
    return {d: (base + 0.001 * i, 20) for i, d in enumerate(donors)}


analyse_stub.donor_means = fake_donor_means
sys.modules["analyse"] = analyse_stub

import null_analysis_combined as nac  # noqa: E402

# Speed only, for this local smoke test: N_BOOT=10,000 x N_PERM_SECONDARY=2,000 inner permutations
# each is ~20M scipy.stats.spearmanr calls (the real per-run cost, acceptable as a one-time CPU-only
# post-processing step, but far too slow for a fixture test run twice -- 200-row + the B6 A5-only
# recompute). Reduce iteration counts on the imported module only; the committed, hash-pinned
# null_analysis.py on disk is untouched (production still runs the real N_PERM=100_000/N_BOOT=10_000).
# Same pattern used for prior Amendment 5 smoke tests (see memory 2026-09-29 12:05Z entry).
nac.na.N_PERM = 500
nac.na.N_PERM_SECONDARY = 50
nac.na.N_BOOT = 50

tmp = tempfile.mkdtemp(prefix="nac_happy_test_")
try:
    donors = [f"D{i}" for i in range(12)]
    root_a5 = os.path.join(tmp, "phase8_null", "out")
    root_ext = os.path.join(tmp, "phase8_null_ext", "out")

    def write_markers(root, gene, npos=15):
        for op in ("delete", "overexpress"):
            d = os.path.join(root, op, gene)
            os.makedirs(d, exist_ok=True)
            for donor in donors:
                json.dump({"gene": gene, "op": op, "donor": donor, "status": "done",
                           "n_token_cells": npos, "seconds": 1.0},
                          open(os.path.join(d, f"{donor}.complete.json"), "w"))

    a5_genes = [f"GA{i}" for i in range(10)]
    ext_genes = [f"GE{i}" for i in range(10)]
    for g in a5_genes:
        write_markers(root_a5, g)
    for g in ext_genes:
        write_markers(root_ext, g)

    all_genes = a5_genes + ext_genes
    frozen = {"rows": [{"gene": g, "estimable": True, "npos_by_donor": {d: 15 for d in donors}} for g in all_genes]}
    frozen_path = os.path.join(tmp, "frozen_200.json")
    json.dump(frozen, open(frozen_path, "w"))

    ng_path = os.path.join(tmp, "null_genes_200.json")
    # a5_gene_set = draw_order[:100], so a5_genes must occupy positions within the first 100 and
    # ext_genes within [100:200] -- pad EACH half separately, don't just append all padding at the end
    # (that bug put ext_genes inside the first-100 slice and made them route to the wrong root).
    draw_order = (a5_genes + [f"PADA{i}" for i in range(90)]
                  + ext_genes + [f"PADB{i}" for i in range(90)])
    json.dump({"n": 200, "draw_order": draw_order}, open(ng_path, "w"))
    # note: main() only asserts len==200 and n==200; pad to 200 for the assertion, but only the real 20
    # have markers -- the padding genes will come back stopped_not_analysed, which is fine for this smoke test.

    design_path = os.path.join(tmp, "design.json")
    json.dump({"donors": donors, "strata": {}}, open(design_path, "w"))

    ovx_a5_path = os.path.join(tmp, "ovx_a5.json")
    ovx_ext_path = os.path.join(tmp, "ovx_ext.json")
    json.dump({"n_total": 100, "positions": {}}, open(ovx_a5_path, "w"))
    json.dump({"n_total": 100, "positions": {}}, open(ovx_ext_path, "w"))

    # R-C manifest: hash every A5 marker file actually written. Paths are relative to phase8_null/
    # (dirname(root_a5)), matching Stanley's real manifest convention (it also covers
    # noop_spotchecks_<host>/ alongside out/, so it can't be rooted at out/ itself).
    manifest_path = os.path.join(tmp, "a5_manifest.txt")
    with open(manifest_path, "w") as f:
        for g in a5_genes:
            for op in ("delete", "overexpress"):
                for donor in donors:
                    p = os.path.join(root_a5, op, g, f"{donor}.complete.json")
                    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
                    rel = os.path.relpath(p, os.path.dirname(root_a5))
                    f.write(f"{h}  {rel}\n")
    # pin the manifest FILE to this fixture's own hash, as if it were the real committed one
    import null_analysis_combined as _nac_patch  # already imported below too; safe to re-import
    _nac_patch.EXPECTED_A5_OUTPUT_MANIFEST_SHA256 = hashlib.sha256(open(manifest_path, "rb").read()).hexdigest()

    noop_a5 = os.path.join(tmp, "noop_a5.jsonl")
    noop_ext = os.path.join(tmp, "noop_ext.jsonl")
    open(noop_a5, "w").close()
    open(noop_ext, "w").close()

    # B6 (Stanley's re-gate-2): --a5-primary-result is now authenticated by a fresh recomputation
    # inside main() from rows[:100] with a freshly-seeded rng(SEED) -- so the fixture must supply a
    # GENUINELY computed value (the same call chain, same inputs) rather than a hand-picked status, or
    # main() will correctly HALT (exit 6) on the mismatch. Reuse the SAME frozen_npos_by_gene content
    # main() itself builds from frozen_200.json, and the SAME first-100-of-draw_order gene list, so this
    # exactly reproduces what main()'s internal B6 recompute will independently arrive at.
    frozen_npos_by_gene = {r["gene"]: r["npos_by_donor"] for r in frozen["rows"]}
    a5_rows_only, a5_mism = nac.na.collect_gene_rows(root_a5, draw_order[:100], donors, {}, 100,
                                                      frozen_npos_by_gene)
    assert not a5_mism, a5_mism
    a5_primary_only, _, _ = nac.na.primary_from_rows(a5_rows_only, np.random.default_rng(nac.SEED))
    if a5_primary_only.get("status") != "not_estimable":
        a5_primary_only = nac.na.apply_noop_gate(a5_primary_only, noop_a5, a5_rows_only)
    a5_primary_path = os.path.join(tmp, "a5_primary.json")
    json.dump({"primary": a5_primary_only}, open(a5_primary_path, "w"))
    print("genuine A5-only recompute status (fixture, should equal what main() reproduces):",
          a5_primary_only["status"])

    out_path = os.path.join(tmp, "out.json")
    rc = nac.main([
        "--phase8-root-a5", root_a5, "--phase8-root-ext", root_ext,
        "--null-genes-200", ng_path, "--frozen-200", frozen_path,
        "--noop-results-a5", noop_a5, "--noop-results-ext", noop_ext,
        "--a5-output-manifest", manifest_path, "--a5-primary-result", a5_primary_path,
        "--design", design_path, "--ovx-index-a5", ovx_a5_path, "--ovx-index-ext", ovx_ext_path,
        "--out", out_path,
    ])
    print("return code:", rc)
    result = json.load(open(out_path))
    print("pair_report:", result.get("pair_report"))
    print("extension_n200 status:", result.get("extension_n200", {}).get("status"))
    n_stopped = sum(1 for r in result["rows"] if r["status"] == "stopped_not_analysed")
    print("n_stopped_not_analysed (should be 180, the padding genes):", n_stopped)
    assert n_stopped == 180
    # B6: the headline's primary status must equal the fixture's genuinely-recomputed A5 status (proves
    # main() reproduced the same authentication value main() itself is required to independently derive
    # -- it is NOT hardcoded, so whichever status the synthetic data happens to produce is the pass
    # condition, not a specific string).
    assert result["pair_report"]["primary_n100_status"] == a5_primary_only["status"], \
        (result["pair_report"]["primary_n100_status"], a5_primary_only["status"])
    assert result["pair_report"]["headline"]  # non-empty, s.6.7 combine_headline always returns a string
    assert rc in (0, 1, 6)  # 6 would mean B6 wrongly rejected a genuinely-matching recompute -- see below
    assert rc != 6, "B6 wrongly HALTED on a fixture built to match main()'s own recompute exactly"
    print("PASS: happy-path main() runs end-to-end, B6 authentication matches, headline computed:",
          result["pair_report"]["headline"])
finally:
    shutil.rmtree(tmp, ignore_errors=True)
