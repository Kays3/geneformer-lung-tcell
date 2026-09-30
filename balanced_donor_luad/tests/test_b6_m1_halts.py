"""Stanley's B6 and m1 (run-package re-gate-2): dedicated negative tests.

B6: --a5-primary-result must be authenticated by recomputation, not trusted as-is. Test 1 supplies a
deliberately WRONG a5_primary.json (status/rho that cannot match a fresh recompute of the same rows)
and requires main() to HALT at exit 6, before any headline is computed.

m1: --ovx-index-a5 and --ovx-index-ext must have disjoint position-key sets. Test 2 supplies two
indexes that share one key and requires main() to HALT at exit 7, before the merge silently overwrites
either entry.
"""
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

nac.na.N_PERM = 500
nac.na.N_PERM_SECONDARY = 50
nac.na.N_BOOT = 50


def write_markers(root, gene, donors, npos=15):
    for op in ("delete", "overexpress"):
        d = os.path.join(root, op, gene)
        os.makedirs(d, exist_ok=True)
        for donor in donors:
            json.dump({"gene": gene, "op": op, "donor": donor, "status": "done",
                       "n_token_cells": npos, "seconds": 1.0},
                      open(os.path.join(d, f"{donor}.complete.json"), "w"))


def build_common_fixture(tmp):
    donors = [f"D{i}" for i in range(12)]
    root_a5 = os.path.join(tmp, "phase8_null", "out")
    root_ext = os.path.join(tmp, "phase8_null_ext", "out")
    a5_genes = [f"GA{i}" for i in range(10)]
    ext_genes = [f"GE{i}" for i in range(10)]
    for g in a5_genes:
        write_markers(root_a5, g, donors)
    for g in ext_genes:
        write_markers(root_ext, g, donors)
    all_genes = a5_genes + ext_genes
    frozen = {"rows": [{"gene": g, "estimable": True, "npos_by_donor": {d: 15 for d in donors}} for g in all_genes]}
    frozen_path = os.path.join(tmp, "frozen_200.json")
    json.dump(frozen, open(frozen_path, "w"))
    draw_order = a5_genes + [f"PADA{i}" for i in range(90)] + ext_genes + [f"PADB{i}" for i in range(90)]
    ng_path = os.path.join(tmp, "null_genes_200.json")
    json.dump({"n": 200, "draw_order": draw_order}, open(ng_path, "w"))
    design_path = os.path.join(tmp, "design.json")
    json.dump({"donors": donors, "strata": {}}, open(design_path, "w"))
    manifest_path = os.path.join(tmp, "a5_manifest.txt")
    with open(manifest_path, "w") as f:
        for g in a5_genes:
            for op in ("delete", "overexpress"):
                for donor in donors:
                    p = os.path.join(root_a5, op, g, f"{donor}.complete.json")
                    h = hashlib.sha256(open(p, "rb").read()).hexdigest()
                    rel = os.path.relpath(p, os.path.dirname(root_a5))
                    f.write(f"{h}  {rel}\n")
    nac.EXPECTED_A5_OUTPUT_MANIFEST_SHA256 = hashlib.sha256(open(manifest_path, "rb").read()).hexdigest()
    noop_a5 = os.path.join(tmp, "noop_a5.jsonl")
    noop_ext = os.path.join(tmp, "noop_ext.jsonl")
    open(noop_a5, "w").close()
    open(noop_ext, "w").close()
    return dict(donors=donors, root_a5=root_a5, root_ext=root_ext, frozen_path=frozen_path,
                draw_order=draw_order, ng_path=ng_path, design_path=design_path,
                manifest_path=manifest_path, noop_a5=noop_a5, noop_ext=noop_ext)


tmp = tempfile.mkdtemp(prefix="nac_b6_m1_test_")
try:
    fx = build_common_fixture(tmp)

    # --- test 1: B6 -- a deliberately WRONG a5_primary.json, cannot match a fresh recompute ---
    a5_primary_path = os.path.join(tmp, "a5_primary_wrong.json")
    json.dump({"primary": {"status": "positive", "rho": -0.999, "p_lower_tail": 0.0001}},
              open(a5_primary_path, "w"))
    ovx_a5_path = os.path.join(tmp, "ovx_a5.json")
    ovx_ext_path = os.path.join(tmp, "ovx_ext.json")
    json.dump({"n_total": 100, "positions": {}}, open(ovx_a5_path, "w"))
    json.dump({"n_total": 100, "positions": {}}, open(ovx_ext_path, "w"))
    out_path = os.path.join(tmp, "out_b6.json")
    rc = nac.main([
        "--phase8-root-a5", fx["root_a5"], "--phase8-root-ext", fx["root_ext"],
        "--null-genes-200", fx["ng_path"], "--frozen-200", fx["frozen_path"],
        "--noop-results-a5", fx["noop_a5"], "--noop-results-ext", fx["noop_ext"],
        "--a5-output-manifest", fx["manifest_path"], "--a5-primary-result", a5_primary_path,
        "--design", fx["design_path"], "--ovx-index-a5", ovx_a5_path, "--ovx-index-ext", ovx_ext_path,
        "--out", out_path,
    ])
    print("B6 test return code:", rc)
    assert rc == 6, f"expected exit 6 (B6 authentication halt), got {rc}"
    result = json.load(open(out_path))
    assert "HALTED" in result["about"] and "recomputation" in result["about"], result
    assert "mismatch_fields" in result and "status" in result["mismatch_fields"], result
    assert "rows" not in result, "B6 halt must happen before a headline/rows payload is written"
    print("PASS B6: main() halts at exit 6 when --a5-primary-result disagrees with a fresh recompute")

    # --- test 2: m1 -- ovx-index-a5 and ovx-index-ext share one position key ---
    ovx_a5_collide = os.path.join(tmp, "ovx_a5_collide.json")
    ovx_ext_collide = os.path.join(tmp, "ovx_ext_collide.json")
    json.dump({"n_total": 100, "positions": {"GA0|D0": 1}}, open(ovx_a5_collide, "w"))
    json.dump({"n_total": 100, "positions": {"GA0|D0": 2}}, open(ovx_ext_collide, "w"))
    # use the genuinely-matching a5_primary from test 1's own rows so this test isolates m1, not B6
    frozen_npos_by_gene = {r["gene"]: r["npos_by_donor"] for r in json.load(open(fx["frozen_path"]))["rows"]}
    a5_rows_only, mism = nac.na.collect_gene_rows(fx["root_a5"], fx["draw_order"][:100], fx["donors"], {},
                                                   100, frozen_npos_by_gene)
    assert not mism
    a5_primary_ok, _, _ = nac.na.primary_from_rows(a5_rows_only, np.random.default_rng(nac.SEED))
    if a5_primary_ok.get("status") != "not_estimable":
        a5_primary_ok = nac.na.apply_noop_gate(a5_primary_ok, fx["noop_a5"], a5_rows_only)
    a5_primary_ok_path = os.path.join(tmp, "a5_primary_ok.json")
    json.dump({"primary": a5_primary_ok}, open(a5_primary_ok_path, "w"))
    out_path2 = os.path.join(tmp, "out_m1.json")
    rc2 = nac.main([
        "--phase8-root-a5", fx["root_a5"], "--phase8-root-ext", fx["root_ext"],
        "--null-genes-200", fx["ng_path"], "--frozen-200", fx["frozen_path"],
        "--noop-results-a5", fx["noop_a5"], "--noop-results-ext", fx["noop_ext"],
        "--a5-output-manifest", fx["manifest_path"], "--a5-primary-result", a5_primary_ok_path,
        "--design", fx["design_path"], "--ovx-index-a5", ovx_a5_collide, "--ovx-index-ext", ovx_ext_collide,
        "--out", out_path2,
    ])
    print("m1 test return code:", rc2)
    assert rc2 == 7, f"expected exit 7 (m1 ovx position-key collision), got {rc2}"
    result2 = json.load(open(out_path2))
    assert "HALTED" in result2["about"] and "collisions" not in result2["about"], result2
    assert result2["n_collisions"] == 1 and result2["example"] == ["GA0|D0"], result2
    print("PASS m1: main() halts at exit 7 on an ovx-index position-key collision between A5 and extension")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
