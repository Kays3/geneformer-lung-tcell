import hashlib
import json
import os
import shutil
import sys
import tempfile
import types

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
analyse_stub.donor_means = lambda *a, **k: {}
sys.modules["analyse"] = analyse_stub

import null_analysis_combined as nac  # noqa: E402

tmp = tempfile.mkdtemp(prefix="nac_main_test_")
try:
    root_a5 = os.path.join(tmp, "phase8_null", "out")
    os.makedirs(root_a5, exist_ok=True)
    dummy_marker = os.path.join(root_a5, "delete", "GENE_A")
    os.makedirs(dummy_marker, exist_ok=True)
    marker_path = os.path.join(dummy_marker, "D0.complete.json")
    json.dump({"gene": "GENE_A", "n_token_cells": 15}, open(marker_path, "w"))

    manifest_path = os.path.join(tmp, "a5_manifest.txt")
    with open(manifest_path, "w") as f:
        f.write("0" * 64 + "  out/delete/GENE_A/D0.complete.json\n")  # deliberately WRONG content hash

    # test 1: pin the manifest FILE to this fixture's own hash (as if it were the real committed one),
    # so we exercise the LISTED-FILE verification path (exit 3), independent of the hardcoded real-world
    # pin (which is tested separately below, unpatched).
    nac.EXPECTED_A5_OUTPUT_MANIFEST_SHA256 = hashlib.sha256(open(manifest_path, "rb").read()).hexdigest()

    out_path = os.path.join(tmp, "out.json")
    rc = nac.main([
        "--phase8-root-a5", root_a5,
        "--phase8-root-ext", os.path.join(tmp, "phase8_null_ext", "out"),
        "--null-genes-200", os.path.join(tmp, "does_not_matter.json"),
        "--frozen-200", os.path.join(tmp, "does_not_matter2.json"),
        "--noop-results-a5", os.path.join(tmp, "noop_a5.jsonl"),
        "--noop-results-ext", os.path.join(tmp, "noop_ext.jsonl"),
        "--a5-output-manifest", manifest_path,
        "--a5-primary-result", os.path.join(tmp, "does_not_matter3.json"),
        "--design", os.path.join(tmp, "does_not_matter4.json"),
        "--ovx-index-a5", os.path.join(tmp, "does_not_matter5.json"),
        "--ovx-index-ext", os.path.join(tmp, "does_not_matter6.json"),
        "--out", out_path,
    ])
    print("return code:", rc)
    assert rc == 3, f"expected exit 3 (R-C halt), got {rc}"
    result = json.load(open(out_path))
    assert "HALTED" in result["about"], result
    assert not os.path.exists(os.path.join(tmp, "does_not_matter.json")), \
        "R-C halt must happen BEFORE any other file is touched/required to exist"
    print("PASS: main() halts with exit 3 on R-C manifest failure, before reading anything else")

    # test 2: the manifest FILE itself doesn't match the hardcoded pin (simulating a swapped/edited
    # manifest) -- must halt at exit 5, before even attempting to verify any listed file.
    nac.EXPECTED_A5_OUTPUT_MANIFEST_SHA256 = "0" * 64  # a hash the fixture manifest cannot match
    out_path2 = os.path.join(tmp, "out2.json")
    rc2 = nac.main([
        "--phase8-root-a5", root_a5, "--phase8-root-ext", os.path.join(tmp, "phase8_null_ext", "out"),
        "--null-genes-200", os.path.join(tmp, "does_not_matter.json"),
        "--frozen-200", os.path.join(tmp, "does_not_matter2.json"),
        "--noop-results-a5", os.path.join(tmp, "noop_a5.jsonl"),
        "--noop-results-ext", os.path.join(tmp, "noop_ext.jsonl"),
        "--a5-output-manifest", manifest_path,
        "--a5-primary-result", os.path.join(tmp, "does_not_matter3.json"),
        "--design", os.path.join(tmp, "does_not_matter4.json"),
        "--ovx-index-a5", os.path.join(tmp, "does_not_matter5.json"),
        "--ovx-index-ext", os.path.join(tmp, "does_not_matter6.json"),
        "--out", out_path2,
    ])
    print("return code 2:", rc2)
    assert rc2 == 5, f"expected exit 5 (manifest-file pin mismatch), got {rc2}"
    result2 = json.load(open(out_path2))
    assert "HALTED" in result2["about"] and "pinned hash" in result2["about"], result2
    print("PASS: main() halts with exit 5 when the manifest FILE itself doesn't match the pinned hash")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
