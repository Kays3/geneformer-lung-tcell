"""Stanley's B5 regression test: apply_noop_gate must parse the REAL noop_spotcheck.py output format --
concatenated json.dump(indent=1) objects with no separator ('}{' glued), no trailing newline between
objects, plus the compact one-line error-fallback format the drivers also append on a failed spot check.
Uses a byte copy of A5's actual results.jsonl (5 real spot checks, all pass=true) plus one appended
compact fallback line for a deliberately-failed 6th gene, exercising both formats in one file."""
import json
import os
import shutil
import sys
import tempfile
import types

SCRIPTS_DIR = sys.argv[1]
REAL_RESULTS = sys.argv[2]      # byte copy of A5's real results.jsonl
NULL_GENES_100 = sys.argv[3]    # byte copy of A5's real null_genes_100_estimable.json
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

import null_analysis as na  # noqa: E402

draw_order = json.load(open(NULL_GENES_100))["draw_order"]
assert len(draw_order) == 100
rows = [{"gene": g, "status": "estimable"} for g in draw_order]

tmp = tempfile.mkdtemp(prefix="noop_real_fmt_")
try:
    # test 1: byte-identical copy of the real file (5 checked, all pass, due == the 5 real gene ids,
    # missing/failed empty) -- proves the streaming raw_decode parser handles the real glued format.
    real_copy = os.path.join(tmp, "results_real.jsonl")
    shutil.copyfile(REAL_RESULTS, real_copy)
    out1 = na.apply_noop_gate({"status": "positive"}, real_copy, rows)
    d1 = out1["no_op_detail"]
    assert len(d1["checked"]) == 5, d1["checked"]
    assert set(d1["due"]) == {draw_order[19], draw_order[39], draw_order[59], draw_order[79], draw_order[99]}
    assert d1["missing"] == [], d1
    assert d1["failed"] == [], d1
    assert out1["status"] == "positive", "no missing/failed -> must not force no_op_failed"
    print("PASS 1: real glued-indent format parses, 5/5 checked+pass, due matches positions 20/40/60/80/100")

    # test 2: same file PLUS one appended compact one-line fallback record for a gene that was NOT one
    # of the 5 due positions (so it doesn't affect due/missing) but IS marked pass=false -- must still
    # parse (mixed indented-glued + compact-with-newline formats in one file) and register the failure.
    mixed_copy = os.path.join(tmp, "results_mixed.jsonl")
    shutil.copyfile(REAL_RESULTS, mixed_copy)
    fallback_gene = draw_order[9]  # position 10, not a due position -- isolates the failed-gene effect
    with open(mixed_copy, "a") as f:
        f.write('{"gene": "%s", "pass": false, "error": "noop_spotcheck.py failed, see run.log"}\n'
                % fallback_gene)
    out2 = na.apply_noop_gate({"status": "positive"}, mixed_copy, rows)
    d2 = out2["no_op_detail"]
    assert len(d2["checked"]) == 6, d2["checked"]
    assert d2["checked"][fallback_gene] is False
    assert d2["failed"] == [fallback_gene], d2
    assert d2["missing"] == [], d2  # the 5 due genes are still all present and passing
    assert out2["status"] == "no_op_failed", "a failed spot check must force no_op_failed"
    print("PASS 2: mixed glued-indent + compact-fallback-line format parses, failure correctly detected")
finally:
    shutil.rmtree(tmp, ignore_errors=True)
