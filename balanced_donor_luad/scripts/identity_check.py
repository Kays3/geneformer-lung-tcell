"""Amendment 3h.6 GPU identity check, one host's frozen pairs.

For each pair (control gene, donor) assigned to this host in phase7_prep/identity_pairs.json:
  1. write a token-positive-only copy of the donor's ISP input, dataset order preserved
     (the 2026-09-22 construction, where output order is known);
  2. run run_isp.py for that gene and donor on it (same fold model, goals, code path);
  3. compare its overexpress per-cell output with the FULL Phase 6 overexpress pickle at the
     reconstructed positions (isp_order), and its delete output with the Phase 6 delete pickle.
Pass per pair: every state, max|d| <= 1e-3 and Spearman rho >= 0.999. Bit identity is not expected.
Values of these control genes are read as a technical check only.
"""
import argparse
import glob
import json
import os
import pickle
import subprocess
import sys

import numpy as np
from datasets import load_from_disk
from scipy.stats import spearmanr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import isp_order as io  # noqa: E402

TOL, RHO = 1e-3, 0.999


def load_pickle(d, op, donor):
    fs = glob.glob(os.path.join(d, f"in_silico_{op}_{donor}_cell_embs_dict_*_raw.pickle"))
    assert len(fs) == 1, fs
    return {s: np.asarray(next(iter(v.values())), dtype=np.float64) for s, v in pickle.load(open(fs[0], "rb")).items()}


def compare(x, y):
    if x.shape != y.shape:
        return {"n_ref": int(x.size), "n_new": int(y.size), "pass": False, "reason": "length differs"}
    dmax = float(np.abs(x - y).max())
    rho = float(spearmanr(x, y).correlation)
    return {"n": int(x.size), "max_abs_delta": dmax, "spearman_rho": rho, "bitwise_identical": bool(np.array_equal(x, y)),
            "pass": dmax <= TOL and rho >= RHO}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--host", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--phase6", required=True, help="this host's Phase 6 output root")
    ap.add_argument("--work", required=True, help="fresh directory for subset inputs and outputs")
    ap.add_argument("--run-isp-args", required=True, help="JSON list: python, run_isp.py and fixed args")
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    assert not os.path.exists(a.work), f"{a.work} exists"
    os.makedirs(a.work)
    pairs = [p for p in json.load(open(a.pairs))["pairs"] if p["host"] == a.host and not p.get("empty")]
    manifest = json.load(open(a.manifest))
    tok = pickle.load(open(a.token_dict, "rb"))
    base_cmd = json.loads(a.run_isp_args)
    results = []
    for p in pairs:
        g, d = p["gene"], p["donor"]
        wd = os.path.join(a.work, f"{g}_{d}")
        ds = load_from_disk(manifest[d]["isp_input"])
        keep = [i for i, ids in enumerate(ds["input_ids"]) if tok[g] in ids]
        sub_path = os.path.join(wd, "input_token_positive.dataset")
        ds.select(keep).save_to_disk(sub_path)
        m2 = dict(manifest); m2[d] = dict(manifest[d], isp_input=sub_path)
        man_path = os.path.join(wd, "manifest.json"); json.dump(m2, open(man_path, "w"))
        out = os.path.join(wd, "out")
        cmd = base_cmd + ["--manifest", man_path, "--genes", g, "--donors", d, "--out", out]
        with open(os.path.join(wd, "run.log"), "w") as log:
            subprocess.run(cmd, check=True, stdout=log, stderr=subprocess.STDOUT)
        order_ds = ds.add_column("_orig_idx", list(range(len(ds))))
        positions = io.positive_positions(io.replay_order(order_ds), [set(x) for x in ds["input_ids"]], tok[g])
        assert len(positions) == len(keep) == p["k"], (len(positions), len(keep), p["k"])
        safe = d.replace("/", "_")
        full = load_pickle(os.path.join(a.phase6, "overexpress", g), "overexpress", safe)
        new_o = load_pickle(os.path.join(out, "overexpress", g), "overexpress", safe)
        ref_del = load_pickle(os.path.join(a.phase6, "delete", g), "delete", safe)
        new_del = load_pickle(os.path.join(out, "delete", g), "delete", safe)
        ovx = {s: compare(full[s][positions], new_o[s]) for s in sorted(full)}
        dele = {s: compare(ref_del[s], new_del[s]) for s in sorted(ref_del)}
        results.append({**p, "overexpress": ovx, "delete_reported": dele,
                        "pass": all(v["pass"] for v in ovx.values())})
        print(json.dumps({"gene": g, "donor": d, "category": p["category"], "pass": results[-1]["pass"],
                          "ovx_max_abs_delta": max(v.get("max_abs_delta", float("inf")) for v in ovx.values()),
                          "ovx_min_rho": min(v.get("spearman_rho", -1) for v in ovx.values()),
                          "delete_max_abs_delta": max(v.get("max_abs_delta", float("inf")) for v in dele.values())}))
    json.dump({"host": a.host, "tolerance": {"max_abs_delta": TOL, "spearman_rho": RHO}, "pairs": results,
               "PASS": bool(results) and all(r["pass"] for r in results)}, open(a.out, "w"), indent=1)
    return 0


if __name__ == "__main__":
    sys.exit(main())
