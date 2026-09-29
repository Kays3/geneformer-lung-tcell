"""Amendment 5 P8: no-op spot check (ISP-STD-1 B.3) -- two independent forward passes on identical
input, never a reused tensor (two separate run_isp.py subprocesses, each loading its own model
instance into a fresh output directory), for a spot-checked gene/donor pair. Shift must be exactly 0.
"""
import argparse
import glob
import json
import os
import pickle
import subprocess
import sys

import numpy as np


def load_pickle(d, op, gene, donor):
    fs = glob.glob(os.path.join(d, op, gene, f"in_silico_{op}_{donor.replace('/', '_')}_cell_embs_dict_*_raw.pickle"))
    assert len(fs) == 1, fs
    return {s: np.asarray(next(iter(v.values())), dtype=np.float64) for s, v in pickle.load(open(fs[0], "rb")).items()}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gene", required=True)
    ap.add_argument("--donor", required=True)
    ap.add_argument("--host", required=True)
    ap.add_argument("--run-isp-args", required=True, help="JSON list: python, run_isp.py, its fixed flags")
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    assert not os.path.exists(a.work), a.work
    os.makedirs(a.work)
    base = json.loads(a.run_isp_args)
    results = {}
    for run in ("run1", "run2"):
        out = os.path.join(a.work, run)
        cmd = base + ["--host", a.host, "--genes", a.gene, "--donors", a.donor, "--out", out]
        subprocess.run(cmd, check=True, stdout=open(os.path.join(a.work, f"{run}.log"), "w"), stderr=subprocess.STDOUT)
        results[run] = {op: load_pickle(out, op, a.gene, a.donor) for op in ("delete", "overexpress")}
    cmp = {}
    for op in ("delete", "overexpress"):
        for state in sorted(results["run1"][op]):
            d = np.abs(results["run1"][op][state] - results["run2"][op][state]).max()
            cmp[f"{op}/{state}"] = float(d)
    passed = all(v == 0.0 for v in cmp.values())
    json.dump({"gene": a.gene, "donor": a.donor, "host": a.host, "max_abs_delta_by_state": cmp, "pass": passed},
               open(a.out, "w"), indent=1)
    print(json.dumps({"gene": a.gene, "donor": a.donor, "pass": passed, "max_delta": max(cmp.values())}))
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
