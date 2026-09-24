import glob, pickle, sys, json
import numpy as np
from scipy.stats import spearmanr
P = sys.argv[1]; G = "ENSG00000166710"
def load(run, op):
    fs = glob.glob(f"{P}/{run}/{op}/{G}/*_raw.pickle"); assert len(fs) == 1, fs
    d = pickle.load(open(fs[0], "rb"))
    return {s: np.asarray(next(iter(v.values())), dtype=np.float64) for s, v in d.items()}
def cmp(a, b):
    out = {}
    for s in sorted(a):
        x, y = a[s], b[s]; assert x.shape == y.shape
        d = np.abs(x - y)
        out[s] = {"n": int(x.size), "bitwise_identical": bool(np.array_equal(x, y)), "max_abs_delta": float(d.max()),
                  "n_cells_differing": int((d > 0).sum()), "spearman_rho": float(spearmanr(x, y).correlation)}
    return out
res = {}
for op in ("delete", "overexpress"):
    r = {k: load(k, op) for k in ("ts1_run1", "ts1_run2", "ts2_run1", "ts2_run2")}
    res[op] = {"same_host_ts1": cmp(r["ts1_run1"], r["ts1_run2"]), "same_host_ts2": cmp(r["ts2_run1"], r["ts2_run2"]),
               "cross_host_run1": cmp(r["ts1_run1"], r["ts2_run1"]), "cross_host_run2": cmp(r["ts1_run2"], r["ts2_run2"])}
    # Registered test = ts1_run1 vs ts2_run1; fixed before seeing numbers to require ALL three state shifts (the strictest reading).
    res[op]["registered_PASS"] = all(v["max_abs_delta"] <= 1e-3 and v["spearman_rho"] >= 0.999
                                     for v in res[op]["cross_host_run1"].values())
print(json.dumps(res, indent=1))
