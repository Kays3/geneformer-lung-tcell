"""Phase 4 classifier gate (registration s.2): must PASS before any ISP.

Pooled over the 5 folds' held-out donors (each scored by the fold model that
never saw it), on the analysis cells (exactly 100 tumour + 100 normal per donor):
  (i)  pooled balanced accuracy >= 0.60, AND
  (ii) exact two-sided sign test over donors of "per-donor balanced accuracy
       > 0.5": p <= 0.05 with MORE donors above 0.5 than below (a two-sided
       test significant in the wrong direction cannot pass). Donors at exactly
       0.5 are dropped (registered sign-test convention).
Per-fold values and tie counts are reported; a fold below chance is flagged
(reported, not a separate gate). Exit status is non-zero on FAIL.
"""
import argparse
import json
import os
import pickle
import sys
from fractions import Fraction
from math import comb

import numpy as np
import pandas as pd

GATE_BA = 0.60
ALPHA = Fraction(5, 100)


def sign_test_two_sided(n_pos, n_neg):
    n = n_pos + n_neg
    if n == 0:
        return Fraction(1)
    k = min(n_pos, n_neg)
    tail = sum(comb(n, i) for i in range(k + 1))
    return min(Fraction(1), Fraction(2 * tail, 2 ** n))


def balanced_accuracy(y_true, y_pred, classes):
    recalls = [np.mean(y_pred[y_true == c] == c) for c in classes if np.any(y_true == c)]
    return float(np.mean(recalls))


def load_fold(work, k):
    rec = json.load(open(os.path.join(work, f"fold{k}", "fold_record.json")))
    pdict = pickle.load(open(rec["pred_dict"], "rb"))
    id_class = pickle.load(open(rec["id_class_dict"], "rb"))
    meta = pdict["prediction_metadata"]
    df = pd.DataFrame({"donor": meta["individual"], "cell_id": meta["cell_id"],
                       "y_true": [id_class[i] for i in pdict["label_ids"]],
                       "y_pred": [id_class[i] for i in pdict["pred_ids"]]})
    df["fold"] = k
    return df, rec



NONDET_BAND = 0.045   # Amendment 3e.2: largest per-donor BA change in the same-seed fold-0 repeat
NONDET_PROVENANCE = {"per_donor_max_abs_delta_observed": 0.045, "pooled_delta_observed": 0.0017,
                     "band_choice": "per-donor max, the conservative choice; pooled observed is ~26x smaller",
                     "basis": "one same-seed repeat of fold 0 (9 donors); observed magnitudes, not bounds"}


def nondeterminism_margins(pooled, per_donor):
    """Amendment 3e.2 reporting fields. They never change PASS; they say whether it sat inside the noise."""
    pooled_within = abs(pooled - GATE_BA) <= NONDET_BAND
    flippable = per_donor[(per_donor - 0.5).abs() <= NONDET_BAND]
    firm_pos = int(((per_donor > 0.5) & ~per_donor.index.isin(flippable.index)).sum())
    firm_neg = int(((per_donor < 0.5) & ~per_donor.index.isin(flippable.index)).sum())
    nf = len(flippable)
    worst = (firm_pos, firm_neg + nf)
    best = (firm_pos + nf, firm_neg)
    ok = lambda pos, neg: sign_test_two_sided(pos, neg) <= ALPHA and pos > neg
    sign_within = ok(*worst) != ok(*best)
    return {"nondet_band": NONDET_BAND, "nondet_band_provenance": NONDET_PROVENANCE, "pooled_within_nondet_of_threshold": bool(pooled_within),
            "n_flippable_donors": nf,
            "sign_test_p_worst_case": float(sign_test_two_sided(*worst)),
            "sign_test_p_best_case": float(sign_test_two_sided(*best)),
            "sign_test_within_nondet": bool(sign_within),
            "clean_result": bool(not pooled_within and not sign_within)}

def evaluate(df):
    classes = ["tumor_primary", "normal_adjacent"]
    pooled = balanced_accuracy(df.y_true.to_numpy(), df.y_pred.to_numpy(), classes)
    per_donor = df.groupby("donor").apply(
        lambda g: balanced_accuracy(g.y_true.to_numpy(), g.y_pred.to_numpy(), classes), include_groups=False)
    n_pos, n_neg = int((per_donor > 0.5).sum()), int((per_donor < 0.5).sum())
    p = sign_test_two_sided(n_pos, n_neg)
    passed = pooled >= GATE_BA and p <= ALPHA and n_pos > n_neg
    per_fold = df.groupby("fold").apply(
        lambda g: balanced_accuracy(g.y_true.to_numpy(), g.y_pred.to_numpy(), classes), include_groups=False)
    return {"pooled_balanced_accuracy": pooled, "gate_ba": GATE_BA,
            "n_donors": int(len(per_donor)), "donors_above_0.5": n_pos, "donors_below_0.5": n_neg,
            "donors_at_0.5": int((per_donor == 0.5).sum()),
            "sign_test_p_two_sided": str(p), "sign_test_p_float": float(p),
            "per_fold_balanced_accuracy": {int(k): float(v) for k, v in per_fold.items()},
            "folds_below_chance": [int(k) for k, v in per_fold.items() if v < 0.5],
            "per_donor_balanced_accuracy": {d: float(v) for d, v in per_donor.items()},
            "PASS": bool(passed), **nondeterminism_margins(pooled, per_donor)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    frames, recs = zip(*[load_fold(a.work, k) for k in range(5)])
    df = pd.concat(frames)
    res = evaluate(df)
    res["n_cells"] = int(len(df))
    res["tie_stats"] = {r["fold"]: {"eval": r["eval_tie_stats"], "test": r["test_tie_stats"]} for r in recs}
    res["fold_seconds"] = {r["fold"]: r["total_seconds"] for r in recs}
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k != "per_donor_balanced_accuracy"}, indent=1, default=str))
    sys.exit(0 if res["PASS"] else 1)


if __name__ == "__main__":
    main()
