"""Step 0.2: study-stratified 5-fold donor cross-fitting assignment (seed 20260924).

Every donor is held out (test) in exactly one fold. Within each fold, the
training donors are all donors of the other folds, from which N_EVAL donors
are set aside as eval (reported only; never used for model selection). The
unit of assignment is the DONOR, so both tissues of a donor always sit in the
same partition. Written before any training; committed as cohort/folds.json.

Stratification: donors are grouped by study; within a study they are sorted
by ID and shuffled with a study-specific seed, then dealt round-robin into
folds, continuing the deal position across studies so fold sizes differ by
at most one.
"""
import argparse
import hashlib
import json

import numpy as np
import pandas as pd

SEED = 20260924
K = 5
N_EVAL = 4


def _rng(tag):
    return np.random.default_rng(int(hashlib.sha256(f"{SEED}|{tag}".encode()).hexdigest()[:16], 16))


def assign_folds(donor_study: pd.Series, k=K):
    """donor_study: index donor_id -> study. Returns dict donor -> fold (0..k-1)."""
    fold_of, pos = {}, 0
    for study in sorted(donor_study.unique()):
        donors = sorted(donor_study[donor_study == study].index)
        order = _rng(f"study|{study}").permutation(len(donors))
        for i in order:
            fold_of[donors[i]] = pos % k
            pos += 1
    return fold_of


def pick_eval(train_donors, donor_study, fold, n_eval=N_EVAL):
    """Deterministic, study-spread choice of eval donors from a fold's training donors."""
    train_donors = sorted(train_donors)
    order = _rng(f"eval|{fold}").permutation(len(train_donors))
    chosen, seen = [], set()
    for i in order:                       # first pass: one per study
        d = train_donors[i]
        if donor_study[d] not in seen:
            chosen.append(d); seen.add(donor_study[d])
        if len(chosen) == n_eval:
            return sorted(chosen)
    for i in order:                       # fill if fewer studies than n_eval
        d = train_donors[i]
        if d not in chosen:
            chosen.append(d)
        if len(chosen) == n_eval:
            break
    return sorted(chosen)


def build(donor_study: pd.Series, k=K):
    fold_of = assign_folds(donor_study, k)
    folds = []
    for f in range(k):
        test = sorted(d for d, x in fold_of.items() if x == f)
        rest = sorted(d for d, x in fold_of.items() if x != f)
        ev = pick_eval(rest, donor_study, f)
        train = sorted(set(rest) - set(ev))
        folds.append({"fold": f, "train": train, "eval": ev, "test": test})
    return {"seed": SEED, "k": k, "n_eval_per_fold": N_EVAL, "unit": "donor",
            "donor_test_fold": {d: int(x) for d, x in sorted(fold_of.items())},
            "donor_study": {d: donor_study[d] for d in sorted(donor_study.index)},
            "folds": folds}


def check(split):
    """Structural invariants; raises on violation."""
    donors = set(split["donor_test_fold"])
    tested = [d for f in split["folds"] for d in f["test"]]
    if sorted(tested) != sorted(donors):
        raise ValueError("every donor must be tested exactly once")
    for f in split["folds"]:
        tr, ev, te = set(f["train"]), set(f["eval"]), set(f["test"])
        if tr & ev or tr & te or ev & te:
            raise ValueError(f"fold {f['fold']}: partitions overlap")
        if tr | ev | te != donors:
            raise ValueError(f"fold {f['fold']}: partitions do not cover all donors")
    sizes = [len(f["test"]) for f in split["folds"]]
    if max(sizes) - min(sizes) > 1:
        raise ValueError(f"unbalanced fold sizes {sizes}")
    return True


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cohort", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    c = pd.read_csv(a.cohort)
    ds = c.drop_duplicates("donor_id").set_index("donor_id").study
    split = build(ds)
    check(split)
    json.dump(split, open(a.out, "w"), indent=1)
    for f in split["folds"]:
        print(f"fold {f['fold']}: train {len(f['train'])} eval {len(f['eval'])} test {len(f['test'])} "
              f"test studies {pd.Series([ds[d] for d in f['test']]).value_counts().to_dict()}")


if __name__ == "__main__":
    main()
