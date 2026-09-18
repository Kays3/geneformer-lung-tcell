"""
Fix for a real Geneformer library bug, found live during the 316M bf16
fine-tune's post-training eval (2026-09-18): `evaluation_utils.vote()`
returns the STRING "tie" when two or more class logits are exactly equal,
otherwise an int class index. When any tie occurs across an eval set,
sklearn.metrics.confusion_matrix crashes:

    ValueError: Mix of label input types (string and number);
    Got [0 1 2] and ['0' '1' '2' 'tie']

because numpy upcasts the whole y_pred array to strings once any element
is a string, while y_true stays int -- a real type-mixing bug in the
library itself, not something specific to our code.

Root cause of WHY this surfaces now: this is the first classifier we have
TRAINED in bf16 (the 104M classifier was trained in fp32; only its ISP
forward passes used bf16, via the separate dtype_cast.py patch). bf16
logits have far fewer significant bits than fp32, so exact ties between
two classes' output logits -- vanishingly rare in fp32 -- become likely
enough to actually occur on a real eval set. This is a genuine, reportable
side effect of bf16 TRAINING (distinct from bf16 INFERENCE, which is what
the rest of this task's precision comparison is about).

Fix: monkeypatch vote() to break ties deterministically (first tied index
wins, matching numpy/torch argmax's own first-index tie-breaking
convention) instead of returning a string. This never changes a
non-tied prediction, and it's committed here as a monkeypatch (not a
patch to the vendored Geneformer source) for the same reason
dtype_cast.py monkeypatches perturber_utils.load_model instead of
patching it in place -- classifier.py resolves `vote` as a bare name
inside evaluation_utils.py's own module globals, so reassigning the
module attribute from outside takes effect at call time.

Counts ties (god's ruling 2026-09-18, item 2): the tie count is the actual
measurement behind the "bf16 training makes exact logit ties
non-negligible" observation -- report `tie_count` / `total_count` next to
any macro_f1/acc computed while this patch is installed, so the sanity
check carries a precision-sensitivity caveat when ties are frequent.
"""
from __future__ import annotations

_installed = False
tie_count = 0
total_count = 0


def install_eval_tie_fix() -> None:
    global _installed
    if _installed:
        return
    from geneformer import evaluation_utils

    def vote_no_ties(logit_list):
        global tie_count, total_count
        total_count += 1
        m = max(logit_list)
        indices = [i for i, x in enumerate(logit_list) if x == m]
        if len(indices) > 1:
            tie_count += 1
        return indices[0]

    evaluation_utils.vote = vote_no_ties
    _installed = True


def reset_tie_counter() -> None:
    global tie_count, total_count
    tie_count = 0
    total_count = 0


def tie_stats() -> dict:
    return {"tie_count": tie_count, "total_count": total_count}
