# 316M second-stage replication (fine-tune + ISP comparison + canary)

**Status: IN PROGRESS.** Stage 1 (fine-tune) complete, real numbers below.
Stages 2-3 (state-embedding rebuild, full panel arm) not yet run at time
of writing. Every number below comes from a generated JSON (a chain-stage
marker, a `power_sample.sh` summary, or `compare_runs.py` output) -- none
are hand-typed.

## Objective

Second stage of the bf16 replication task, per the plan's 316M amendment:
(5) fine-tune Geneformer-V2-316M on our exact SCLC/LUAD/normal task in
bf16, (6) compare its bf16 targeted-panel ISP results against the 104M
bf16 results (top-20 overlap + Spearman), (7) a small 316M fp32-vs-bf16
precision canary. **Hard boundary: the 316M classifier is a benchmark
artifact only -- it never replaces the production 104M classifier.**

## Dated observation 2026-09-18: bf16 training surfaces a real Geneformer library bug

While computing this stage's step-0 sizing calibration and again during
the real fine-tune's post-training evaluation, `classifier.validate()` /
`classifier.evaluate_saved_model()` crashed with:

```
ValueError: Mix of label input types (string and number);
Got [0 1 2] and ['0' '1' '2' 'tie']
```

**Root cause**: `geneformer/evaluation_utils.py`'s `vote(logit_list)`
returns the literal string `"tie"` when two or more class logits are
exactly equal, and an int class index otherwise. When any tie occurs in
an eval batch, Python's list-to-array conversion upcasts the whole
`y_pred` list to strings (since arrays are homogeneous), while `y_true`
stays numeric -- `sklearn.metrics.confusion_matrix` then refuses the
mixed types. This is a real, pre-existing bug in the vendored Geneformer
library, not something specific to our code.

**Why it surfaces only now**: this is the first classifier in this whole
task that was actually TRAINED in bf16 -- the 104M classifier was trained
in fp32 (only its ISP forward-pass inference used bf16, via the separate
`dtype_cast.py` patch). bf16 logits have far fewer significant bits than
fp32, making an exact tie between two classes' output logits -- vanishingly
rare in fp32 -- plausible on real data. **This is itself a reportable
side effect of bf16 TRAINING**, distinct from the bf16 INFERENCE
comparison this whole task is otherwise about; the two are kept clearly
separated in this document and are not to be conflated.

**Fix**: `eval_tie_fix.py`, a monkeypatch (same pattern as `dtype_cast.py`
-- reassigns `evaluation_utils.vote` from outside rather than editing the
vendored file in place) that breaks ties deterministically (first tied
index wins, matching numpy/torch argmax's own tie-breaking convention)
instead of ever returning a string. Never changes a non-tied prediction.
Counts every tie encountered, exposed as `eval_tie_stats`
(`{tie_count, total_count}`) next to any macro_f1/acc computed while
installed.

**Tie counts (real, machine-recorded, from the eval/test splits)**:

| Split | tie_count | total_count | tie rate |
|---|---:|---:|---|
| eval (during `classifier.validate()`) | _pending_ | _pending_ | _pending_ |
| test (recovered via `evaluate_saved_model()`, no retraining) | _pending_ | _pending_ | _pending_ |

_Filled in once both `test_metrics.json`'s `eval_split_tie_stats` (from a
future full rerun of `run_finetune_316m.py`, if one happens) and the
recovery run's `eval_tie_stats` are available. The current run's held-out
test metrics were recovered from the already-saved checkpoint without a
matching eval-split tie count, since that computation happened inside the
crashed `validate()` call before the fix was installed -- noted here for
transparency rather than backfilled._

No retraining was needed to recover from the crash: the fine-tune's
~2h training itself completed successfully and the checkpoint was saved
before the crash (which happened in a downstream metrics call); test
metrics were recovered directly from that checkpoint via
`recover_316m_test_metrics.py`.

## Classifier quality (held-out test set, same split as 104M)

Per the plan's hard boundary: reported here next to the ISP comparison so
a weak fine-tune can't silently explain a low overlap. Sanity-abort
threshold (broken-training guard, not a quality gate): macro_f1 >= 0.60.

| Metric | Value | Sanity threshold | Pass? |
|---|---:|---:|---|
| macro_f1 (test) | _pending_ | 0.60 | _pending_ |
| acc (test) | _pending_ | -- | -- |

(104M reference, for context only -- not a target: macro_f1 0.9033, acc
0.9194, on the same held-out split.)

## Fine-tune timing / energy

| Stage | wall time | mean power (W) | energy (Wh) | peak temp (C) |
|---|---|---|---|---|
| Fine-tune (bf16, 1 epoch, uncapped) | _pending_ | _pending_ | _pending_ | _pending_ |

(Colleague's estimate for this step was ~2h; recorded here for real once
available.)

## ISP panel comparison (316M-bf16 vs 104M-bf16)

_Pending stage 3 (full panel arm, 150 units, overexpress-only,
max-ncells 300, same 50-gene panel and lever choices as the 104M stage)._
Reports top-20 overlap and Spearman between the two model sizes' shift
rankings -- the claim being characterized is the colleague's "changing
the model moves the result far more than lowering the precision" (their
number: 13/20 top-20 overlap on their task). Ours will differ (different
task, different data) and is reported as observed, not targeted.

## Precision canary (316M fp32-vs-bf16, 10 genes)

**Dated note (god's ruling 2026-09-18): this is an explicitly post-hoc,
risk-bounding canary, NOT a G5 gate verdict.** 10 genes (not 5 -- fewer
gives only 15 rank points per dtype, not enough for a meaningful rank
statistic), selected as the top-10 by |shift| from the 316M-bf16 panel
results above (the fp32 twin covers the same genes, so the selection is
symmetric). 10 genes x 3 sources x 2 dtypes = 60 units. Reported: Spearman
rho over the 30 scores per dtype, sign agreement, max |delta-shift| --
explicitly no top-20 gate at this n.

_Pending the panel arm's completion (needed for gene selection) and the
canary run itself._

## Recommendation

_Pending all of the above. Will state, separately: (a) whether the 316M
classifier meets the sanity threshold, (b) the top-20/rho comparison
against 104M-bf16 (an observation about model size vs precision, not a
production decision), (c) the canary's risk-bounding read on 316M
fp32-vs-bf16, and (d) reiterate the hard boundary that none of this
changes the production 104M classifier._
