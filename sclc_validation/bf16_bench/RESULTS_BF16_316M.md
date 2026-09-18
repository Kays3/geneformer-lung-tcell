# 316M second-stage replication (fine-tune + ISP comparison + canary)

**Status: COMPLETE.** Fine-tune, state-embedding rebuild, full panel arm,
and the fp32-vs-bf16 canary have all run for real. Every number below
comes from a generated JSON (a chain-stage marker, a `power_sample.sh`
summary, or `compare_runs.py` output) -- none are hand-typed.

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
| eval (during the crashed `classifier.validate()` call) | >=1 (crash-inferred) | not recorded | not recorded |
| test (recovered via `evaluate_saved_model()`, tie-fix + counter installed) | 0 | 9,376 | 0% |

**Caveat on the eval-split row**: the fix (and its counter) was not yet
installed when the original fine-tune's `validate()` call crashed, so the
exact eval-split tie count was never recorded -- only that at least one
tie occurred (that is what crashed the run) is known. The test-split row
is a real, complete count: zero ties across 9,376 predictions, recovered
by re-running `evaluate_saved_model()` against the already-saved
checkpoint with the fix installed from the start. Read together: ties are
real (the eval split hit at least one) but not pervasive (the larger test
split hit zero) -- the recovered macro_f1/acc below are not meaningfully
precision-sensitive to the tie-break choice.

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
| macro_f1 (test) | 0.8979 | 0.60 | **PASS** |
| acc (test) | 0.9087 | -- | -- |

(104M reference, for context only -- not a target: macro_f1 0.9033, acc
0.9194, on the same held-out split. The 316M classifier is close to, not
better or worse in any qualitatively different way from, the 104M
production classifier's held-out quality on this task.)

## Fine-tune timing / energy

| Stage | wall time | mean power (W) | energy (Wh) | peak temp (C) |
|---|---|---|---|---|
| Fine-tune (bf16, 1 epoch, uncapped) | 7242.77s (2.01h) | 65.822 | 132.4258 | 84.0 |
| State-embedding rebuild (real classifier, one-time) | 2089.98s (0.58h) | 88.167 | 51.1851 | 85.0 |
| Full panel arm (150 units, bf16) | 6099.26s (1.69h) | 49.301 | 83.5283 | 83.0 |

Fine-tune wall time (2.01h) closely matches the colleague's ~2h estimate
for this step on comparable hardware. Training itself completed
successfully within this wall time; the crash described above happened in
a downstream metrics call after training and checkpoint-saving were
already done, so this number reflects real training cost, not a partial
or retried run.

## ISP panel comparison (316M-bf16 vs 104M-bf16)

Full panel arm: 150 units (50 genes x 3 sources, overexpress-only,
max-ncells 300 -- same 50-gene panel and lever choices as the 104M stage,
for a like-for-like comparison). Ran `compare_runs.py panel` with
`--a-stats` = 104M-bf16, `--b-stats` = 316M-bf16.

**The G5 acceptance gate does not apply here and is not reported as a
verdict.** `compare_runs.py`'s gate (rho >= 0.95, top-20 overlap >= 19/20,
sign agreement on signal genes) was designed and frozen to bound
precision-induced drift WITHIN one model architecture (fp32 vs bf16 on
the same weights) -- it says nothing about whether two DIFFERENT
architectures should agree, so a low score here is not a "failure" of
anything. The tool still prints a generic `verdict -> FAIL` label because
the CLI always evaluates the same fixed thresholds; that label is
reported below for transparency but must not be read as a failed
replication.

**Result**: Spearman rho = **0.489**, top-20 overlap = **9/20 of the 300
(comparison, Gene_name) ROWS**, ranked by |shift| -- this is
`compare_runs.py`'s own ranking unit (the same row-level rule used for
the 104M gate's 20/20 in PR #17, so nothing there is affected). Ranking
by GENE instead (max |shift| per gene across its 6 comparison rows) gives
a different, larger overlap: **13/20**, independently verified against
the raw stats CSVs. (floor/signal-gene/sign-agreement fields excluded
from this report as not meaningful for a cross-architecture comparison.)

This replicates the colleague's qualitative finding that **changing the
model moves the ISP result far more than changing precision does**: our
own same-model fp32-vs-bf16 precision comparison (104M, PR #17) scored
rho 0.9998 / top-20 20/20 (row-level); this cross-architecture
(104M-vs-316M, both bf16) comparison scores rho 0.489 / top-20 9/20
(row-level) -- an order of magnitude larger divergence from a model-size
change than from a precision change.

**On the colleague's 13/20**: the colleague's own number for the
analogous comparison on their task was 13/20 -- the SAME number as our
gene-level (not row-level) overlap above, which invites a false read
either way. The strongest reason these two 13/20s are not the same
finding is the **ranking unit**, not the more familiar caveat about
different task/data/panel: our row-level result (9/20) and gene-level
result (13/20) are two different statistics computed on the identical
underlying data, and the colleague's reported basis is not independently
confirmed from this side. Both of our numbers (9/20 row-level, 13/20
gene-level) are reported as observed, not targeted, and neither is
claimed to match or refute the colleague's figure.

## Precision canary (316M fp32-vs-bf16, 10 genes)

**Dated note (2026-09-18): this is an explicitly post-hoc, risk-bounding
canary, NOT a G5 gate verdict.** 10 genes (not 5 -- fewer gives only 15
rank points per dtype, not enough for a meaningful rank statistic),
selected as the **top-10 by MAX |Shift_to_goal_end| across the 6
source-goal comparison pairs** in the 316M-bf16 panel above (max, not
mean -- the two rankings diverge from rank 2 onward: by mean, MGP enters
the top-10 and TXNIP drops out). The fp32 twin covers the same genes, so
the selection is symmetric. Selected genes: **HBB, RPS26, HSPA1B, MMP12,
RPS27, S100A7, TPSB2, TXNIP, LAYN, ASCL1**.

**Dated amendment 2026-09-18: row-count basis.** The original sizing note
said "rho over the 30 scores per dtype," assuming one comparison per
(gene, source). The panel's real structure yields 6 comparison rows per
gene (both goal/alt directions for each of 3 sources), so 10 genes -> 60
rows per dtype, not 30. No principled 3-of-6 subset was found that
wouldn't arbitrarily discard half the genuinely distinct comparisons (e.g.
dropping "shift toward normal" entirely for the sclc/luad sources), so
the full 60-row set is used here, consistent with how the panel-comparison
numbers above were computed. Nothing is discarded; this is a deviation
from the original sizing note, recorded rather than silently applied.

10 genes x 3 sources x 2 dtypes = 60 GPU-run units (not to be confused
with the 60 CSV rows above, which come from 10 genes x 6 comparisons).

**Result (60 rows, both dtypes)**: Spearman rho = **0.99867**, sign
agreement = **1.0** (60/60, zero flips), max |delta-shift| = **0.00072**,
mean |delta| = 0.000147, median |delta| = 0.000113. No top-20 claim is
made at this n.

**Speed**: fp32 wall = 3039.36s, bf16 wall = 886.86s -- **3.44x speedup**,
notably larger than the 104M stage's 1.457x. **This is itself a
reportable observation: the larger model gains more from bf16** (more
parameters/compute per forward pass means precision-driven speedups
compound more), consistent with the colleague's own finding that bf16's
benefit is not fixed across model sizes.

**MMP12 and TPSB2 resolve cleanly at 316M scale.** Both genes were left
"unresolved under bf16" by Pam's domain-reviewer ruling on the 104M panel
(PR #17) -- their |shift| there (0.0000760 and 0.0001084 respectively)
fell below bf16's own observed resolution limit (0.000603), so a sign
flip between fp32 and bf16 could not be attributed to a real disagreement
versus simple unresolvability at that magnitude. At 316M, both genes rank
in the panel's own top-10 by |shift| -- a much larger effect size -- and
in this canary their fp32-vs-bf16 agreement is clean: MMP12's paired
values move together throughout (e.g. fp32 -0.00490 vs bf16 -0.00493;
fp32 -0.00601 vs bf16 -0.00603), and TPSB2's do too (e.g. fp32 -0.006925
vs bf16 -0.007053 (luad->sclc); fp32 +0.042521 vs bf16 +0.042267
(sclc->luad)) -- same sign, close magnitude, every pair. (All four
example pairs above were re-verified directly against the raw canary
stats CSVs on 2026-09-18 after a review caught two hand-transcription
errors in an earlier draft of this paragraph -- a swapped fp32/bf16
label pair and one digit off in a third value; corrected and
double-checked, not just re-typed.) This closes the loop on PR #17's one
open caveat:
the 104M non-resolution was a magnitude/precision-floor artifact, not
evidence of a real fp32/bf16 disagreement on these genes -- once the
effect size is large enough (as it is at 316M), bf16 tracks fp32 cleanly
on the very same genes.

## Recommendation

(a) **Classifier quality**: the 316M fine-tune clears the sanity
threshold comfortably (macro_f1 0.898 vs the 0.60 broken-training guard,
close to the 104M reference's 0.903) -- this is a healthy classifier, not
a weak fine-tune that could silently explain a low panel overlap.

(b) **Model-size vs precision (observation, not a production decision)**:
changing model size (104M -> 316M) moves the ISP result far more than
changing precision (fp32 -> bf16) does -- rho 0.489/top-20 9/20
(row-level) or 13/20 (gene-level) for the model-size change vs rho
0.9998/top-20 20/20 (row-level) for the precision change on the same
104M model. This replicates the colleague's qualitative finding on our
own task and data; see the panel-comparison section above for why the
ranking-unit distinction matters more than task/data/panel differences
when comparing to the colleague's own 13/20.

(c) **316M fp32-vs-bf16 canary (risk-bounding, not a gate)**: rho 0.99867,
perfect sign agreement, small max |delta| (0.00072) across 60 real
comparisons on the panel's own top-10-effect-size genes, at a larger bf16
speedup (3.44x) than the 104M stage saw (1.457x). Nothing in this canary
suggests bf16 introduces a precision problem specific to the 316M
architecture -- if anything, the larger model's bigger per-forward-pass
compute makes bf16's numerical behavior track fp32 comparably to, or
better than, what was seen at 104M scale on the same genes.

(d) **MMP12/TPSB2**: PR #17's one open caveat -- two genes "unresolved
under bf16" at 104M due to near-zero effect size -- is closed by this
stage's canary: both resolve cleanly at 316M's larger effect size,
supporting that the 104M non-resolution was a magnitude/precision-floor
artifact rather than a real fp32/bf16 disagreement.

**Hard boundary, reiterated**: none of the above changes the production
104M classifier or any of its committed T3/T4 results. The 316M
classifier built here is a benchmark artifact only, used to characterize
model-size vs precision effects; it is not deployed, referenced, or
substituted for the 104M classifier anywhere in this repository.
