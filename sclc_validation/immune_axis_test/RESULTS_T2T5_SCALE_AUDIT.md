# Results — T2/T5 pseudobulk scale audit

**Status: root cause confirmed; T2 correction regenerated.** This audit
does not change a biological ordering, a gene list, a rank, or a model result. It
identifies an implementation-level normalization mismatch between two summaries of
the same held-out cells and genes.

## Finding

On the prepared H5AD's 9,377 held-out test cells, the current `X` row sum divided by
stored `obs['n_counts']` is 1.017–3.224 (median **1.766**); no cell has equal values.
T2 divided the selected-gene numerator by stored `n_counts`, while T5 used Scanpy
`normalize_total`, which divides by the current full-`X` row sum.

The audit reproduces T2's published cell-weighted seven-gene exhaustion means exactly:
Normal 0.18015, SCLC 0.18634, LUAD 0.24152. Recomputing only the denominator with
the T5 rule yields 0.14524, 0.14818, and 0.19117 respectively: T5/T2 ratios 0.806,
0.795, and 0.792. This is the observed approximately 0.79 scale mismatch.

The reason the stored metadata differs from current `X` is not established by this
audit. It must not be described as a particular preprocessing event without the
object's provenance.

## Decision

**The current full-`X` row sum is canonical.** It is the library size of the matrix
being summarized and matches T5's established Scanpy normalization. T2 now computes
its CP10k denominator from that row sum rather than `obs['n_counts']`.

## Regeneration and historical record

- `baseline_expression_{pooled,per_donor}.csv`, `t2_program_summary.csv`, T2 figures,
  and all T2 log1p(CP10k) values were regenerated from the corrected source. Detection
  rates and raw-count columns are unaffected.
- T6 consumes T2's per-donor log1p table. Its test-only values, weighting table,
  leave-one-donor-out results, figures, and manifest were regenerated after T2. Its
  historical scale-separation assertion was replaced with an equality regression on
  the same test cells and genes.
- T5's existing outputs use the canonical row-sum denominator and need no numerical
  change. Its complete-versus-test donor-composition conclusion remains a separate
  question from this scale fix.

The pre-correction T2/T6 values recorded above are reproducible **legacy-`n_counts`**
outputs; the committed T2/T6 tables use the canonical scale and agree numerically with
matched T5 test cells and genes. The audit data are in
[`results/t2t5_pseudobulk_scale_audit.json`](results/t2t5_pseudobulk_scale_audit.json).
