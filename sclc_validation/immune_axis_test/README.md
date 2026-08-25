# Immune axis test — is Normal < SCLC < LUAD a real ordering?

The poster and talk carry one unresolved claim: the model orders the three
states **Normal < SCLC < LUAD** on the checkpoint axis, which appears to
contradict the clinical picture of SCLC as an immune-cold, ICI-resistant tumour.
That claim rests on two shifts measured from a single source state.

This module tests it with the **overexpression arm**, which is census-complete
(every held-out source cell perturbed for every gene, N constant per state) and
therefore the only arm in which the three source states are directly comparable.

[**PLAN.md**](PLAN.md) — the full plan: what is being tested, what has been run,
the three analyses still needed, the replacement figure spec, and the results
that would overturn the conclusion.

## What runs today

```bash
python3 sclc_validation/immune_axis_test/axis_consistency.py
```

Laptop only — reads the committed stats tables from
[`../perturbation_workflow/targeted_panel/results/`](../perturbation_workflow/targeted_panel/results),
needs no GPU and no raw perturbation pickles. Writes six tables to `results/`.

Headline: **1 gene of 50 is consistent with a one-dimensional ordering**, and for
the seven-gene exhaustion program the "SCLC in the middle" prediction fails from
the Normal source for every gene. See PLAN.md §4 for the full result and §5 for
the working interpretation.

```bash
# on the compute host, needs the prepared H5AD:
python3 sclc_validation/immune_axis_test/measure_baseline_expression.py
# back on the laptop, needs the two CSVs measure_baseline_expression.py wrote:
python3 sclc_validation/immune_axis_test/analyze_baseline_expression.py
```

Headline: SCLC sits measurably **below** LUAD on the exhaustion program (both
detection rate and expression), so T2 does **not** trigger the pre-registered
falsification criterion — but it also does not confirm the specific "SCLC ≲
Normal" claim: measured exhaustion is Normal < SCLC < LUAD, the poster's
original ordering, reached by direct measurement rather than the invalid
two-shift argument. See [RESULTS_T2.md](RESULTS_T2.md) for the full picture.

```bash
# on the compute host, needs the prepared H5AD -- run once per POPULATION:
POPULATION=complete  python3 sclc_validation/immune_axis_test/differential_expression.py  # T5a
POPULATION=test_only python3 sclc_validation/immune_axis_test/differential_expression.py  # T5b
# back on the laptop, needs both populations' output pulled back:
python3 sclc_validation/immune_axis_test/analyze_differential_expression.py
python3 sclc_validation/immune_axis_test/check_donor_composition.py
```

**Run — see [RESULTS_T5.md](RESULTS_T5.md).** The first real run surfaced a gene-identity
bug the synthetic-fixture smoke tests couldn't catch (DE output was keyed on Ensembl IDs,
downstream analysis matched on gene symbols); fixed in `differential_expression.py` and
re-run. Headline: 4 of 7 curated exhaustion genes rank in the top 20% of genome-wide DE
significance (SCLC vs LUAD); ISP hits are 10–21× enriched for genome-wide DE significance
in all three state pairs; and the curated exhaustion program's SCLC-vs-LUAD ordering
**reverses** between the complete population and T2's test-only population — traced to
SCLC's held-out test split being 75% one atypically low-exhaustion donor, not a pipeline
error. See RESULTS_T5.md for the full picture.

## Status

Tests 1, 2, and 5 (T5a/T5b) run. Test 2 (measured baseline expression) does **not**
trigger PLAN.md's pre-registered falsification criterion (SCLC sits measurably below LUAD
on the exhaustion program, both by detection rate and expression) — see
[RESULTS_T2.md](RESULTS_T2.md) for the full picture, including per-gene
heterogeneity and the single-donor Normal caveat. Test 5 both extends that picture
genome-wide and complicates it: see [RESULTS_T5.md](RESULTS_T5.md) — the curated
exhaustion program's SCLC-vs-LUAD ordering is **not the same** on the complete
population as on T2's test-only population, and the difference traces to donor
composition (a single donor carrying 75% of SCLC's held-out test cells), not a
computation error. Test 3 (embedding geometry) and 4 (program-level overexpression,
GPU) remain specified but not executed. **No poster or talk wording has been
changed** — see PLAN.md §10, which requires both T2 and T3 before that, and T5's
new disagreement is an additional reason to wait for T3 rather than less.

```text
PLAN.md                            the plan; read this first
RESULTS_T2.md                      T2 results
RESULTS_T5.md                      T5 (T5a/T5b) results
axis_consistency.py                test 1, runs on the laptop
measure_baseline_expression.py     test 2, runs on the compute host (needs the atlas)
analyze_baseline_expression.py     test 2 analysis and figures, runs on the laptop
differential_expression.py         test 5 (T5a/T5b), runs on the compute host (needs the atlas)
analyze_differential_expression.py test 5 analysis and figures, runs on the laptop
check_donor_composition.py         test 5 donor-composition diagnostic, runs on the laptop
results/                           output tables for all tests
figures/                           test 2 and test 5 figures
```
