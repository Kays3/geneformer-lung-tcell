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

### T3 — embedding geometry (run; see [RESULTS_T3.md](RESULTS_T3.md))

T3 is implemented in [`embedding_geometry.py`](embedding_geometry.py). It needs the
three aggregate training-donor centroids, which are deliberately stored outside Git
with the larger experiment. On a compute host, export the trusted pickle to a portable,
reviewable `.npz` (three vectors only; no cell- or donor-level records):

```bash
python3 sclc_validation/immune_axis_test/export_state_centroids.py \
  --input /srv/lab/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/state_embeddings/training_donor_disease_centroids.pkl \
  --output sclc_validation/immune_axis_test/results/training_donor_disease_centroids.npz
```

Copy the `.npz` and its adjacent `.json` manifest back to the same path in this
checkout, then run:

```bash
python3 sclc_validation/immune_axis_test/embedding_geometry.py
python3 sclc_validation/immune_axis_test/test_embedding_geometry.py
```

The analysis measures pairwise centroid distances and triangle angles, reconstructs
the unique **in-plane** displacement consistent with each gene's two retained
overexpression scores, summarizes the four programs, and writes the T3 figure. Two
scores cannot identify an out-of-plane displacement, so every vector output records
that assumption. Likewise, centroid geometry alone does not determine T1c's regression
slope without an assumption about displacement covariance; T3 reports full-space and
centroid-plane isotropic references instead of presenting either as a unique prediction.
See [RESULTS_T3.md](RESULTS_T3.md) for the measured triangle and interpretation.

### T4 — program-level overexpression

T4 is implemented but its GPU phases have **not** been run. The deterministic
[`t4_program_manifest.json`](t4_program_manifest.json) contains 4 primary programs,
7 nested exhaustion sets, and 20 expression-matched null sets per program. Rebuild
and validate it on the laptop with:

```bash
python3 sclc_validation/immune_axis_test/build_t4_manifest.py
python3 sclc_validation/immune_axis_test/test_t4_manifest.py
python3 sclc_validation/immune_axis_test/test_t4_runner.py
python3 sclc_validation/immune_axis_test/test_t4_analysis.py
```

Inspect the compute plan without importing Geneformer or touching a GPU:

```bash
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --dry-run --phase program
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --dry-run --phase nested
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --dry-run --phase null
```

On a configured GPU host, run the phases separately so the 12 primary units and
21 nested units finish before the 240 matched-null units. Every unit has its own
completion marker and raw directory, so rerunning the same command resumes:

```bash
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --phase program
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --phase nested
python3 sclc_validation/immune_axis_test/run_t4_overexpression.py --phase null
python3 sclc_validation/immune_axis_test/analyze_t4_overexpression.py
```

Override external paths with `SCLC_PERTURBATION_ROOT`, `HTAN_FINETUNE_ROOT`,
`GENEFORMER_ROOT`, `GENEFORMER_TOKEN_DICT`, and `T4_RUN_DIR`. The analyzer writes
cell-descriptive and donor-level summaries, CD4/CD8/Treg strata, matched-null tests,
and nested-set monotonicity. Normal still has one held-out donor, so its donor-level
uncertainty is not estimable.

The laptop-only ambient/stress/ribosomal sensitivity arm **is complete**:
19 flagged genes were excluded, leaving 31. The 1-D consistency result did not
reverse (1/50 genes before, 1/31 after); PC1 variance changed from 0.740 to 0.727.
See [`results/axis_sensitivity_summary.json`](results/axis_sensitivity_summary.json).

## Status

Tests 1, 2, 3, and 5 (T5a/T5b) run. Test 2 (measured baseline expression) does **not**
trigger PLAN.md's pre-registered falsification criterion (SCLC sits measurably below LUAD
on the exhaustion program, both by detection rate and expression) — see
[RESULTS_T2.md](RESULTS_T2.md) for the full picture, including per-gene
heterogeneity and the single-donor Normal caveat. Test 5 both extends that picture
genome-wide and complicates it: see [RESULTS_T5.md](RESULTS_T5.md) — the curated
exhaustion program's SCLC-vs-LUAD ordering is **not the same** on the complete
population as on T2's test-only population, and the difference traces to donor
composition (a single donor carrying 75% of SCLC's held-out test cells), not a
computation error. Test 3's centroid geometry is complete and rejects the
collinear-centroid premise; see [RESULTS_T3.md](RESULTS_T3.md). Test 4's pipeline
and ambient sensitivity are implemented; its GPU phases remain unexecuted. **No
poster or talk wording has been changed.**
```text
PLAN.md                            the plan; read this first
RESULTS_T2.md                      T2 results
RESULTS_T3.md                      T3 centroid geometry and identifiability results
RESULTS_T5.md                      T5 (T5a/T5b) results
axis_consistency.py                test 1, runs on the laptop
measure_baseline_expression.py     test 2, runs on the compute host (needs the atlas)
analyze_baseline_expression.py     test 2 analysis and figures, runs on the laptop
differential_expression.py         test 5 (T5a/T5b), runs on the compute host (needs the atlas)
analyze_differential_expression.py test 5 analysis and figures, runs on the laptop
check_donor_composition.py         test 5 donor-composition diagnostic, runs on the laptop
export_state_centroids.py          aggregate-vector export, runs on the compute host
embedding_geometry.py              test 3 analysis and figure, runs on the laptop
test_embedding_geometry.py         synthetic T3 geometry regression tests
RESULTS_T3.md                      measured T3 geometry and identifiability narrative
build_t4_manifest.py               deterministic primary/nested/null T4 design
t4_program_manifest.json           pinned T4 run definitions and matched genes
run_t4_overexpression.py           resumable GPU runner with CPU-only dry run
analyze_t4_overexpression.py       donor/subtype/null/titration analysis
axis_consistency_sensitivity.py    T1 ambient/stress/ribosomal sensitivity
test_t4_*.py                       CPU-only T4 regression tests
results/                           output tables for all tests
figures/                           test 2 and test 5 figures
```
