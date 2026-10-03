# ISP-BAL-1: balanced criteria for Geneformer in silico perturbation

Version 1.0, 2026-10-03. **Proposed, not yet adopted.** It replaces nothing until the cross-evaluation in
[`EXECUTION_PLAN.md`](EXECUTION_PLAN.md) has scored it, and a dated amendment adopts it next to
ISP-STD-1.

## Why a balanced set

The two sources fail in opposite directions:

- **The Geneformer protocol** (Zhang, Venkatesh & Theodoris, *Nature Protocols* 2026; criteria P1–P10 in
  [`CRITERIA_INVENTORY.md`](CRITERIA_INVENTORY.md)) is practical and sensitive. Its gene-level statistic
  treats cells as independent, so with thousands of cells from a few donors, random genes can look
  significant.
- **This repository's criteria** (O1–O19) protect against that, and against study confounding, ambient
  RNA and same-sign hits. They have only been tested on cancer T cells and never against a validated
  target. With few donors they can make a positive impossible: 9 donors cap an exact Wilcoxon test at
  p = 0.0039, so no Holm family above 12 genes can pass.

ISP-BAL-1 keeps what each set does well and binds the two together:

- **The protocol's screen ranks; donor-level tests decide.** The screen (D1) still produces the ranked
  list the protocol is built for, but it promotes a gene only when the run shows it is calibrated (D2).
- **Claims rest on donors, with power checked first.** A claim needs the donor-level test with matched
  controls (D3–D5). That test uses Holm when the family can pass and a labelled BH fallback when it
  cannot, so a design with few donors is weaker but not hopeless.
- **Both directions plus the null.** Deletion and overexpression must disagree (D6), and the effect must
  beat the random-gene null (D7). Opposed shifts in the two arms are already the baseline for random
  genes, so concordance alone is not evidence.
- **The protocol's alternative state counts.** When a design has one, a candidate must move toward the
  goal and not toward the alternative (D6, from P5).
- **Evidence has levels,** from a screen hit to experimental validation (E1), so a ranked gene and a
  replicated finding are never written the same way.

## The criteria

Every criterion has an executable check in [`criteria/isp_criteria.py`](criteria/isp_criteria.py). The
unit tests in [`tests/test_isp_criteria.py`](tests/test_isp_criteria.py) build the failure each one
exists to catch and confirm it is caught. Run them with
`python -m pytest protocol_criteria_e3/tests -q`. The last column says where the cross-evaluation tests
the criterion on real data.

### A. Design (before any model is trained)

| ID | Rule | Threshold | From | Check | Tested on real data |
|---|---|---|---|---|---|
| A1 | Every donor sits in exactly one partition; both samples of a donor stay together | no leaking donor | P1, O3 | `check_donor_disjoint` | all runs |
| A2 | No pair of classes is separated by study (or site, batch) alone: some source must contribute both | every pair bridged | O1 | `check_class_confounding` | July classes (should fail), LUAD, E2, Chaffin |
| A3 | Covariates (disease, age, sex) do not differ across partitions | test p > 0.05 | P1 | `check_confounder_balance` | Chaffin |
| A4 | Equal analysis cells per donor | every donor at the cap | O2 | `check_equal_donor_weight` | all runs |

### B. Classifier (must pass before any perturbation)

| ID | Rule | Threshold | From | Check | Tested on real data |
|---|---|---|---|---|---|
| B1 | Held-out macro F1 with a confusion matrix, 95% CI by bootstrap over donors | lower bound > 1/K | P3 | `classifier_gate` | Phase 1 (LUAD, E2), Phase 3 |
| B2 | Pooled balanced accuracy and an exact sign test over held-out donors | ≥ 0.60, and p ≤ 0.05 with more donors above chance than below | O4 | `classifier_gate` | as B1 |
| B3 | A result within 0.045 of 0.60 is reported as within run-to-run noise; the decision follows the rule | band 0.045 | O5 | `classifier_gate` | Phase 5 seed repeat |
| B4 | A classifier trained on donor-permuted labels fails B1–B2 | must fail | new | `permutation_calibration` | Phase 5 |

### C. Perturbation integrity

| ID | Rule | Threshold | From | Check | Tested on real data |
|---|---|---|---|---|---|
| C1 | A perturbation with no edit gives zero shift | exactly 0 | O6 | `noop_gate` | every run |
| C2 | Shifts agree across precisions or model files | Pearson r ≥ 0.99, Spearman ≥ 0.95, top-20 overlap ≥ 19 | P7, O16 | `precision_agreement` | Phase 5 |
| C3 | One probe agrees across GPU hosts | max abs difference ≤ 1e-3, Spearman ≥ 0.999 | O15 | `cross_host_equivalence` | Phase 5 |
| C4 | State embeddings and perturbation use the same layer, mode and cell filter | identical | P8 | `config_consistency` | every run |

### D. Gene-level evidence

| ID | Rule | Threshold | From | Check | Tested on real data |
|---|---|---|---|---|---|
| D1 | Screen: each gene's per-cell `Shift_to_goal_end` against random genes' shifts, rank-sum test, BH, with `N_Detections` | `Sig` at FDR < 0.05 | P4, P6 | `screen_goal_state_shift` | Phase 1, 3, 4 |
| D2 | Calibration of D1 on random null genes, each screened against the others | fails only if the exact binomial lower 95% bound of the false-positive rate exceeds 0.05 | new | `screen_calibration` | Phase 1 (LUAD, E2), Phase 4 |
| D3 | Per-donor value from the raw per-cell output | mean over the donor's perturbed cells | O7 | `donor_values` | all |
| D4 | Adjusted against matched controls: a = s minus the median of 20 controls matched on detection and rank | ≥ 20 controls, else `NOT_ESTIMABLE_CONTROLS` | O8 | `control_adjust` | all |
| D5 | Exact two-sided Wilcoxon over donors. Holm within the family when its strictest step is attainable (family size ≤ 0.05 / (2 / 2^n)); otherwise BH, labelled `BH_POWER_ADAPTED` | adjusted p ≤ 0.05 | O9, O10, balanced | `donor_test`, `family_correction`, `max_holm_family` | Phase 3 (9 DCM donors) |
| D6 | Deletion and overexpression both significant with opposite signs (`COHERENT`). With an alternative state: toward the goal and not toward the alternative | `COHERENT`; alternative-state shift ≤ 0 | O11, P5 | `dose_concordance`, `alt_state_ok` | all; P5 in Chaffin |
| D7 | The control-adjusted effect beats random genes. The null deletion/overexpression rank correlation is reported for the run | effect above the 95th percentile of null absolute effects | O13 | `beats_null`, `null_pair_correlation` | all |
| D8 | The sign holds across folds | at most 1 fold flips | O14 | `fold_stable` | all |
| D9 | Ambient RNA: flagged or undetermined genes cannot reach SUPPORTED | non-circular leave-one-anchor-out classifier, AUC ≥ 0.8 | O12 | the LUAD ambient pipeline; passed to E1 as `ambient_clear` | tissue contrasts |

### E. Evidence levels and reporting

| Level | Needs |
|---|---|
| `NOT_ESTIMABLE` | the gene cannot be tested (too few donors, controls or tokens); reason stated |
| `OPEN` | tested, no level reached. Never written as "no effect" (O17) |
| `SCREEN` | D1 `Sig` = 1 in a run that passed D2. A ranking for follow-up, not a claim |
| `SUPPORTED` | D5 significant, D6 `COHERENT` (and the alternative-state rule when it applies), D7, D8 and D9 all pass. SCREEN is not required: the screen ranks, the donor test decides |
| `REPLICATED` | SUPPORTED in the same direction in an independent cohort run unchanged (O19) |
| `VALIDATED` | SUPPORTED or REPLICATED, and confirmed experimentally (P10) |

The levels are assigned by `evidence_level`, and `evaluate_genes` applies D1–D8 and E1 to one run's
per-cell table. Every gene keeps a row with its level (O17). Every study is registered with sha256 pins
and a gatekeeper sign-off before its first GPU output, and changed only by dated amendments (O18).

## What changed relative to each source

| Relative to the protocol | Relative to this repository (ISP-STD-1, O1–O19) |
|---|---|
| `Sig` alone is a SCREEN, and only in a calibrated run (D1, D2) | The protocol's screen is kept and reported (D1), not dropped |
| Claims need donor-level, control-adjusted, two-operation evidence (D3–D7) | Holm is replaced by labelled BH where Holm cannot pass with the donors available (D5) |
| Confounding audit of the classes (A2) and equal donor weight (A4) | The protocol's confounder balance (A3), alternative state (P5 in D6), `N_Detections` (D1) and same-filter rule (C4) are added |
| Donor-level classifier gate and permuted-label check (B2, B4) | Macro F1 and the confusion matrix are reported with a donor bootstrap CI (B1) |
| Run-integrity gates: no-op, cross-host (C1, C3) | The protocol's quantization check joins the bf16 canary in one rule (C2) |
| Evidence levels separate a ranked gene from a claim (E1) | Experimental validation (P10) is a level of its own above replication |

## How the set itself is tested

1. **Unit tests** (now). Synthetic data reproduces each failure mode:
   - a leaking donor;
   - study-confounded classes (the July defect);
   - a skewed covariate;
   - a guessing classifier;
   - a nonzero no-op;
   - diverging model sizes;
   - pseudo-replicated null genes, which D2 must flag while independent null genes pass;
   - a same-sign gene;
   - a fold flip;
   - a gene shifted in 3 of 12 donors, which the screen calls and the donor test must not promote;
   - a real, donor-consistent, coherent gene, which must reach SUPPORTED.

   Across 20 seeds, the pseudo-replication, few-donor and real-gene checks gave the intended verdict
   every time.
2. **On real data** (the cross-evaluation, [`EXECUTION_PLAN.md`](EXECUTION_PLAN.md)). ISP-BAL-1 is scored
   beside P and O on the same four properties:
   - calibration, on random null genes and permuted labels;
   - sensitivity, on GSN (validated) and ASCL1/NEUROD1;
   - reproducibility, across seed, host and precision;
   - transfer, from LUAD to colon.

   It is adopted only if two things hold:
   - its false-positive rate on null genes is no higher than O's;
   - it is at least as sensitive: GSN passes under ISP-BAL-1 whenever it passes under O.
