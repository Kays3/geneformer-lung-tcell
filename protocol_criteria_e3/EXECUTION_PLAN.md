# Execution plan: cross-evaluating the protocol's criteria and ours

**Status: plan, not registered, nothing run.** Written 2026-10-03. It extends
[`E3_DESIGN.md`](E3_DESIGN.md) from "the protocol's criteria in T cells of the protocol's atlas" to a full
cross-evaluation: each set of criteria is judged on each other's data. The criteria are listed with
IDs in [`CRITERIA_INVENTORY.md`](CRITERIA_INVENTORY.md): P1 to P10 from the Geneformer protocol (Zhang,
Venkatesh & Theodoris, *Nature Protocols* 2026), O1 to O19 from this repository.

## 1. The question

Since July this repository has built a set of rules for reading a Geneformer in silico perturbation
(ISP) screen (O1 to O19). They were written because protocol-style results on our T-cell data did not
hold up:

- **The July screen.** It used the protocol's `goal_state_shift` statistics (P4) and turned up
  same-sign hits such as S100A8/A9. Its classes also turned out to be study-confounded.
- **Lung and colon under our rules.** The classifier criterion held. Random genes already gave opposed
  deletion and overexpression shifts. The lung candidate genes did not carry to colon.

That record has two gaps:

1. Our rules have only been tested on our own data, in cancer T cells, which the protocol itself names
   as a weak setting. They may be too strict and reject real effects. Nothing has yet checked them
   against a target known to work.
2. The protocol's rules have never been scored on our data with the same controls we apply to ours. So
   the lung and colon results do not yet say how often the protocol's own statistic is wrong.

The plan closes both gaps with a 2 × 2 cross-evaluation:

| | Our T-cell cohorts (LUAD 43 donors, colon E2 19 donors) | Protocol's atlas (Chaffin 2022 hearts: cardiomyocytes, then T cells) |
|---|---|---|
| **Protocol's criteria (P)** | **Cell A**: not done cleanly. The July screen used P4 on confounded classes. Phase 1 re-scores the existing LUAD and E2 outputs with P on CPU | **Cell C**: the protocol's own result. Phase 3 reproduces it on cardiomyocytes, and Phase 4 extends it to T cells |
| **Our criteria (O)** | **Cell B**: done (LUAD, E2 final reports) | **Cell D**: never done. Phase 3 asks whether the protocol's validated target, GSN, survives our criteria; Phase 4 does the same in T cells |

Cell D in cardiomyocytes is the sensitivity test our criteria have never had. GSN was validated in the
wet lab (P10), so if our criteria reject it, they are too strict for that design. The plan then
reports which component rejects it.

## 2. What "evaluating a criterion" means

Each criterion is scored on four properties, with a benchmark fixed in advance for each.

| Property | Question | Benchmark | Pass |
|---|---|---|---|
| Calibration | Does it call things that should not be called? | Random null genes (100 per cohort, already run in LUAD and E2; a fresh draw in Chaffin). For classifier criteria, a donor-level label permutation | Empirical positive rate on null genes ≤ 0.05, upper exact binomial 95% bound reported; permuted-label classifier fails the criterion |
| Sensitivity | Does it call what should be called? | GSN in DCM cardiomyocytes, toward nonfailing (validated, P10). ASCL1 and NEUROD1 in SCLC to LUAD T cells (our internal positive control, `sclc_validation/perturbation_workflow`) | The positive passes |
| Reproducibility | Does the verdict survive a repeat? | Seed repeat of one fold; second GPU host; bf16, fp32 and int8 on probe genes; 104M against 316M | Verdict unchanged; noise measured, not assumed |
| Transfer | Does a hit in one cohort hold in another? | LUAD to colon (both criteria sets); cardiomyocytes to T cells only for the classifier criteria | Sign agreement above chance, exact binomial test |

A criterion can pass calibration and fail sensitivity; that is the expected trade-off, and the scorecard
reports both instead of a single grade.

## 3. Phases

The phases are ordered by dependency. Phase 1 needs no download and no GPU, so it can start as soon as
Phase 0 is signed off. It also gives the first results.

### Phase 0: freeze (CPU, about 1 day; needs the human)

- **Copy and pin the standard.** Copy ISP-STD-1 v1.1 from `hive/standards/` into this folder and pin it
  by sha256.
- **Freeze the inventory.** [`CRITERIA_INVENTORY.md`](CRITERIA_INVENTORY.md) becomes the registered
  list.
- **Write the registration** from this plan and `E3_DESIGN.md`, with the predictions in s.5 written in.
  It needs a gatekeeper sign-off before any output.
- **The human's decisions:**
  - GO for Phase 1 (CPU, existing outputs);
  - GO for the Chaffin download (size and destination stated first);
  - the GPU budget for Phases 3 to 5 (s.6).

### Phase 1: the protocol's criteria on our existing outputs (Cell A; CPU only, about 2 to 3 days)

- **Inputs.** The raw per-cell ISP outputs of the LUAD Phase 6 and Amendment 5 runs and of E2, on
  thinkstation1 and thinkstation2 at the paths in their run records. Also the held-out predictions of
  their fold classifiers.
- **What is new.** No new perturbation. Every number comes from outputs that already exist, so the two
  criteria sets are compared on identical cells and genes.

| Step | What | Criteria |
|---|---|---|
| 1a | For every gene run in LUAD and E2 (Panel B, matched controls, null genes): cell-pooled `Shift_to_goal_end`, Wilcoxon against the shifts of random genes, BH, `Sig`, with `N_Detections` | P4, P6 |
| 1b | Calibration of P4. Each null gene is tested against the other 99, and the share with `Sig` = 1 is P4's false-positive rate. The same null genes go through O7 to O11 against their matched controls, giving ours | P4 against O7–O11 |
| 1c | Classifier criteria on the same predictions: macro F1 and the confusion matrix next to the O4 gate and the O5 band | P3 against O4, O5 |
| 1d | Agreement on Panel B and controls: a 2 × 2 table of P4 `Sig` against an O positive status, with Cohen's kappa | P4 against O10, O11 |
| 1e | Transfer under each set: the LUAD to colon sign agreement of P4 `Sig` genes, next to the O result (3 of 10) | P4, O19 |
| 1f | The July all-gene table (`allgene_delete_overexpress_shift.csv`, 82,702 rows): of the genes significant in both arms, the share that are same-sign (O11) or ambient-flagged (O12) | O11, O12 applied to P4 output |

- **Registered deviation.** The protocol's random-gene reference is every gene in the screen
  (`genes_to_perturb="all"`). Here it is the 100 random null genes, and as a sensitivity, the null genes
  pooled with the controls. Both are stated with every P4 number.
- **Not applicable.** P5 (alternative state) does not apply to two-class LUAD and E2.
- **Output.** `phase1_results/` tables and a short report.

### Phase 2: the Chaffin atlas, Stage 0 (CPU, about 1 day after the download)

As in `E3_DESIGN.md` s.3:

- the download with sha256;
- donor and nucleus counts per disease class for cardiomyocytes and for T cells (or lymphocytes);
- the feasibility rule (GO, GO-REDUCED or NO-GO), fixed before the counts are read;
- the cohort draw, the folds, tokenisation, matched controls, the null-gene draw and the run order, all
  committed and pinned before any GPU call, as in E2.

### Phase 3: cardiomyocytes, the protocol's result under both sets (Cells C and D; GPU)

- **3a. Reproduction (P1–P3, P4, P5, P7).** Follow the protocol on both of its models, because its
  numbers come from two different models:
  - **Single-task model:** NF, HCM and DCM classification of cardiomyocytes by held-out patients,
    reported as macro F1 (protocol: 0.85).
  - **Multi-task model** (cell type and disease, on cardiomyocytes, fibroblasts and macrophages): the
    deletion ISP of DCM cardiomyocytes, goal NF, alternative HCM, `goal_state_shift`. This is the
    protocol's source of Fig. 2b and its GSN result (protocol: 99% macro F1 on disease, 85% on cell
    type).
  - A quantized rerun for P7, on the multi-task model as in the protocol.
  The gate is E3's H3a: macro F1 ≥ 0.80 and GSN with `Sig` = 1 toward NF. If it fails, the plan stops
  at `PIPELINE_MISMATCH`.
- **3b. Our criteria on the same fine-tuned folds (O3–O14),** for both models:
  - the per-donor classifier gate;
  - the no-op gate;
  - the null-gene test (O13);
  - donor-level, control-adjusted, two-operation tests for GSN and for the genes named in the protocol's
    Fig. 2b.
  The Fig. 2b genes are CADPS2, FBN2, HOPX, GPC5, SPDYE2, FGF12, ANGPTL4, GRIP1, MYH6, HMGB1, ESRRG,
  FAM118A, H3F3B, NOVA1, SNX17, ARL17B, TMEM176B, ITM2B, IL31RA, ARSJ and GSN.
- **3c. Which component decides.** For GSN and each Fig. 2b gene, the verdict is recorded after each
  of our criteria is added in turn:
  1. cell-pooled (P4);
  2. donor-level (O7);
  3. control-adjusted (O8);
  4. both operations (O11);
  5. fold-stable (O14).
  The step where a gene drops out is the answer to "which of our criteria disagrees with the protocol,
  and on which gene".

**The between-donor design forces two adaptations, registered before any output.** Disease in the
Chaffin atlas is a between-donor label, so our same-donor goal (O2, O7) cannot apply.

- **Goal state.** The goal is the NF centroid of the training donors in each fold; S1 in the LUAD
  registration used the same idea.
- **Power.** With 9 DCM donors, the smallest two-sided exact Wilcoxon p is 2/2^9 = 0.0039. Holm can then
  reach significance only in families of at most 12 genes, so O10 is underpowered by construction for
  any longer list. The plan therefore registers:
  - GSN as a single pre-specified test (family size 1);
  - the 20 other Fig. 2b genes as a separate family under BH, labelled as a deviation from O10;
  - the eligibility rule O9 recomputed for the actual donor counts.

Without this, our criteria would reject everything in the protocol's atlas for a reason that has nothing
to do with biology.

### Phase 4: T cells of the same hearts (Cells C and D; GPU)

E3's H3b to H3e as designed, on the Stage 0 cohort:

- the classifier under P3 and O4;
- arm P (`goal_state_shift` on all genes);
- arm S (our criteria on the P4 `Sig` = 1 genes, a detection-matched set of non-significant genes, and
  100 null genes).

Phase 1's calibration is repeated here (P4 against O on null genes), so the false-positive rates come
from three cohorts. If Stage 0 returns NO-GO for T cells, Phase 4 is skipped, the counts are reported,
and the cross-evaluation rests on Phases 1 and 3.

### Phase 5: reproducibility block (GPU, small)

Run on Chaffin cardiomyocytes, and on T cells if Phase 4 runs:

- a same-seed repeat of one fold (O5);
- the cross-host probe (O15);
- bf16, fp32 and int8 shifts on 20 probe genes (O16, P7);
- a donor-level label-permutation fine-tune for classifier calibration;
- 104M against 316M on the probe genes.

### Phase 6: synthesis (CPU)

- **The scorecard.** Every P and O criterion is scored on calibration, sensitivity, reproducibility and
  transfer, with the numbers and the benchmark behind each cell.
- **A report in this repository's IMRaD format.**
- **An updated graphical abstract,** built from the result files as now.
- **A dated amendment to ISP-STD-1** proposing changes, if the scorecard shows that one of our criteria
  rejects a validated target or lets null genes through.

## 4. Order and dependencies

```
Phase 0 (freeze, human GOs)
   ├── Phase 1 (CPU, existing outputs) ─────────────────┐
   └── Phase 2 (download, Stage 0) ── Phase 3 (cardiomyocytes, gate)
                                          ├── Phase 4 (T cells, if GO)
                                          └── Phase 5 (reproducibility)
                                                         └── Phase 6 (synthesis) ◄─┘
```

- Phase 3a is a hard gate. If our pipeline cannot reproduce the protocol's cardiomyocyte result,
  nothing after it is interpretable.
- Phase 1 stands on its own: it is reported even if Phase 3 stops.

## 5. Predictions, written down before any output

These are recorded so the result cannot be read selectively afterwards. Each will be scored right or
wrong in the report.

1. **P4 is miscalibrated.** Cells, not donors, are its unit, and screens have thousands of cells, so a
   large share of random null genes will reach `Sig` = 1. The prediction is above 0.05 in at least two of
   LUAD, E2 and Chaffin. Our donor-level criteria (O7 to O11) will stay at or below 0.05.
2. **Our criteria are conservative.** They will pass GSN as a single pre-specified test, but fail most
   of the protocol's other Fig. 2b genes, mainly at the donor-level step (O7) or the both-operations step
   (O11).
3. **Agreement is low.** Fewer than half of the P4 `Sig` genes in our cohorts will reach a positive
   status under O, giving a kappa below 0.4.
4. **The classifier criteria agree** (P3 with O4) in every cohort.
5. **Transfer is poor under either set.** Lung-to-colon sign agreement of P4 `Sig` genes will not exceed
   chance by more than it does under O.

## 6. Compute (estimates; fixed at registration)

Costs are taken from E2, measured on thinkstation1:

- five fine-tunes took 1,266 to 1,318 s each;
- the perturbation screen took 14,022 gene × donor × operation calls on 100 cells each, from
  30 September 17:42 to 2 October 01:20 UTC (31.6 h). That is about 8 s per call, or about 0.08 s per
  perturbed cell.

| Phase | Work | Estimate |
|---|---|---|
| 1 | CPU re-analysis of existing outputs | hours of CPU, no GPU |
| 3a | 5 fine-tunes (plus tuning trials, if the protocol's Ray Tune search is used) and a deletion screen over all expressed genes in DCM cardiomyocytes, capped by `max_ncells` | fine-tunes about 2 GPU-h each pass. The all-gene screen is the largest item: at 0.08 s per perturbed cell it scales with cells × expressed genes per cell, likely tens to over a hundred GPU-h. It is measured on a 50-cell pilot before launch |
| 3b | about 250 genes (GSN, 20 named, controls, 100 null) × 9 DCM donors × 2 operations | about 4,500 calls, about 10 GPU-h |
| 4 | as 3a and 3b on T cells, scaled by the Stage 0 counts | fixed after Stage 0 |
| 5 | probes and one extra fine-tune set | under 10 GPU-h |

**Stop rule.** If a pilot projects a phase past its registered ceiling, the phase stops and the human
decides. Nothing is dropped to fit the budget: no gene, donor, control or direction, as in the LUAD
registration s.8.

## 7. Deliverables

1. `phase1_results/`: P4 statistics, the calibration and agreement tables, and a short report (first
   results; CPU only).
2. The Stage 0 counts and the frozen Chaffin design files.
3. `phase3_results/` and `phase4_results/`: both criteria sets on the same fine-tuned folds.
4. The scorecard, the IMRaD report and the graphical abstract.
5. If warranted, a proposed ISP-STD-1 amendment with the evidence for each change.
