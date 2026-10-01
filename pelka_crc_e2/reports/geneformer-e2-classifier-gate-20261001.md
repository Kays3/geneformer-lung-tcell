# E2 interim report: the colorectal tumour/normal T-cell classifier gate (H2a)

**Study:** E2, the balanced-donor Geneformer design re-run on Pelka et al. 2021 colorectal cancer (GSE178341), 19 donors, unsorted cells only.  
**Status:** interim. The classifier gate (H2a) is complete; the perturbation screen (H2b, H2c) is still running and is not reported here.  
**Date:** 1 October 2026, written 11:10 JST (02:10 UTC); revised after review 11:13 JST (02:13 UTC).  
**Sources (read only):** commit `2c104ab` on branch `analysis/e2-pelka-crc-20261001`, files `pelka_crc_e2/phase4_results/classifier_gate.json`, `noop_gate.json` and `goal_embeddings_record.json`; `pelka_crc_e2/registration/E2_REGISTRATION.md` (passed by Stanley at `8b11d5d`). LUAD context: `balanced_donor_luad/phase4_results/classifier_gate.json` on `origin/main`.

## Summary

We fine-tuned Geneformer-V2-316M, once per cross-validation fold, to tell tumour-infiltrating T cells from T cells of the same donor's normal colon. Each donor was scored only by a model that never saw any of that donor's cells. Pooled held-out balanced accuracy was 0.904 against a registered bar of 0.60. All 19 donors scored above chance (range 0.655 to 0.970), and the exact two-sided sign test gave p = 1/262,144 (3.8 × 10⁻⁶). The gate passed with margin, and no plausible run-to-run variation would change that reading. The result shows that tumour and normal-colon T cells of one donor are separable after fine-tuning on other donors from the same study. It does not show that the LUAD classifier transfers to colon (that question was E1), and it says nothing yet about the perturbation results.

## Introduction

The balanced-donor LUAD study found that a fine-tuned 316M classifier separates tumour from adjacent-normal T cells of the same donor: pooled held-out balanced accuracy 0.825, with 43 of 43 donors above chance. E2 asks whether that result and the two perturbation results built on it hold outside lung. H2a is the first and gating question. If a classifier trained by the same recipe could not separate colorectal tumour from normal-colon T cells, the perturbation screen would have no tumour-to-normal axis to measure along, so a FAIL was registered to stop E2 before any perturbation. This report covers H2a only, because the human asked for it before the screen finishes.

## Methods

**Cohort.** The cells are CD4 and CD8 T cells (author annotation `clMidwayPr` TCD4 or TCD8), restricted to unsorted samples in both tissues, from primary tumour (T) and normal colon (N). Every donor with at least 100 such cells in both tissues was kept, giving 19 donors (9 MMR-deficient (MMRd), 10 MMR-proficient (MMRp)). Restricting to unsorted cells removes the tumour-only CD45+ sorting imbalance found in 11 of the 25 donors when all processing types were allowed. A fixed-seed draw, unchanged from the LUAD study, gives exactly 100 analysis cells per donor per tissue (3,800 cells) and up to 300 training-pool cells per donor per tissue (9,039 cells). Chemistry (10x 3' v2 or v3) never differs between a donor's two tissues. Raw tumour counts per donor ranged from 137 to 2,220 cells and normal counts from 104 to 1,457. Because of the 100-cell cap, each donor carries exactly 1/19 of the pooled score.

**Model.** The model is fine-tuned, not zero-shot, and the LUAD fold models are not reused. For each of five folds, a new classifier was fine-tuned from the base Geneformer-V2-316M weights (sha256 `965cecce…`) with the LUAD recipe unchanged:

- 1 epoch, learning rate 5e-5, batch size 8;
- the first 6 layers frozen, seed 43;
- label `origin` (tumour vs normal), bf16 precision, with the registered tie-handling fix.

Training used the pool cells of the fold's training donors.

**Folds.** Five-fold donor cross-fitting (seed 20260924) kept both tissues of a donor in the same partition. The folds held out 4, 4, 4, 4 and 3 donors. In each fold, 2 of the remaining donors were set aside as an evaluation set that was never used for model selection, leaving 13 training donors (14 in fold 4) (LUAD set aside 4 and trained on about 30). Each held-out donor was scored on its own 200 analysis cells by its fold's model.

**Gate (registration s.6.1, unchanged from LUAD).** PASS required both of the following:

- pooled held-out balanced accuracy ≥ 0.60;
- an exact two-sided sign test over donors of "per-donor balanced accuracy > 0.5", with p ≤ 0.05 and more donors above chance than below.

With 19 donors the sign test needs at least 15 above chance (p = 0.019; 14 gives 0.064). We also report per-fold values, tie counts, and a non-determinism band. The band is the largest per-donor change seen between two same-seed runs of LUAD fold 0, 0.045. It was measured on LUAD, is reported as is, and was never allowed to change PASS.

**No-op gate.** Before any perturbation, the pipeline ran two forward passes with no edit on the first donor (C107), with B2M as the nominal gene. Every per-cell shift had to be exactly 0.

**Compute.** All runs were on thinkstation1 (GB10) in the frozen environment recorded in the registration. That environment differs from the LUAD one only in pandas (2.3.3 instead of 3.0.5). The five fine-tunes ran from 00:41 to 02:29 JST on 1 Oct (15:41 to 17:29 UTC, 30 Sep), taking 1,266 to 1,318 s each. Goal embeddings took a further 637 s.

## Results

**H2a: PASS.** Pooled held-out balanced accuracy was 0.904 over 3,800 cells. All 19 donors scored above chance and none was at or below 0.5, giving a sign-test p of 1/262,144 = 3.8 × 10⁻⁶. No predictions were tied in any fold (0 of 800, or 600 in fold 4, test cells).

| Fold | Held-out donors | Balanced accuracy |
|---|---|---|
| 0 | C111, C125, C137, C140 | 0.911 |
| 1 | C110, C129, C132, C162 | 0.889 |
| 2 | C107, C130, C134, C170 | 0.866 |
| 3 | C123, C142, C155, C157 | 0.935 |
| 4 | C126, C135, C143 | 0.925 |

No fold fell below chance.

| Donor | MMR | BA | Donor | MMR | BA |
|---|---|---|---|---|---|
| C107 | MMRp | 0.945 | C137 | MMRd | 0.955 |
| C110 | MMRd | 0.880 | C140 | MMRp | 0.855 |
| C111 | MMRd | 0.865 | C142 | MMRd | 0.920 |
| C123 | MMRd | 0.965 | C143 | MMRd | 0.950 |
| C125 | MMRp | 0.970 | C155 | MMRp | 0.960 |
| C126 | MMRp | 0.950 | C157 | MMRp | 0.895 |
| C129 | MMRp | 0.825 | C162 | MMRp | 0.885 |
| C130 | MMRd | 0.925 | C170 | MMRd | 0.940 |
| C132 | MMRd | 0.965 | | | |
| C134 | MMRp | 0.655 | | | |
| C135 | MMRp | 0.875 | | | |

![Held-out balanced accuracy per donor, ordered, coloured by fold, with circles for MMRd and squares for MMRp donors. The dashed line is the 0.60 gate, the dotted line chance and the orange line the pooled value. Right: the 43 LUAD per-donor values beside the 19 E2 values, for context only (different cohort, study and fold models).](geneformer-e2-classifier-gate-20261001-fig1.png)

**Robustness of the reading.** The pooled value sits 0.304 above the bar, about seven times the 0.045 non-determinism band. No donor lies within the band of 0.5, so the worst-case sign-test p equals the observed one, and the gate file records the result as clean. C134 is the lowest donor at 0.655. It also had the fewest raw tumour T cells (137), so its 100 tumour analysis cells were drawn from a small sample. Even C134 stays 0.155 above chance, more than three times the band.

**Descriptive, not tested.** The median per-donor balanced accuracy was 0.940 in MMRd donors (n = 9) and 0.890 in MMRp donors (n = 10). The three lowest donors (C134, C129, C140) are all MMRp. With these numbers we draw no conclusion about mismatch-repair status, and no test was registered.

**No-op gate: PASS.** In all 100 tumour cells, 100 own-normal cells and 100 cells scored against the fold's global-normal goal, the largest absolute shift was 0.0, with no non-zero values. The goal embeddings were built and verified for all 19 donors with all 5 fold models.

**Context: LUAD.** On the same gate, LUAD gave a pooled balanced accuracy of 0.825, per-fold values of 0.784 to 0.869 and per-donor values of 0.55 to 0.985 (median 0.85), with 43 of 43 donors above chance. E2 is higher and tighter (per-donor median 0.925). We do not read this difference as colon T cells being more separable than lung T cells, for the reasons set out in the Discussion. It is given only to place the E2 value.

## Discussion

The gate passed by a wide margin, and every donor contributed to the pass rather than a few strong donors carrying a weak majority. The registered reading is that the tumour/normal T-cell signal is learnable in colon under the same recipe. The perturbation screen can therefore proceed on an axis the fold models do separate.

The claim is narrow. Each fold model learned to separate the two tissues from 13 other donors of the same study, processed in the same way, and generalised to held-out donors of that study. That is separability within one colorectal study. It is not evidence that the LUAD models, or the tumour/normal axis they learned, carry over to colon. That transfer question was E1, and E2 trains its own models from the base weights by design.

The comparison with LUAD (0.904 vs 0.825) should not be read biologically. The two cohorts differ in study (one study against five), dissociation and sorting history, T-cell annotation, 10x chemistry mix, donor count (19 against 43), the amount of training data per fold, and the fine-tuned models themselves. Any of these could account for a gap of this size. A single study may also simply be more homogeneous between its own donors, which makes held-out donors easier to score.

What the classifier reads remains open. A high balanced accuracy is consistent with T-cell-intrinsic state differences between tumour and normal colon. It is equally consistent with differences in ambient RNA from the surrounding tissue, which differs between tumour and normal samples within a donor. As in LUAD, no ambient classifier was fitted for colon, so the gate cannot tell these apart. The CD4:CD8 balance also differs between tissues (training pool: tumour 3,215 : 1,893, normal 2,818 : 1,113), and part of the separation may reflect subset composition rather than state within a subset.

## Limitations

- Nineteen donors from one study. The per-donor values are well estimated from 200 cells each, but the cohort-level picture rests on a single centre's processing.
- The non-determinism band comes from one same-seed repeat of one LUAD fold. It is an observed magnitude, not a bound, and was not re-measured on colon. Given the margins above, this does not affect the reading.
- The evaluation set is 2 donors per fold instead of the LUAD 4, as registered, because 19 donors leave fewer to spare. It was not used for selection, so it does not bias the held-out score.
- The training-pool cap is 300 per tissue, but 21 of the 38 donor-tissue pairs had fewer cells (minimum 104), so donors do not contribute equal amounts of training data. Held-out scoring is unaffected because analysis cells are capped at exactly 100.
- This report covers H2a only. Nothing in it predicts the outcome of H2b (the null anti-correlation, primary) or H2c (the Panel B pattern).

## Status of the remaining steps

The perturbation screen started at 02:42 JST on 1 Oct (17:42 UTC, 30 Sep). It covers 369 genes, deletion and overexpression, in all 19 donors. At the 04:44 JST check (19:44 UTC) it was running at about 280 GPU-s per gene, about 2.3 times the per-gene cost behind the registered 12.3 GPU-hour prediction. The projected total is about 31 GPU-hours, within the 52-hour authorisation and the driver's hard ceiling. The expected finish is about 07:30 JST on 2 Oct (22:30 UTC, 1 Oct). Writing this report used CPU only and did not touch the running job. No deviation from the registration has occurred. The ISP runs about 2.3 times slower than predicted, which stays within the registered ceiling; it will be noted in the dated log.
