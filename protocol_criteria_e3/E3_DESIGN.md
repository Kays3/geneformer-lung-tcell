# E3 design: the Geneformer protocol's evaluation criteria in T cells of the protocol's own atlas

**Status: DRAFT DESIGN, not registered. The full cross-evaluation it belongs to is planned in
[`EXECUTION_PLAN.md`](EXECUTION_PLAN.md), with the criteria listed in
[`CRITERIA_INVENTORY.md`](CRITERIA_INVENTORY.md). Nothing has been downloaded, tokenised, fine-tuned or
perturbed for E3.** Three things must happen before any compute:

1. the human's GO for the data download (Stage 0, s.3) and for the GPU budget (s.7);
2. a gatekeeper pass on this design;
3. a registration commit that freezes it.

After registration, nothing in this file may change except a dated deviation log appended at the end.

- Written 2026-10-03 on branch `claude/te-work-analysis-9mg3lz` of `Kays3/geneformer-lung-tcell`.
- Reference protocol: Zhang Y, Venkatesh MS, Theodoris CV. Discovery of candidate therapeutic targets
  with Geneformer. *Nature Protocols* (2026), https://doi.org/10.1038/s41596-026-01364-8. Below, "the
  protocol".
- Reference atlas: Chaffin M, et al. Single-nucleus profiling of human dilated and hypertrophic
  cardiomyopathy. *Nature* 608, 174–180 (2022), the protocol's example dataset (its reference 18).

## 1. Why this study

The protocol's published evaluation of a fine-tuned model and its in silico perturbation (ISP) rests on
four criteria. Their worked example on the Chaffin atlas reports the value given for each.

| # | Protocol criterion | Where in the protocol | Reported on the Chaffin atlas |
|---|---|---|---|
| P1 | Disease classification on held-out patients (split by `individual`, balanced on disease, age and sex), reported as a confusion matrix and macro F1 | Steps 19–21 | Single task, cardiomyocytes, nonfailing (NF, n = 9), hypertrophic (HCM, n = 11) and dilated (DCM, n = 9) cardiomyopathy: macro F1 85% |
| P2 | Per-gene shift of the start state toward the goal state (`Shift_to_goal_end`), compared with perturbing random genes by Wilcoxon test, with BH correction; `Sig` = 1 if `Goal_end_FDR` < 0.05 | Step 26, `InSilicoPerturberStats(mode="goal_state_shift")` | DCM cardiomyocytes, goal NF, alternative HCM: genes that shift toward NF and away from HCM (Fig. 2b) |
| P3 | Quantized and full-precision perturbation shifts agree | Steps 33–34 | Pearson r > 0.99 |
| P4 | Experimental validation of a top-ranked target | Anticipated results | GSN knockout improved contractility in an iPSC model of DCM (from Theodoris et al. 2023) |

This repository's T-cell studies apply a stricter standard (ISP-STD-1). It uses donor-level values,
matched control genes, both deletion and overexpression, and a random-gene null for the two directions.

In lung (43 donors) and colon (19 donors), the classifier criterion held. Its counterpart is the
balanced-donor gate; pooled held-out balanced accuracy was 0.825 and 0.904, every donor above chance.

The random-gene comparison did not hold up. Random genes alone give opposed deletion and
overexpression shifts (rho −0.593 in lung, −0.245 in colon). So a gene with opposite effects in the two
arms is the baseline, and a shift that beats random genes on pooled cells is not by itself evidence of
specific regulation. The curated lung hits also did not replicate in colon: 3 of 10 kept their sign.

Both cohorts are cancer T cells, a setting the protocol itself flags. Its Limitations say malignant cells
were excluded from pretraining, so a cancer-tuned model may do better.

E3 removes that objection. It applies the protocol's criteria, exactly as published, to cells of the
protocol's own atlas. It first checks that our pipeline reproduces the protocol's cardiomyocyte result.
Then, in the T cells of the same hearts, it reports P1 to P3 next to the ISP-STD-1 readings, gene by gene.

## 2. Questions and hypotheses

- **H3a (pipeline check, gating).** On Chaffin cardiomyocytes, the protocol's single-task recipe
  reproduces its reported performance on held-out patients. The reading is `REPRODUCED` if macro F1 is
  at least 0.80 (protocol: 0.85). The perturbation run on DCM cardiomyocytes, with the protocol's multi-task
  model (its source of Fig. 2b), must also place GSN among genes with `Sig` = 1 that shift toward NF. If either fails, E3 stops and is reported as a pipeline mismatch,
  because nothing downstream could then be blamed on T cells.
- **H3b (P1 in T cells).** A classifier fine-tuned on the atlas's T cells (or lymphocytes, s.3)
  separates NF, HCM and DCM on held-out patients. It is reported with the protocol's metrics (confusion
  matrix, macro F1) and with this repository's donor-level gate: per-donor balanced accuracy and an exact
  sign test over held-out donors.
- **H3c (P2 against ISP-STD-1), primary.** Among T-cell genes that the protocol's
  `goal_state_shift` marks `Sig` = 1 toward NF (deletion, DCM start state), the share that also reaches a
  positive ISP-STD-1 status, by the donor-level, control-adjusted, two-operation test, is below one half.
  The share and its exact binomial 95% CI are reported whatever they are.
- **H3d (the built-in null).** Across random eligible T-cell genes, deletion and overexpression shifts
  toward NF are negatively rank-correlated. This is the lung and colon R3 test, run unchanged.
- **H3e (P3, descriptive).** Quantized and full-precision shifts agree (Pearson r) on the T-cell P2 run.

## 3. Stage 0: data and feasibility (CPU only, outcome-free)

- **Source.** The Chaffin et al. 2022 single-nucleus atlas. The deposit and its accession are recorded
  with sha256 in the registration; the protocol's Colab tutorial copy is also checked. Neither is
  reachable from the session that wrote this draft (network policy), so nothing has been read yet.
- **Ask first.** The download size and destination are put to the human before any byte is fetched.
- **Cells.** T cells by the authors' annotation. If the atlas annotates only a lymphocyte class, that
  class is used and every result says "lymphocytes". NK cells are excluded if annotated separately.
- **Feasibility rule (registered before any count is read):**
  - **GO:** at least 6 donors per class (NF, HCM, DCM), each with at least 50 nuclei.
  - **GO-REDUCED:** fewer, but at least 6 donors with at least 50 nuclei in each of NF and DCM. E3 then
    runs as NF against DCM, and HCM is dropped as the alternative state.
  - **NO-GO:** anything less. E3 stops and reports the counts.
- **Cap.** A fixed-seed draw of at most 100 analysis nuclei per donor, so every donor carries equal
  weight. The rest form the training pool, at most 300 per donor.

## 4. Design

- **Model.** Geneformer-V2-316M, base weights pinned by sha256, bf16 with the registered tie-handling
  fix. Fine-tuning follows the protocol's Step 18–19 recipe. Hyperparameters come from its Ray Tune
  search run once on the training folds; the chosen values are frozen before any ISP.
- **Splits.** Donor-disjoint k-fold cross-fitting (k = 5, seed fixed at registration). Splits are by
  `individual`, balanced on disease, age and sex where the metadata allow. Every donor is scored only by
  a model that never saw it, and eval donors are never used for selection.
- **Arm P (protocol, as published).**
  - `InSilicoPerturber` with `perturb_type` delete and overexpress, `genes_to_perturb="all"`, start
    DCM, goal NF, alternative HCM, `emb_layer=0`, `emb_mode="cls"`.
  - `InSilicoPerturberStats(mode="goal_state_shift")` writes `Shift_to_goal_end`,
    `Goal_end_vs_random_pval`, `Goal_end_FDR` and `Sig`.
  - Run with full precision and once with `model_type="CellClassifier"` quantized for P3.
- **Arm S (ISP-STD-1).** The same fine-tuned fold models and the same DCM analysis nuclei.
  - Per-donor mean shift from the raw per-cell output, never from `Shift_to_goal_end`.
  - Adjusted against 20 matched control genes, matched on detection and rank.
  - Exact two-sided Wilcoxon over donors, Holm within family, and dose concordance between deletion
    and overexpression, all as in the balanced-donor LUAD registration.
  - The random-gene null (100 genes, fresh seed) for H3d.
- **Gene families.** Arm S cannot run on all genes; its cost is per gene and donor. It covers:
  - every gene with `Sig` = 1 in Arm P (deletion, toward NF);
  - an equal number of Arm P non-significant genes matched on detection;
  - the 100 null genes.
  If the `Sig` = 1 list exceeds the GPU budget (s.7), a fixed-seed sample is drawn from it, with the size
  fixed at registration.

## 5. Readings (fixed before any output)

| Hypothesis | Reading if met | Reading if not |
|---|---|---|
| H3a | `REPRODUCED` | `PIPELINE_MISMATCH` (stop) |
| H3b | `CLASSIFIER_PASS` (gate as in E2) | `CLASSIFIER_FAIL` (stop before ISP) |
| H3c | `PROTOCOL_SIG_NOT_SUFFICIENT` (share < 0.5, upper CI < 0.5) | `PROTOCOL_SIG_SUPPORTED` (lower CI ≥ 0.5) or `UNDECIDED` (CI spans 0.5) |
| H3d | `positive` / `control_draw_sensitive_open` / `negative`, as in E2 | as in E2 |
| H3e | descriptive r with 95% CI | descriptive |

## 6. Limits stated in advance

- Disease is a between-donor label here, unlike the within-donor tumour and normal tissues of E1 and E2.
  Any donor effect is therefore confounded with disease, and the sign test over held-out donors is the
  guard.
- Single-nucleus T-cell profiles are sparse, and cardiac T cells are few. A NO-GO in Stage 0 is a real
  outcome and will be reported as such.
- The protocol's P4 (experimental validation) cannot be tested computationally. E3 reports whether the
  protocol's validated target (GSN) is reproduced in cardiomyocytes, nothing more.

## 7. Compute and budget (to be confirmed with the human)

- **Hosts.** thinkstation1 and thinkstation2 (GB10), the hosts of E1 and E2.
- **Estimate.** Per-gene ISP cost in E2 ran about 2.5 times its registered prediction (E2 deviation D1),
  so this estimate uses the E2 measurement directly. Arm P with all genes on DCM nuclei is the largest
  item; the protocol puts a full screen at 12 h to 2 weeks. The registration states a GPU-hour ceiling
  and a stop rule once Stage 0 has fixed the cell and gene counts.
- **No new data** beyond the Chaffin atlas.
