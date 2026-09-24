# Balanced-donor Geneformer workstream — Phase 0 dataset selection rule

- Registered: 2026-09-24, before any candidate dataset was screened.
- Author: Phyllis (phyllis-mudz5ozf). Brief: god, 2026-09-24T05-45-00Z.
- What I had seen when I wrote this: the brief, and the lung repo's
  `current_workflow/METHODS.md` donor-split table (structural counts only). I had
  seen no dataset's outcomes, and no dataset's donor table except that one.
- This file is append-only. Any change goes in a dated AMENDMENT block at the
  end, and the original text stays readable.

## 0. What is being selected for

The workstream mirrors the lung T-cell work: a Geneformer disease-state
**classifier on T cells** (CD4 and CD8, with no subtype restriction), then
in-silico perturbation (ISP) of held-out T cells toward another state.

- **Contrast:** T cells from **tumour tissue** vs T cells from **normal tissue
  of the same organ**. That is exactly two groups. A third group is allowed
  only if it passes every criterion on its own.
- **Preferred organ:** lung. Another organ is allowed; see rule R3.

## 1. Definitions

- **Qualifying donor (for a group):** a donor who has at least **C_MIN = 100**
  T cells in that group, after QC (singlets, non-missing donor ID, raw counts
  available), and whose sample site matches the group's single site (F2).
- **Group donor count n_g:** the number of qualifying donors in group g.
- **The cohort is built balanced.** From each qualifying donor, draw exactly
  **C_MIN cells** for the ISP/analysis cohort. The draw is a fixed-seed random
  draw of cell IDs, recorded before tokenisation. Training may use up to
  **C_TRAIN = 3 x C_MIN = 300** cells per donor, with the same seed rule.

## 2. Criteria. Every one is structural: donor tables, cell counts, labels. No outcome is used.

| # | Criterion | Number | Reason |
|---|---|---|---|
| A1 | Qualifying donors per contrasted group | **n_g >= 12** hard floor; **>= 16 preferred** | See "Where I differ", point 1. Held-out evaluation and ISP use donors that were not trained on, so the donors that reach the statistics are a subset. 12 lets a donor-disjoint split keep >= 6 per group for analysis (clears the exact-test floor, min two-sided p 0.0022) and >= 6 for training. |
| A2 | Cells per qualifying donor | >= C_MIN = 100 T cells | Your (d) asked for >= 50. At 100, a per-donor mean over a subset of cells (for example one expression stratum) stays above 50. It also matches the S100 100-per-donor cap already in use. |
| A3 | Largest donor share of a group's cells (analysis cohort) | <= 25% | Your (b). Under the equal cap it holds automatically: the share is 1/n_g <= 8.3%. It stays as a backstop in case the construction breaks. |
| A4 | Max/min donor cells within a group | analysis cohort = 1.0 (equal cap); training cohort <= 3x | Your (c). Equal counts make donor weighting and cell weighting identical, which is the point of the workstream. |
| F1 | Cell labels | T cells labelled (CD4/CD8 or T cell), plus disease and tissue labels, in the source metadata | Your (e). |
| F2 | One sample site per group | Tumour group = primary tumour tissue only. Normal group = normal or adjacent-normal tissue of the same organ only. No effusions, blood, lymph nodes or metastases mixed into a group. | PleuralEffusion carried 74.9% of SCLC cells. That was a site mixture as much as a donor imbalance, and donor balance alone does not remove it. |
| F3 | No study-by-group confounding | Every contributing study supplies donors to **both** groups, **or** the groups are paired within donor. A group drawn from studies the other group never uses is disqualifying. | Otherwise a classifier learns study, chemistry and site batch, and balanced donors do not help. |
| F4 | One assay family | 10x 3' or 10x 5' across both groups. Mixed chemistries are allowed only if each chemistry appears in both groups. | Same reason as F3. |
| F5 | Raw integer counts available | Required | Geneformer tokenisation needs raw counts (METHODS s.3). |
| F6 | Access | Public, with no data-use application needed to download | Phase 2 must not stall on an access request. |

## 3. Ranking among the candidates that pass (applied mechanically, in order)

- **R1.** Paired design first: tumour and normal from the same donors, with at
  least 8 such donors qualifying in both groups. A signed test over d donors has
  a minimum two-sided p of 2/2^d, which is 0.0078 at d = 8. Pairing also removes
  between-donor variation from the contrast.
- **R2.** Larger min(n_g) across the groups.
- **R3.** Lung before any other organ (continuity with the existing classifier
  and the reusable code).
- **R4.** Smaller download size.

## 4. Stop conditions

- If no candidate passes A1-A4 and F1-F6, I report that and stop. Thresholds are
  not relaxed to rescue a candidate. A relaxation needs a dated amendment
  approved by god **before** re-screening.
- Every candidate I screen is recorded with its numbers, pass or fail, and the
  first criterion it failed.

## 5. What the screen may and may not look at

- **May:** donor IDs, per-donor cell counts by group, cell-type, disease,
  tissue, assay and study labels, file sizes.
- **May not:** expression values, differential expression, marker genes,
  classifier accuracy, embeddings, or any published result about these donors'
  T-cell states. This holds until the dataset is chosen and god accepts the
  choice.

## Where I differ from the brief's defaults

1. **(a) >= 6 per group is necessary but not sufficient.** It has to hold for
   the donors that enter the statistics, and those are held-out donors. With a
   donor-disjoint split, >= 6 in the test partition needs roughly >= 12 in
   total. The alternative is k-fold donor cross-fitting, where every donor is
   scored by a model that never saw it. That lets all n_g donors enter the
   statistics, but it costs k fine-tunes of GPU time, which is a Phase 4
   decision.
2. **(b) and (c) are properties of the cohort I build, not of the raw
   dataset.** Any dataset with enough donors can be balanced by capping cells
   per donor. What cannot be fixed afterwards are donor count (A1), cells per
   donor (A2), and site and study confounding (F2-F4). F2-F4 are my addition.
   They are the part of the PleuralEffusion failure that balancing does not
   touch.
3. **Known cost of the equal cap:** excluding donors with fewer than 100 T cells
   selects for T-cell-rich samples. I record it here as a limitation. It is not
   corrected for.

## AMENDMENT 1 — 2026-09-24T08:56:21Z — panel-size floor (requested by god, 2026-09-24T08-55-00Z)

- **Timing, stated plainly:** this was added **after** the Phase 1 structural
  screen finished (PHASE1_SURVEY.md, sent before god's request reached me) and
  **before** any expression value, gene-eligibility count or panel exists. No
  structural criterion above changes.
- **P1. The primary statistic is donor-level.** With the paired design (R1),
  the Phase 5 primary is registered as a test across **held-out paired
  donors** (one within-donor value per donor). Its n is donors, which is what
  A1 buys power for.
- **P2. Held-out floor:** >= 6 paired donors in the analysis partition, or
  all donors under k-fold donor cross-fitting. This cannot be traded for
  training donors. If a split would leave fewer than 6, that split is not
  used.
- **P3. Gene-panel floor for any gene-level statistic** (primary or
  secondary): **>= 12 eligible genes**. For a two-group split of genes, that
  means >= 6 per side, where the minimum two-sided p is 2/C(12,6) = 0.002165.
  For a rank correlation across genes, n >= 12 as well. God's table:
  8 genes gives 0.0286 (perfect split only); 12 gives 0.0022.
- **P4. Order:** choose the dataset on structural grounds, then run gene
  eligibility on the chosen data. If fewer than 12 genes are eligible, the
  gene-level analysis is **not run**. It is reported as "not run: n_eligible
  = k < 12", never as "no effect". There is no substitute gene and no lowered
  floor without a new dated amendment approved by the human.

## AMENDMENT 1a — 2026-09-24T09:03:29Z — dependency note (requested by god, 2026-09-24T09-21-00Z)

- **One part of Amendment 1 did depend on the screen.** P1 (donor-level primary)
  was chosen **because the Phase 1 screen found a paired design** (43 donors
  with both tumour and adjacent-normal samples).
- That dependency is structural: the test is matched to the design. It was not
  chosen to maximise significance, and no outcome existed when it was chosen.
- P2-P4 do not depend on the screen. P3 is arithmetic, P2 follows god's
  held-out ruling, and P4 is an ordering rule.
