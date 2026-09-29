# Phase 4 costing: single donor-disjoint split vs k-fold donor cross-fitting (2026-09-24)

This is a costing for god to take to the human. Nothing here has run on a GPU.

## The one measured rate

- **Measurement:** Geneformer-V2-**316M**, bf16, 1 epoch, batch 8, on the
  NVIDIA GB10: **7,242.77 s for 3,401 steps**, which is about 27,208 training
  cells (`sclc_validation/bf16_bench/RESULTS_BF16_316M.md`, 2026-09-18).
- **Rate:** about **3.76 training cells/s**, or 0.47 steps/s.
- **No measurement exists for the 104M model** that the July workflow used.
  - The July 104M run had 1,786 steps × 8 ≈ 14,288 cells.
  - Its start time is lost: the August migration rewrote file mtimes, and the
    trainer state records no runtime.
- **How the rate is used below:** as an **upper bound** for 104M. The model is
  about a third of the size, so it is *inferred* to be faster. That is not
  measured.
- **Uncertainty not covered:** token length. Our median is 755 tokens per
  cell. The SCLC benchmark's median is not recorded here, and attention cost
  grows with length.

## Cells per fine-tune (from the drawn cohort: <= 300 cells per donor per side)

The training pool has 25,700 cells over 43 donors, about 598 per donor across
both sides.

| Design | Train donors | Train cells | Eval donors | Held-out donors reaching ISP stats |
|---|---|---|---|---|
| A. single split | 22 | ~13,150 | 6 | **15** |
| B. k = 5 cross-fitting | ~30 per fold (43 - ~9 held out - ~4 eval) | ~17,950 per fold | ~4 per fold | **43** (every donor once) |
| B'. k = 3 cross-fitting | ~25 per fold | ~14,950 per fold | ~4 per fold | **43** |

## Fine-tune GPU hours

The "upper bound" column uses the 316M rate. The "inferred" column assumes
104M runs about 3× faster, which is **not measured**.

| Design | Fine-tunes | Upper bound (316M rate) | Inferred for 104M |
|---|---|---|---|
| A. single split | 1 | **~1.0 h** (+ ~0.1 h held-out eval) | ~0.35 h |
| B. k = 5 | 5 | **~6.6 h** (+ ~0.3 h eval) | ~2.2 h |
| B'. k = 3 | 3 | **~3.3 h** (+ ~0.2 h eval) | ~1.1 h |

## What k-fold actually buys (precision, not significance)

- **The p-value floor does not bind in either design.** A has a floor of
  2/2^15 = 6.1e-5; B has 2/2^43. Neither is a reason to choose k-fold.
- **Effect-size precision.** The standard error of a mean within-donor shift
  scales as 1/√n_donors. At 43 vs 15 donors, the interval is
  √(43/15) = **1.69× narrower, about 41% narrower**.
- **Generalisation.** The estimate covers all 5 studies in their real
  proportions, instead of whichever donors landed in one test split.
  - Under A, Leader_Merad (23 of 43 donors) could dominate or vanish from a
    15-donor test set depending on the split.
  - A could be stratified by study instead; I recommend that if A is chosen.

## What k-fold costs beyond the fine-tunes

1. **Phase 6 ISP cost also scales, by about 43/15 = 2.9×.** Every donor's
   analysis cells must be perturbed with its own fold's model, so 43 donors
   are perturbed instead of 15. The ISP cost is per perturbed cell, so this
   multiplier probably matters more than the fine-tune multiplier. Phase 6
   will be costed separately.
2. **Between-model variance.** Donors are scored by k different models, so
   fold-to-fold model differences enter the donor-level values.
   - Mitigation: identical recipe and seed per fold, and fold recorded as a
     column. Phase 5 can register a fold-stratified check.

## Recommendation

1. **First, a measured calibration.** Before either number goes to the
   human, run one **104M fine-tune on ~1,000 cells, 1 epoch**, in the
   thinkstation1 environment, to measure the real 104M rate at our token
   lengths.
   - Cost: well under **0.1 GPU-h**.
   - It is still GPU time, so it needs the human's approval like everything
     else.
   - It replaces the "inferred" column with a measurement.
2. **Then decide A vs B on precision.** If the human accepts a roughly 3×
   Phase 6 cost, I recommend **B with k = 5**: all 43 donors, and a 41%
   narrower interval.
   - If Phase 6 cost dominates, **A with a study-stratified split** is the
     honest cheaper option.
   - B' (k = 3) buys the same all-donor coverage for less fine-tune time.
     Its larger held-out folds leave each model about 25 training donors.

## Fixed recipe (to be registered before any fine-tune, so nothing is tuned on held-out donors)

- **Model:** Geneformer-V2-104M, matching the July workflow.
- **Training settings:** 1 epoch, lr 5e-5, batch 8, freeze_layers 6,
  seed 43, fp32 as in July. bf16 is an alternative only if registered in
  advance.
- **Label:** `origin` (tumor_primary vs normal_adjacent).
- **Split:** donor-disjoint, so both tissues of a donor are always in the same
  partition.
- **Selection:** no hyperparameter search, and the eval split is not used for
  model selection.

---

## v2 — 2026-09-24: one table for the human (fine-tune + ISP, after god's Phase 5 rulings)

This version adds Panel B (36 genes) and overexpression, and it costs the
perturbation (ISP) stage for the first time.

### The only measured ISP rate

- **Measurement:** 316M, bf16. **6,099 s for 150 units**, where a unit is one
  gene on one source set of up to 300 cells (overexpression). That is
  **40.7 s per unit**.
- **No measurement exists for 104M or for fp32.** For 316M, bf16 was 3.44×
  faster than fp32. So 104M-fp32 could land near the 316M-bf16 rate. That is
  *inferred*, not known.

### Genes and cells per run

- **Genes perturbed:** Panel A 15 + Panel B 36 + about 120 matched controls
  (six shared strata of 20) = **171**, in both directions.
- **Cells:** 100 tumour analysis cells per held-out donor, counted as
  300-cell units.
  - Single split: 15 donors = 1,500 cells = 5 units per gene-direction.
  - k-fold: 43 donors = 4,300 cells = 14.3 units per gene-direction.

### Upper bounds

These assume every analysis cell is perturbed for every gene. Deletion only
touches cells that express the gene, so low-detection genes (Panel A,
roughly half the controls) will cost less. Measuring that is what the
calibration run is for.

| Configuration | Fine-tune | ISP | **Total GPU-h (upper bound)** |
|---|---|---|---|
| **Calibration (ask first)**: 104M fine-tune on ~1,000 cells, plus 1 gene × 300 cells in both directions | ~0.08 | ~0.03 | **<= 0.15** |
| Recommended if cost allows: k = 5 cross-fit, Panels A+B, del+ovx | ~6.6 | ~55.4 | **~62** |
| **Fallback 1 (my recommendation on current numbers):** single study-stratified split, Panels A+B, del+ovx | ~1.0 | ~19.3 | **~20** |
| Fallback 2: single split, Panel A only, del+ovx | ~1.0 | ~15.2 | **~16** |

### What the table shows that the earlier framing did not

- **Controls dominate the ISP cost.** About 120 of the 171 genes are
  controls. So cutting Panel B, god's first cut, saves only about 21%
  (20 h → 16 h), because the control strata remain.
- **The cut that saves real time is k-fold vs single split** (about 3×).
- **Overexpression and controls are non-negotiable** for the claims made. The
  concordance check needs overexpression, and without matched controls there
  is no control-adjusted value.
- **Never cut:** controls, overexpression, or the donor floor.
- **What k-fold buys:** precision (1.69× narrower interval, all 5 studies).
  Its costs: about 3× Phase 6, plus between-model variance.

### Order of asks to the human

1. The <= 0.15 GPU-h calibration. It replaces every inferred rate here with a
   measured one.
2. With the measured rates in hand, Fallback 1 (about 20 GPU-h upper bound),
   or k-fold if the measured cost makes it affordable.

The budget is a hard stop: on a projected overrun, stop and report. Nothing
is dropped to fit.

---

## v3 — 2026-09-24: MEASURED rates (calibration run, human-approved), Geneformer-V2-316M bf16

The human directed that all work use Geneformer-V2-316M in bf16 (registration
Amendment 1). The approved calibration was run on the thinkstation1 GB10 with
that exact stack. Record: `provenance/calibration_rates_316m_bf16.json`.

### Calibration run

| What | Measured |
|---|---|
| Fine-tune, 1 epoch, 900 cells, median 745 tokens | 109.95 s, so **8.19 training cells/s**. The v1 upper bound assumed 3.76 from the 316M benchmark. |
| ISP deletion, 1 gene × 300 cells | **20.4 s** (the 316M benchmark gave 40.7 s per unit) |
| ISP overexpression, 1 gene × 300 cells | **18.8 s** |
| GPU time used by the calibration | 153.8 s successful run, plus about 16 s across three failed starts, so about **0.05 GPU-h of the 0.15 approved** |

The three failed starts were:
1. a prepare_data file name;
2. the base-model path, where the `/srv/lab` 316M file is an LFS pointer;
3. the corrected run's first attempt.

Each failed start stopped before any training step.

### Costing from measured rates

- The ISP figures use 20.4 s per gene-direction per 300 cells, which is
  conservative because each timed call also included a model load.
- They are still upper bounds: every analysis cell is assumed to be perturbed
  for every gene. Deletion only touches cells that express the gene, so real
  cost is lower.
- Embeddings (per-donor goal centroids and the held-out evaluation) are about
  0.5 h.

| Configuration | Fine-tune | ISP | **Total GPU-h (upper bound)** |
|---|---|---|---|
| **Recommended: single study-stratified split, Panels A+B, del+ovx** | 0.45 | 9.7 | **~10.6** |
| k = 5 cross-fit, Panels A+B, del+ovx | 3.05 | 27.8 | **~31** |
| Single split, Panel A only, del+ovx | 0.45 | 7.7 | **~8.6** |

**Never cut:** controls, overexpression, the donor floor. The budget is a hard
stop. On a projected overrun, stop and report. Nothing is dropped to fit.

---

## v4 — 2026-09-24: conditional contingency (registration Amendment 2)

| Contingency | When | GPU-h (upper bound, measured 316M rates) |
|---|---|---|
| 104M Panel A control arm | Only if Panel A on 316M has at least one eligible non-replication (OPEN / DELETION_ONLY / DOSE_INCOHERENT / REVERSED). Needs its own approval when triggered. | ~8.6 (fine-tune 0.45 + ISP 7.7 + embeddings ~0.5). 104M is smaller, so this is expected to be lower. |

This is not part of the ~10.6 GPU-h recommended run. It is decided by the
Panel A outcome, under the rule registered in advance, not by whoever is
looking at the result.
