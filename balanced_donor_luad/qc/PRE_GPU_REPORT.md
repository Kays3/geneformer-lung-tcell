# Pre-GPU report (2026-09-24): everything knowable before any GPU is used

The fine-tunes have **not** started. Every number below comes from CPU-only
work on committed inputs.

## 1. Fold assignment

`cohort/folds.json`: study-stratified 5-fold donor cross-fitting, seed
20260924.

| Fold | Held-out donors | Eval donors | Train donors | Held-out by study |
|---|---|---|---|---|
| 0 | 9 | 4 | 30 | Leader 5, Kim 2, He 1, Lambrechts 1 |
| 1 | 9 | 4 | 30 | Leader 5, Kim 2, He 1, Lambrechts 1 |
| 2 | 9 | 4 | 30 | Leader 5, Kim 2, He 1, Lambrechts 1 |
| 3 | 8 | 4 | 31 | Leader 4, Kim 2, He 1, Laughney 1 |
| 4 | 8 | 4 | 31 | Leader 4, Kim 2, He 1, Laughney 1 |

Every donor is held out exactly once, and both tissues of a donor are always
in the same partition. The tests fail on a broken split.

## 2. Ambient validity (the line that caps Panel A)

- **Leave-one-anchor-out AUC = 0.977**, which clears 0.8, so the diagnostic
  is VALID.
  - Threshold (25th percentile of the ambient anchors' LOAO scores) = 0.852.
  - Median LOAO score: ambient anchors 0.945; T-cell anchors 0.026.
- **Panel A labels:**
  - 7 are `CURATED_LINEAGE_FOREIGN`;
  - 6 are `AMBIENT_FLAGGED`: PIGR, EEF1G, AGR2, SLC34A2, KRT7, MDK;
  - 2 are `NOT_FLAGGED`: BTG1, POLR2J3.
- **Panel B:** unaffected (`CURATED_T_CELL`).

## 3. Eligibility (gene token in >= 10 of a donor's 100 tumour cells; eligible if >= d_min donors)

### Panel A: 2 of 15 eligible (d_min 10)

| Gene | Estimable donors | Result |
|---|---|---|
| POLR2J3 | 40 | eligible, stratum A00, controls OK |
| BTG1 | 43 | eligible, but stratum A01 is **NOT_ESTIMABLE_CONTROLS** (fewer than 20 genes match its very high detection and top rank) |
| WFDC2 6, SCGB3A2 8, EEF1G 8, KRT19 5, KRT18 4, KRT8 4, MUC1 4, AGR2 2, MDK 2, KRT7 2, SLC34A2 2, PIGR 1, HBA1 1 | 1–8 | **NOT_RUN** |

- **So exactly one Panel A gene (POLR2J3) can produce a tested outcome.**
- **What the 13 NOT_RUN genes mean:** they are detected in the T cells of too
  few donors to test. This is itself registered information. The genes that
  topped July's screen are barely present in T cells at all.
- **The panel-level Panel A secondary** is "not run: n_eligible = 2 < 12".

### Panel B: 34 of 36 eligible (d_min 11)

- FOXP3 and IKZF2: 0 donors, because Tregs were excluded from the cohort by
  design, the same `cell_type_major` CD4/CD8 choice as July. They are
  NOT_RUN.
- All 34 eligible genes have estimable controls.
- The panel-level Panel B secondary runs, since 34 >= 12.

## 4. Control strata (the proposed rule, sent to god before any gene statistic existed)

- **19 strata:** 2 for Panel A, 17 for Panel B.
- **Short strata:** 1, A01 (BTG1), so BTG1 is NOT_ESTIMABLE_CONTROLS.
- **318 unique control genes** (20 per stratum, some shared).
- The registration's "six strata / ~120 controls" was an estimate. The real
  panels span too wide a detection range for six strata to match every
  member within tolerance.

## 5. Gene → host (Amendment 3.1: alternate by sorted Ensembl ID within each family)

- 354 genes are mapped (36 eligible panel genes + 318 controls):
  thinkstation1 182, thinkstation2 172.
- BTG1 is not run, because it has no controls. That leaves 353 genes, each
  run in **both** directions on its own host.

## 6. Cost (upper bound, measured 316M bf16 rates)

| Stage | GPU-h |
|---|---|
| Fine-tunes | 3.05 |
| Embeddings | ~0.5 |
| ISP: 353 genes × 2 operations × 43 donors × 100 cells | **~57.4** |
| **Total** | **~61** (about 30 h per host running in parallel) |

The registration estimated ~31 GPU-h. The increase is almost entirely
controls (318, not ~120).

## 7. Found and fixed during prep (no GPU affected)

- **thinkstation2's fresh Geneformer clone** had its **gene dictionaries as
  git-LFS pointers**, because only the 316M weights had been pulled. They are
  now pulled and verified byte-identical to thinkstation1:
  - token dictionary `67c445f4`;
  - gene-median dictionary `a51c53f6`;
  - Ensembl mapping `0819bcbd`;
  - gene name/ID dictionary `fabfa0c2`.
  - Any perturbation run on thinkstation2 would have failed or misbehaved
    before this fix.

## 8. What this means before spending the GPU

- **Panel A** can at most give a result for **one gene (POLR2J3)**. The
  July-replication question is already answered in advance for 13 of 15
  genes: *not testable in T cells*.
- **Panel B** (34 T-cell genes, both directions, concordance) keeps its full
  value.
- **The cost doubles** because of controls.
