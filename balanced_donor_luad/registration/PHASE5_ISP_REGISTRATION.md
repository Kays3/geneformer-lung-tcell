# Phase 5 — ISP pre-registration

**REGISTERED 2026-09-24 by god (message 2026-09-24T11-35-00Z), from draft v2 (sha256 prefix 265d114b). Append-only from registration: changes only as dated amendments at the end.**

- **Status:** DRAFT v2. It incorporates god's rulings of
  2026-09-24T11-05-00Z: Panel B added, deletion AND overexpression,
  concordance as an outcome dimension, and the ambient limitation stated
  as a premise limit.
- **Registration:** it is registered only when god accepts it. After that it
  is append-only; changes go in dated amendments.
- **What exists when this is written:**
  - a tokenised cohort;
  - no classifier;
  - no embedding;
  - no perturbation;
  - no expression statistic for any panel gene in this cohort.
- **What was read to write it:**
  - the July panel source table (gene IDs and ranks);
  - `current_workflow/METHODS.md` s.6;
  - the s100 analysis-layer README;
  - the ambient diagnostic's code;
  - `reports/s100a8-s100a9-crossref-20260922.md`, for the diagnostic's
    circularity;
  - the cohort's gene table and the V2 token dictionary. These were used to
    map symbols to Ensembl IDs, not to read expression.

## 0. Questions (two panels, two separate families)

- **Panel A (replication):** do the genes whose deletion moved T cells from
  LUAD toward normal in the July three-state screen move **tumour T cells
  toward the same donor's adjacent-normal T-cell state**, in a design where
  study, chemistry, site and between-donor variation cannot contribute? And
  are those moves dose-coherent?
- **Panel B (T-cell-intrinsic):** does the classifier read **any
  T-cell-intrinsic signal** in this contrast? This is tested by perturbing the
  repo's curated T-cell anchors directly. The July screen found none of them
  in its top 120.

> **Not an independent-donor replication.** The July LUAD class and this
> cohort draw on overlapping LuCA donors. The design changes, not the people.

## 0a. LIMITATION OF THE WHOLE WORKSTREAM, registered with the claims

**Pairing removes between-donor, between-study, between-site and
between-chemistry variation. It does NOT remove within-donor, between-tissue
contamination.**

Tumour vs adjacent normal is precisely a between-tissue contrast. Tumour
tissue sheds more epithelial and tumour RNA than adjacent normal tissue in
the same donor, so ambient RNA can differ systematically across the contrast
even when every donor is compared with themselves.

- Any result on a gene whose signal could come from ambient RNA is reported
  under the ambient rules below (s.1c). It is never reported as T-cell
  biology.
- Balanced, paired construction does not license an ambient-free
  interpretation.

## 1. Panels (fixed now, by rule; every lookup by Ensembl ID)

### 1a. Panel A — July LUAD→normal top genes (15)

- **Source:** `current_workflow/perturbation_statistics/source_tables/top_goal_shift_genes.csv`
  - sha256 `f7c915c9171fff38917835ff8920f90bd0092c8c68cecdd6b5a90dd9d2ce5416`;
  - last changed c783e3e;
  - read at origin/main 2efcfc7.
- **Rule:** every row with `comparison == "luad_to_normal"`.
- **Checks already done:** all 15 are in the cohort's gene set and in the V2
  token dictionary.

| July rank | Ensembl ID | Symbol | Diagnostic anchor? |
|---|---|---|---|
| 1 | ENSG00000162896 | PIGR | — |
| 2 | ENSG00000254772 | EEF1G | — |
| 3 | ENSG00000111057 | KRT18 | ambient anchor |
| 4 | ENSG00000185499 | MUC1 | ambient anchor |
| 5 | ENSG00000106541 | AGR2 | — |
| 6 | ENSG00000101443 | WFDC2 | ambient anchor |
| 7 | ENSG00000171345 | KRT19 | ambient anchor |
| 8 | ENSG00000157765 | SLC34A2 | — |
| 9 | ENSG00000164265 | SCGB3A2 | ambient anchor |
| 10 | ENSG00000170421 | KRT8 | ambient anchor |
| 11 | ENSG00000135480 | KRT7 | — |
| 12 | ENSG00000133639 | BTG1 | — |
| 13 | ENSG00000168255 | POLR2J3 | — |
| 14 | ENSG00000206172 | HBA1 | ambient anchor |
| 15 | ENSG00000110492 | MDK | — |

### 1b. Panel B — curated T-cell-intrinsic anchors (36 perturbable of 39 listed)

- **Source:** `KNOWN_TCELL` in
  `sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py`
  - file sha256 `6be4bebe69a21a6fc1f4f4a9a76743ae2ccdcc219cbd4cc2fc56943f5862be0c`;
  - last changed d69f56b;
  - the same bytes at origin/main 2efcfc7.
- **Rule:** each of the 39 symbols is mapped to an Ensembl ID through the
  LuCA `raw/var.feature_name`. All 39 map uniquely.
  - **Perturbable** = the mapped Ensembl ID is in the V2 token dictionary
    (`token_dictionary_gc104M.pkl`). That gives **36 genes**.
  - **TRAC** (ENSG00000277734), **TRBC1** (ENSG00000211751) and **TRBC2**
    (ENSG00000211772) are not in the token dictionary. They are kept as rows
    with status `NOT_RUN: not in V2 token dictionary`.
- **Full list:** the 36 Ensembl IDs are frozen in `panel_B.json`, committed
  with the registration.

### 1c. Ambient classification (computed on this cohort before any GPU run; frozen, hashed, committed)

- **Circularity, noted before it can mislead.** The existing diagnostic
  trains a classifier on its own anchor lists. Its recorded
  `cv_auc_mean = 1.0 ± 0.0` on the previous dataset is what a circular score
  looks like (crossref report s.2b(d)). Its validity check (CV AUC >= 0.8) is
  therefore close to vacuous.
- **Replacement, non-circular:**
  - **Validity:** the diagnostic is refit **leaving each anchor out** and that
    anchor is scored by the model that did not see it. The diagnostic is
    valid on this cohort only if the **leave-one-anchor-out AUC is >= 0.8**.
  - If it is not valid, every non-anchor gene gets `ambient: UNDETERMINED`.
    An undetermined gene is treated like a flagged one: it can never reach a
    plain positive status.
- **Non-anchor genes:** flagged when `ambient_risk` >= the 25th percentile of
  the ambient anchors' **leave-one-out** scores. This is the diagnostic's own
  threshold rule, computed without circularity.
- **Anchor genes:** anchors are labelled by curation, not by score.
  - **Panel A ambient anchors** (7 genes above) are
    `ambient: CURATED_LINEAGE_FOREIGN`, meaning a T cell does not transcribe
    them. They are treated as flagged.
  - **Panel B genes** are `ambient: CURATED_T_CELL`, the negative anchor set.
    They are not assessed by the diagnostic, which could not assess them
    independently anyway.

## 2. Classifier gate (Phase 4 output; must pass before any ISP)

- **Model and recipe:** Geneformer-V2-104M, with the recipe fixed in
  PHASE4_COSTING.md. The label is `origin`.
- **Load-bearing rules:** both tissues of a donor are always in the same
  partition, and eval is never used for selection.
- **Gate:** held-out **balanced accuracy >= 0.60**, **and** an exact two-sided
  sign test over held-out donors, testing whether per-donor balanced accuracy
  exceeds 0.5, gives p <= 0.05.
- **If it fails:** stop. There is no ISP, and the result is reported as
  "classifier gate failed", with the numbers.

## 3. Perturbation and per-donor value (both directions)

- **Source cells:** each held-out donor's **100 analysis `tumor_primary`
  cells**. Under k-fold, every donor is scored by the model of the fold that
  held that donor out.
- **Operations:** **deletion** and **overexpression**, as separate runs for
  every panel and control gene.
- **Goal state (primary, paired):** the mean CLS embedding of the **same
  donor's** `normal_adjacent` pool cells (<= 300). These cells are held out
  from the scoring model and are never perturbed.
- **Per-cell shift:** cosine similarity to the goal after perturbation, minus
  before (METHODS s.6). It is computed in the cells the operation applies to,
  which are the cells where the gene's token is present.
- **Donor value:**
  - s_op(g,d) = mean per-cell shift over donor d's perturbed cells.
  - It is computed from the raw per-cell output, **never** from
    `Shift_to_goal_end`.
- **Control-adjusted value:**
  - a_op(g,d) = s_op(g,d) - median over c in C(g) of s_op(c,d), using the
    same operation for the controls.
  - C(g) is **20 matched controls** from `matched_controls.py` logic, matched
    on detection fraction and median token rank in the held-out tumour
    analysis cells.
  - Controls are drawn **per stratum and shared** across genes in that
    stratum: six strata of 20, which gives about 120 control genes. Seed is
    20260924. Panel genes and anchors are excluded.
  - Below 20 candidates in a stratum, the gene is `NOT_ESTIMABLE_CONTROLS`.
- **No-op gate:** before the real run, a no-op perturbation (`run_noop_gene`)
  on 1 donor must give a shift of exactly 0.

## 4. Eligibility (evaluated after dataset selection and before ISP; P4)

- **Estimable (g, d, op):** the gene's token is present in **>= 10 of donor
  d's 100** tumour analysis cells.
- **Eligible (g, op):** estimable in **>= d_min held-out donors**.
  - d_min is found by scanning attainable minimum p = 2/2^d against Holm's
    strictest step at family size m.

| Family | m | Holm strictest step | Attainable p near the step | d_min |
|---|---|---|---|---|
| Panel A, each operation | 15 | 0.003333 | d = 9: 0.003906 fails; d = 10: 0.001953 passes | **10** |
| Panel B, each operation | 36 | 0.001389 | d = 10: 0.001953 fails; d = 11: 0.000977 passes | **11** |
| (not used) A+B combined | 51 | 0.000980 | d = 11: 0.000977 passes | 11 |

- **Why two families, not one:** Panel A replicates a result; Panel B asks a
  different question. A single family of 51 would give the same d_min (11)
  for Panel B and raise Panel A's from 10 to 11. So the separate-family
  choice is on meaning, not on cost.
- **Ineligible (g, op):** kept as rows with status
  `NOT_RUN: estimable donors = k < d_min`. They are never dropped and never
  reported as "no effect".

## 5. Tests

### 5a. Per gene and operation (n = donors)

- **Test:** the exact two-sided Wilcoxon signed-rank test on
  {a_op(g,d) : d estimable}.
  - Exact null by dynamic programming over integer (doubled) ranks.
  - Zeros are dropped and ties get midranks.
  - Critical values are found by scanning attainable values.
- **Multiplicity:** Holm within each of the four families (A-del, A-ovx,
  B-del, B-ovx), each at familywise alpha = 0.05, two-sided.
- **Floor:** stated so it is not dropped silently. The minimum attainable p
  (6.1e-5 for design A, 2/2^43 for k-fold) **does not bind**.
  **Multiplicity and effect size bind.**
- **Effect reported:** median a_op(g,d) with an exact 95% CI (Walsh
  averages), plus the count of donors with a > 0.

### 5b. Dose concordance (an outcome dimension, not a sensitivity)

For each gene, from its two arm results:

| Concordance | Condition |
|---|---|
| `COHERENT` | both arms Holm-significant with **opposite** signs of median a |
| `INCOHERENT` | both arms Holm-significant with the **same** sign (the S100A8/A9 signature) |
| `UNRESOLVED` | exactly one arm Holm-significant |
| `NONE` | neither arm significant |

A positive claim needs **both** arms significant, so it is an
intersection-union test at alpha = 0.05. The conjunction does not inflate the
familywise error of the claim.

### 5c. Outcome rows (every gene always reported, both panels)

Here "toward" means the deletion arm's median a > 0, i.e. deletion moves
tumour T cells toward the donor's normal-adjacent state.

**Panel A:**

| Status | Condition |
|---|---|
| `REPLICATED` | deletion Holm-significant and toward; `COHERENT`; ambient not flagged, not curated-foreign and not undetermined |
| `REPLICATED_AMBIENT` | as above, but ambient-flagged, `CURATED_LINEAGE_FOREIGN` or `UNDETERMINED`. Reported as *consistent with ambient RNA, not T-cell biology* |
| `DOSE_INCOHERENT` | deletion Holm-significant (either sign); `INCOHERENT`. Not a driver |
| `DELETION_ONLY` | deletion Holm-significant; overexpression not significant (`UNRESOLVED`). **Not REPLICATED** |
| `REVERSED` | deletion Holm-significant and away; `COHERENT` |
| `OPEN` | deletion not significant. **Mandatory row, not a failure** |
| `NOT_RUN` | not eligible, or not in the token dictionary (reason stated) |
| `NOT_ESTIMABLE_CONTROLS` | fewer than 20 matched controls |

**Panel B** uses the same logic with neutral names, because either direction
is informative:
- `T_CELL_SIGNAL_TOWARD` / `T_CELL_SIGNAL_AWAY` (Holm-significant and
  `COHERENT`);
- `DOSE_INCOHERENT`, `DELETION_ONLY`, `OPEN`, `NOT_RUN`,
  `NOT_ESTIMABLE_CONTROLS`.

A Panel B result of **zero T-cell signal among the eligible genes** is a
registered finding: the model reads no T-cell-intrinsic signal in this
contrast. It is not a failed experiment.

## 6. Panel-level secondaries (n = genes; P3 floor of 12 applies)

- **Test:** for each panel, an exact two-sided sign test over eligible genes,
  on whether the median deletion effect is toward normal.
- **Floor:** runs only if **>= 12** genes are eligible. Otherwise the status
  is **`not run: n_eligible = k < 12`**.
- **Ambient split:** Panel A is reported separately for ambient-flagged and
  non-flagged genes. Each subset needs >= 12 on its own; with 15 genes in
  total, **at most one subset can reach 12**, and this is stated now.
- **Also reported, count only:** the number of `COHERENT` genes per panel.

## 7. Sensitivities (reported; only S3 can change a status, and only downward)

- **S1.** Goal = global `normal_adjacent` centroid from training donors (the
  July-style goal).
- **S2.** Leader_Merad (23 donors) vs the other four studies.
- **S3.** *If k-fold:* per-fold median a. A positive status (`REPLICATED`,
  `REPLICATED_AMBIENT`, `T_CELL_SIGNAL_*`) whose deletion median has the
  opposite sign in more than one fold is **downgraded to `OPEN`** and
  labelled "fold-heterogeneous". It never upgrades.
- **S4.** Cell-weighted `Shift_to_goal_end`, shown next to the donor value.

## 8. Stop conditions

- The classifier gate fails.
- The no-op gate gives a nonzero shift.
- The ambient computation cannot be produced (as opposed to being invalid,
  which is the `UNDETERMINED` branch).
- The GPU budget is projected to overrun. Nothing is dropped to fit: no
  control, gene, donor or direction.

A partial arm is never analysed.

## AMENDMENT 1 — 2026-09-24T09:47:20Z — model and precision changed by the human's directive

- **Human directive, 2026-09-24, verbatim:** "all future work use bf16 model
  geneformer 316M". It was given directly in my session, and god was
  informed.
- **Replaces in s.2** "Geneformer-V2-104M, with the recipe fixed in
  PHASE4_COSTING.md" with the following. Nothing else changes: not the
  thresholds, the panels, the outcome rows, the tests or the gates.
  - **Model:** Geneformer-V2-316M. Weights come from
    `~/workspace/geneformer-uv-starter/Geneformer/Geneformer-V2-316M/model.safetensors`,
    sha256 `965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3`.
    That hash equals the git-LFS pointer's oid. The `/srv/lab/geneformer`
    copy is a 135-byte LFS pointer and is **not** used.
  - **Code:** the vendored Geneformer checkout `f45a6c7` plus
    `sclc_validation/bf16_bench/geneformer-bf16-export.patch`. The checkout
    was verified to carry exactly that patch (reverse dry-run clean, 9
    changed lines).
  - **Precision:** bf16.
    - Training uses the Trainer's `bf16: true`.
    - ISP uses `dtype_cast.install_dtype_cast("bf16")`.
    - Held-out evaluation for the s.2 classifier gate uses
      `eval_tie_fix.install_eval_tie_fix`, as in `run_finetune_316m.py`, and
      reports the tie counts.
  - **Unchanged recipe:** 1 epoch, lr 5e-5, batch 8, freeze_layers 6, seed
    43, label `origin`, no search, eval not used for selection, and both
    tissues of a donor always in the same partition.
- **Tokenisation stays valid.** The existing tokenised cohort was made with
  `/srv/lab/geneformer` `04c2b2e`. Its token, gene-median and
  Ensembl-mapping dictionaries are byte-identical to `f45a6c7`'s (sha256
  prefixes 67c445f4, a51c53f6, 0819bcbd), so it is not re-tokenised.
