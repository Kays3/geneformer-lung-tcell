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

## AMENDMENT 2 — 2026-09-24T09:51:03Z — Panel A model confound: limitation + CONDITIONAL 104M control arm (requested by god, 2026-09-24T13-20-00Z)

This amendment is registered before any Phase 6 number exists. Nothing in
s.1–s.8 or Amendment 1 changes. This adds one limitation and one conditional
arm.

### Limitation (registered, applies to Panel A)

- **Model of record for the July table:** the Panel A source table comes from
  the July screen, whose base model is **Geneformer-V2-104M**
  (`current_workflow/METHODS.md`, s.4 table, "Base model").
- **Model for this run:** under Amendment 1, Panel A runs on
  **Geneformer-V2-316M**, bf16. That is a different model, not only a
  different design.
- **Quantified bound on the nuisance term,** from this programme's own
  measurement (`sclc_validation/bf16_bench/RESULTS_BF16_316M.md`, "Result",
  origin/main):
  - 104M-vs-316M, both bf16, on an overlapping ISP task: **Spearman
    rho = 0.489**, top-20 row overlap 9/20.
  - Same-model precision change: rho = 0.9998
    (`sclc_validation/bf16_bench/RESULTS_BF16.md`).
  - That was a different task (a 50-gene overexpression panel on the
    SCLC/LUAD/normal classifier). So it bounds the size of a model effect;
    it does not measure it here.
- **Consequence:** a Panel A non-replication on 316M alone could be caused by
  the design change or by the model change, and cannot be attributed.

### Conditional 104M control arm (pre-registered; triggered by outcome, not by choice)

- **Trigger:** after the 316M Phase 6 statuses exist, **if at least one
  eligible Panel A gene has a status other than `REPLICATED` or
  `REPLICATED_AMBIENT`**, meaning `OPEN`, `DELETION_ONLY`,
  `DOSE_INCOHERENT` or `REVERSED`, then a **104M Panel A control arm is
  required before any Panel A non-replication is reported.**
  - Genes that are `NOT_RUN` or `NOT_ESTIMABLE_CONTROLS` do not trigger it.
  - If every eligible Panel A gene is `REPLICATED` or `REPLICATED_AMBIENT`,
    the arm is **not run**, and the reason is recorded.
- **Specification:**
  - **Model:** Geneformer-V2-104M, weights sha256 prefix `fff5cba29ddd8792`.
    The `/srv/lab` copy and the vendored copy are byte-identical, real
    weights.
  - **Precision:** bf16, per the human's directive. Precision moves results
    by rho 0.9998, so bf16 does not reintroduce a nuisance.
  - **Same everything else:** the same held-out donors, cells, split, goal
    construction, matched controls, tests, thresholds, eligibility, Holm
    family and outcome rows as the 316M Panel A run. No re-selection of
    anything.
- **Attribution table:** applied per gene, for genes that are non-replicated
  on 316M.

| 104M arm status for that gene | Reported as |
|---|---|
| `REPLICATED` / `REPLICATED_AMBIENT` | non-replication **attributable to the model change** (104M reproduces the July finding under the paired design; 316M does not) |
| any other eligible status | non-replication **under the paired design, robust to model** (neither model replicates) |
| `NOT_RUN` / `NOT_ESTIMABLE_CONTROLS` on 104M | **unattributable**, reported as such |

- **If the arm is triggered but not approved or not run:** Panel A
  non-replications are reported **only** as "non-replication on 316M;
  model-vs-design attribution not established".
- **Cost** (upper bound, from the measured 316M rates; 104M is smaller):
  about **8.6 GPU-h**. It needs its own human approval when triggered.
- **Directive note:** this arm uses 104M, not 316M. Whether "all future work
  use bf16 model geneformer 316M" allows a same-model control is the human's
  call. God is raising it with them.

## AMENDMENT 3 — 2026-09-24T10:24:17Z — human directives, 3-class target, LUSC screening rule, cross-host gate (accepted by god 2026-09-24T14-30-00Z with four additions, incorporated; host split chosen by the human)

This amendment is registered before any Phase 6 number exists, and before
any LUSC candidate has been screened. The read-only LuCA check that showed
the 10x LUSC shortfall (9 donors) used structural fields only.

### 3.1 Human directives (2026-09-24, verbatim, given in my session)

- **"do only 316M model avoid 104M".**
  - The conditional 104M Panel A control arm of Amendment 2 **will not be
    run**.
  - Amendment 2's pre-registered fallback applies. Any Panel A
    non-replication is reported **only** as "non-replication on 316M;
    model-vs-design attribution not established".
- **Host use.** The directive was "use thinkstation1 for deletion and
  thinkstation2 for overexpression". God pointed out that an operation-based
  split confounds host with the delete/overexpress axis, which is exactly what
  the registered dose-concordance test compares. Asked to choose, the human
  chose **split by gene** (2026-09-24).
  - **Each host runs BOTH operations** for its own genes.
  - The concordance comparison is therefore always within one host, so a host
    difference cannot enter it.
  - **Assignment (deterministic, fixed now):**
    - within each family (Panel A, Panel B, each control stratum), genes are
      sorted by Ensembl ID;
    - they alternate thinkstation1, thinkstation2, thinkstation1, …, starting
      with thinkstation1;
    - so each host gets half of each family, within one gene.
  - **Stack:** both hosts are identical: Geneformer f45a6c7 + bf16 export
    patch, 316M weights sha `965ceccea81953d3…`, the same `uv.lock`, and
    identical sorted freeze hashes (provenance/ENVIRONMENTS.md). The drivers
    differ (595.84 vs 580.173.02).
- **"no hard limit for GPU, use at will".**
  - The design is **k = 5 donor cross-fitting** (the option registered
    earlier). s.7 S3 is active.
  - GPU hours are metered and reported. They are not a stop. Nothing is
    dropped for cost.
- **"Three classes like July"** (chosen by the human), with a search beyond
  LuCA for LUSC. See 3.2–3.4.

### 3.2 Target design: 3 classes

- **Classes:**
  1. LUAD primary tumour;
  2. LUSC primary tumour;
  3. adjacent-normal tissue.
- **Unchanged:** ISP start state = LUAD tumour; goal = the same donor's own
  adjacent-normal centroid (s.3).
- **LUSC centroid:** computed from training-fold LUSC donors and recorded as
  the alternative state, as in July.
- **Panels, tests, eligibility, outcome rows and the ambient rules:** all
  unchanged.

### 3.3 LUSC screening rule (fixed before screening)

- **Criteria:** Phase 0 A1–F6 unchanged, applied to LUSC tumour donors
  (>= 100 T cells, sites per F2, public raw counts per F5–F6). Plus:
- **F3-LUSC.** LUAD vs LUSC cannot be paired within a donor, so:
  - every study contributing LUSC donors must also contribute qualifying
    LUAD tumour donors on the same chemistry;
  - otherwise the LUSC class is a study label, and that study is excluded;
  - adjacent-normal donors for the normal class must come from studies that
    also contribute tumour donors (the existing F3).
- **F4** is applied across all three classes. Every chemistry present must
  appear in every class.
- **A1-LUSC:** >= 12 qualifying LUSC tumour donors (16 preferred) across the
  studies that pass F3-LUSC and F4.
- **Stop rule:** if fewer than 12 pass, the result is "3-class not feasible
  under the registered rule".
  - The human is asked before any relaxation.
  - The default is the registered 2-class design, run unchanged.
- **Recorded per candidate:** donors per histology per tissue with >= 100 T
  cells, chemistry, LUAD overlap by study, download size, and the first
  criterion failed. Nothing is looked at beyond structural fields.

### 3.4 3-class classifier gate (replaces s.2's numbers only if 3-class goes ahead)

- **Balanced accuracy:** held-out, pooled across folds, >= 0.60 (chance
  0.333).
- **Per-donor test:** an exact two-sided sign test over held-out donors,
  where each donor's balanced accuracy over its own true class(es) exceeds
  1/3, with p <= 0.05.
- **If it fails:** STOP.

### 3.5 Cross-host equivalence gate (before AND after Phase 6)

- **Probe:** one fixed (gene, donor), with both operations run on both hosts,
  using identical sha256-verified fold models and centroids. The probe gene
  is B2M (ENSG00000166710, in neither panel), on the first donor by sorted ID
  in fold 1.
- **Pass:** for each operation, per-cell shift max |Δ| <= 1e-3 **and**
  Spearman rho across cells >= 0.999. Two measures are used because rho
  survives a constant offset and the absolute bound does not.
- **When:** run **before Phase 6** (a failure means STOP: no Phase 6) and
  **again at the end** on the same probe, reporting both. The pair bounds
  drift across a long run.
- **If the end gate fails:**
  - pooled results are still reported, but every outcome row carries
    `host_drift: true`;
  - per-host outcome rows are reported alongside;
  - no pooled claim is made without them.
- Because the split is by gene, the gate protects **pooling across genes**.
  It does not protect concordance, which cannot be host-confounded.

### 3.6 The LUSC asymmetry (registered before any number)

> **LUAD↔normal remains within-donor paired. Any contrast involving LUSC is
> between-donor, guarded by study overlap on matched chemistry. Claims
> involving LUSC carry weaker warrant than claims about LUAD↔normal, and are
> to be reported as such.**

- Every outcome row and figure involving LUSC carries the field
  `design_warrant: "between-donor (study-overlap guarded)"`.
- LUAD↔normal rows carry `design_warrant: "within-donor paired"`.
- A table or figure mixing the two must show the field.

### 3.7 Cost: no cap, but an estimate and an account

- **Before Phase 6:** a **3-class k = 5 GPU-hour estimate**, from the measured
  316M bf16 rates and the actual LUSC cell counts, is produced and sent to
  god and the human.
- **During the run:** GPU-hours are metered, with `nvidia-smi` sampled every
  60 s on both hosts.
- **At close:** actual vs estimate is reported per stage and per host.
- The absence of a cap never licenses dropping anything, and it never
  licenses skipping the account.

### 3.8 Panel A qualifier is a field, not prose

- Every Panel A outcome row whose status is `OPEN`, `DELETION_ONLY`,
  `DOSE_INCOHERENT` or `REVERSED` carries the mandatory field:
  - `claim_qualifier: "non-replication on 316M; model-vs-design attribution not established"`.
- The row renderer (outcome table, figures, slides, RESULTS.md) **fails** if a
  Panel A non-replication row lacks it. A unit test shows that failure.
- This follows Amendment 2's fallback, which binds because the human declined
  the 104M arm.

## AMENDMENT 3a — 2026-09-24T10:28:10Z — timestamp correction for Amendment 3's acceptance (the header is not edited)

- **What was wrong.** Amendment 3's header cites god's acceptance as
  "2026-09-24T14-30-00Z". That is a hive **message id**, which god has
  confirmed was hand-written and runs about 4 h fast. It is **not** system
  UTC.
  - Read naively, the header says acceptance (14:30Z) came *after* this
    registration (10:24:17Z).
  - Hive message ids are not system UTC in general. God's ids after
    2026-09-24T10:26:52Z are generated from `date -u`.
- **Ordering evidence (verified by me, not taken on report).** The
  filesystem mtime of my delivered copy of god's acceptance message is:

  ```
  inbox/.done/2026-09-24T14-30-00Z-god-phyllis-amendment3.json   mtime 2026-09-24T10:19:19Z
  Amendment 3 registration timestamp (system UTC)                2026-09-24T10:24:17Z
  ```

  - Acceptance precedes registration by about 5 minutes.
  - The file's own `created_at` field repeats the wrong 14:30:00Z; only the
    mtime is a real clock.
- **Why the header is not edited.** Replacing one timestamp with another in
  place would hide the error. This dated note is the correction.

## AMENDMENT 3b — 2026-09-24T10:51:53Z — design settled as 2-class by the human (no numbers exist)

- **The LUSC survey** (provenance/PHASE1B_LUSC_SURVEY.md) found the 3-class
  outcome depends on how F3-LUSC's phrase "same chemistry" is read:
  - by chemistry family: 12 LUSC donors, passing exactly at the floor;
  - by exact version: 10, failing.
  - Neither reading was chosen by me.
- **The human chose "Go 2-class (registered)"** on 2026-09-24, in session.
  Therefore:
  - The design is the **registered 2-class LUAD tumour vs same-donor
    adjacent normal**: 43 paired donors, 10x, as built in Phase 3.
  - s.2's classifier gate applies. Amendment 3.2–3.4 (3-class target, LUSC
    rule, 3-class gate) are **not used**.
  - No row, figure or claim involves LUSC. The `design_warrant` field is
    therefore always "within-donor paired". The field and its validator stay
    in place.
  - Unchanged from Amendment 3: 316M only (no 104M arm; Amendment 2's
    fallback binds, with the `claim_qualifier` field), k = 5 cross-fitting,
    host split by gene, the equivalence gate before and after, and cost
    metered with an estimate and an account.
- **Estimate** (upper bound, measured 316M bf16 rates): about 31 GPU-h.
  - fine-tunes about 3.05 h;
  - ISP about 27.8 h;
  - embeddings about 0.5 h.

## AMENDMENT 3c — 2026-09-24T11:07:38Z — control strata rule and three registered limitations (accepted by god; message mtime-ordered before this commit)

### 3c.1 Strata rule

The rule is fixed **before any per-gene statistic of the held-out cells
existed**. It was proposed to god in that state. That timing is what makes
it a rule and not a choice.

- **Gene stats:** `detect_frac` and median 0-based token rank over every
  donor's 100 held-out **tumour** analysis cells (43 × 100 = 4,300 cells),
  via `s100_isp.matched_controls.median_token_rank_and_detection`
  (per-source-state).
- **Grouping, per family, over ELIGIBLE panel genes only:**
  - sort by (log2 detect, rank percentile);
  - add the next gene to the open stratum only if it is within 0.5 log2
    and 5 rank-percentile points of ALL current members;
  - otherwise open a new stratum;
  - K is set by the data.
- **Controls:** 20 per stratum via `build_matched_control_table` (seed
  20260924), excluding every panel gene and every anchor.
  - Fewer than 20 candidates means `NOT_ESTIMABLE_CONTROLS`, never widened.
  - A control drawn by two strata is run once.
- **Replaces** s.3's "six strata / ~120 controls", which was an estimate,
  not a rule.
- **Per-gene matching (20 controls per gene) is rejected.** Stratum controls
  already lie within tolerance of every member, so per-gene matching would
  buy independence, not exactness, at roughly 170 extra GPU-h.
- **Applied result** (committed `controls/pre_gpu_*`):
  - 19 strata: Panel A 2, Panel B 17. All are built over eligible genes
    only; ineligible genes are excluded before grouping.
  - 1 short: A01, BTG1, so `NOT_ESTIMABLE_CONTROLS`.
  - 318 unique controls, and 353 genes to run.

### 3c.2 Limitation: shared-control dependence affects the panel-level secondary only

- Genes in one stratum share their 20 controls, so their control-adjusted
  values are **not independent**. The panel-level sign test (s.6) assumes
  independent signs, so its p-value is reported with this caveat.
- **The per-gene primary (s.5a) is unaffected,** because it is within-gene
  across donors.
- This is registered, not repaired. Repairing it (per-gene controls) is not
  proportionate for a secondary.

### 3c.3 Limitation: the ambient axis cannot separate contamination from ubiquity

- The diagnostic scores genes on detection-based features. **Ambient
  contamination and ubiquitous expression both produce broad detection.**
- The visible case is EEF1G (ambient_risk 0.999), a translation elongation
  factor expressed in every cell type.
- **Its flag stands:** the diagnostic is registered and is not adjudicated
  gene by gene.
- **Meaning of the label:** `AMBIENT_FLAGGED` means *"resembles known
  contaminants on a detection axis"*, **not** *"is contamination"*. Every
  report of the label carries that meaning.

### 3c.4 Limitation: why 13 of 15 Panel A genes are NOT_RUN has two readings, and the data cannot choose

- **The finding:** 13 Panel A genes are present (>= 10 of 100 cells) in the
  tumour T cells of only 1–8 of 43 donors.
- **Two explanations, not separable here:**
  - (a) **ambient RNA**, which is sporadic across cells and donors;
  - (b) **cohort composition**: July's LUAD class drew heavily on
    non-primary tissue, so a gene present in T cells of metastatic or
    adjacent-normal tissue can be genuinely absent from primary-tumour T
    cells. This cohort is pure primary tumour.
- **What both readings imply:** those genes were not measuring primary-tumour
  T-cell biology. The report names both readings and picks neither.
- **Wording, fixed:** `NOT_RUN: eligibility (estimable donors = k < 10)`,
  never "no effect". The claim is "untestable in this cohort", which is a
  statement about presence, not biology.

## AMENDMENT 3d — 2026-09-24T11:11:23Z — a classifier-gate FAILURE is a reportable finding (registered before the gate has a number)

- **Timing:** registered while Phase 4 was still fine-tuning. The gate
  (`scripts/classifier_gate.py`) had **not** been computed.
- **If s.2's gate fails** (pooled held-out balanced accuracy < 0.60, or the
  per-donor sign test not significant in the right direction):
  - ISP still does not run. The stop stands.
  - The failure is **reported as a result**, with all gate numbers: pooled
    and per-fold balanced accuracy, per-donor values, the sign-test counts
    and p, and tie counts.
- **Registered interpretation of a failure:** a balanced, within-donor-paired,
  single-site, study-clean tumour-vs-adjacent-normal T-cell classifier
  **cannot separate the two states**. Because the design removes the confounds
  the July classifier had, this is the cleanest statement this programme can
  make about whether the model reads this contrast.
- **Bound, stated with it:** the statement is about *this model, cohort and
  recipe* (Geneformer-V2-316M bf16, 43 donors, one atlas, the fixed recipe).
  It is **not** a statement about tumour-vs-normal T-cell biology in general.
- **The same applies to a pass:** it is reported with the same numbers,
  whichever way the gate goes.

## Amendment 3e — reading rules against measured nondeterminism (2026-09-24, ~18:15Z)

**Registered before:** the classifier gate has been computed (only fold 0 of 5 exists), and before any ISP or equivalence-probe run. It changes no threshold and no decision. It fixes how a result near a threshold, or a probe failure, is read. The trigger is the same-seed fold-0 comparison in `provenance/PHASE4_RUN_NOTES.md`.

### 3e.1 Measured reference (fine-tuning, same host, code, seed and data)
Attempt 1 vs attempt 2 fold 0, on 1,800 held-out cells:
- label agreement 95.9%;
- |Δ logit| median 0.47, p99 1.76, max 2.89;
- per-donor balanced-accuracy Δ from −0.045 to +0.015;
- fold-pooled balanced-accuracy Δ +0.0017.

This is ONE repeat of ONE fold (9 donors). It is an observed magnitude, not a bound.

### 3e.2 Classifier gate near its threshold
The registered gate and its decision are unchanged: attempt 2 only, pooled BA ≥ 0.60, AND an exact two-sided sign test p ≤ 0.05 with more donors above 0.5 than below. The PASS/FAIL outcome is exactly what that rule gives. The report must ALSO state:
- (a) **Pooled BA margin.** If the pooled BA lies within ±0.045 of 0.60 (0.555 to 0.645), the result is reported as "PASS/FAIL, within run-to-run nondeterminism of the threshold". 0.045 is the largest per-donor change observed. It is deliberately more conservative than the pooled change (0.0017).
- (b) **Sign-test margin.** Donors whose balanced accuracy is within 0.045 of 0.5 are "flippable". The report gives the sign-test p in two cases: all flippable donors moved against the gate, and all moved in its favour. If the gate decision differs between the two, the sign-test component is reported as "within run-to-run nondeterminism".
- (c) If neither (a) nor (b) applies, the result is reported as a clean PASS or FAIL.

In every case the decision (proceed or stop, Amendment 3d) follows the registered rule, not the band.

### 3e.3 Cross-host equivalence probe: which floor applies
The probe runs INFERENCE with the SAME fine-tuned fold model files on both hosts (copied, sha256-verified). Training nondeterminism (3e.1) therefore does not enter it, and 3e.1 is NOT its reference band. Its floor is same-host INFERENCE repeatability. That is measured alongside, not assumed:
- **Runs:** B2M (ENSG00000166710), fold 1's first donor, both operations. ts1 twice, ts2 twice, same inputs, into separate output folders.
- **Registered criterion, unchanged (Amendment 3):** ts1-vs-ts2 per-cell shifts, per operation, max|Δ| ≤ 1e-3 AND Spearman rho ≥ 0.999. PASS or FAIL is decided by this alone. A FAIL stops Phase 6 and goes to god (and the human).
- **Reported alongside, per operation:** the same-host repeat max|Δ| and rho on each host.
- **Reading a FAIL:**
  - It is attributed to the HOST only if the cross-host max|Δ| exceeds the larger same-host repeat max|Δ|.
  - Otherwise it is reported as "inference nondeterminism exceeds the registered tolerance on a single host". That is a different finding, with a different remedy. Neither reading permits proceeding without a decision.
- **Reading a PASS:** it is reported together with the same-host floor. If both same-host repeats are bitwise identical, the pass is against zero; if not, against that floor.

**Erratum (appended 2026-09-24T18:11:53Z):** the Amendment 3e header says "~18:15Z". The authoritative time is commit b9aafc6 at 2026-09-24T18:08:42Z. The header estimate ran ahead of the clock. No other content is affected.

## Amendment 3f — no donor exclusion on classifier accuracy (appended 2026-09-24T20:50:31Z)

**Registered:** after the classifier gate (PASS, clean; phase4_results/classifier_gate.json at 0fffc96) and BEFORE any ISP output exists. No goal embedding had finished and no ISP call had run. Requested by god.

1. **No accuracy-based exclusion.** No donor is excluded from any Phase 6/7 primary or secondary analysis because of its per-donor held-out classifier balanced accuracy. The donors that enter each gene's test are exactly those set by the eligibility rule already registered (s.3: the gene is present in ≥ 10 of the donor's 100 tumour analysis cells), whatever their per-donor BA.
2. **The accuracy/effect relation is descriptive only.** For every Panel A and Panel B gene with at least one significant arm, and for the pooled panel-level donor values, the report gives the Spearman correlation across donors between per-donor BA (classifier_gate.json) and the donor's control-adjusted shift, per operation. It also flags the six donors with BA < 0.65 in the per-donor plots. This is a sensitivity observation. It never changes a status, a test, or a donor set.
3. **Any later exclusion** needs a dated amendment, accepted by god BEFORE the ISP results are unblinded. Any exclusion proposed after results exist is reported as post hoc and is never applied to the primary.

**Erratum to 3f (appended 2026-09-24T20:50:47Z):** "No goal embedding had finished" is wrong. At commit time (2026-09-24T20:50:32Z), fold 0's global and per-donor goal embeddings already existed on thinkstation1. The claim that matters is unaffected: no ISP call had run and no ISP output existed. Goal embeddings contain no gene-perturbation result.

**Precision to 3f.2 (appended 2026-09-24T20:52:10Z):** the gene-level correlation is conditioned on significance and describes those genes only. The panel-level donor correlation is the unconditioned quantity. No ISP output existed at this append.

## Amendment 3g — three analysis rules left open by s.3, s.5a and 3f (appended 2026-09-25T04:04:06Z)

**Registered before any Phase 6 output is read.** Phase 6 was still running at this append. The Phase 7 code (commit ab9af4b) has been run only on synthetic data. Rulings by god (message 2026-09-25, "three rules ruled"). In each case the proposal was chosen, not the alternative.

### 3g.1 Controls per donor (s.3)
- A control counts toward donor d's median only if it is **estimable in d under the same rule as the gene**: its token is present in >= 10 of d's 100 tumour analysis cells.
- a(g,d) is **undefined if fewer than 10 of the 20** controls qualify, and that donor drops from that gene's test. **The 10-of-20 threshold is a convention fixed in advance. No principled basis is claimed for it.**
- Every row reports the **count of qualifying controls per donor** as well as the list of dropped donors.
- Eligibility is unchanged: s.4 counts estimable donors. If donors drop under 3g.1, the test runs on the rest, and the row reports both counts.

### 3g.2 Holm family size (s.5a)
- m = the **full registered panel**: 15 for Panel A and 36 for Panel B, each operation separately. Untested members enter as p = 1.
- Reason: m was registered. The s.4 d_min table was derived from this m, and choosing m after learning which genes fell out would let the data pick the correction.
- Every row reports the raw p, the Holm-adjusted p and m.

### 3g.3 The 3f panel-level donor value
- Per donor: the median of a_op(g,d) over that panel's tested genes where the donor is estimable. This is correlated with per-donor classifier BA by Spearman, per operation, **descriptive only**.
- The **number of contributing genes per donor** is reported, and the write-up states that the value is a median over a varying number of genes.

### 3g.4 Registered sensitivity line for 3g.1 and 3g.2
After unblinding, the report states **whether any registered status would change** under each alternative rule, one at a time:
- **3g.1 alternative:** any control with >= 1 perturbed cell counts, and 1 qualifying control suffices.
- **3g.2 alternative:** m = the genes actually tested.

This is **reported only. No registered or headline status ever changes on it.** A status that would flip under an alternative is reported as fragile, in the same sentence as the status.

### 3g.5 Launch gate
- The rules are frozen in `registration/phase7_rules.json`.
- `registration/required_gates.json` lists every pre-analysis gate, with the registration line that requires it. The Phase 7 launcher refuses to start unless every listed result file exists and shows the required value.
- The end-of-run equivalence result must **exist** before the launch. Per s.3.5, a FAIL does not block the analysis: every row then carries `host_drift: true`.

## Amendment 3h — overexpress cell set, reconstruction and its identity check (appended 2026-09-25T14:43:32Z)

**Registered after Phase 6 finished and before any outcome value is read.** At this append, Phase 7 had been launched once. It refused at the loader (`IncompleteRun`) before computing any statistic. What has been read since: markers, per-cell vector lengths, cell lengths and tokens, and the B2M end-of-run probe. No panel or control gene's shift value has been read. Rulings by god (messages 2026-09-25 14:35Z, 14:39Z and 14:41Z).

### 3h.1 The finding
- s.3 says the per-cell shift is computed "in the cells the operation applies to, **which are the cells where the gene's token is present**". That is true for delete. **It is not true for overexpress with a gene list.** The vendored Geneformer (`in_silico_perturber.py`, ~L524) skips the token filter when `perturb_type == "overexpress"`, and inserts the gene into every start-state cell.
- So every Phase 6 overexpress pickle holds **100 per-cell values**, where the marker's `n_token_cells` = k.
  - This applies to 14,738 of the 15,179 overexpress calls.
  - The rest are 404 calls with k = 100 and 37 `no_token_cells` calls.
  - Every delete pickle holds exactly k values.
- Both operations then pass through `downsample_and_sort` (`Dataset.sort("length", reverse=True)`). So pickle order is length-descending, not dataset order.
- **Prior art in this repository, missed:** card `isp-runner-paired-arms-20260922` documented this library behaviour and removed it by construction, pre-filtering both operations to token-positive cells. The earlier whole-genome screen (`genes_to_perturb="all"`, PR #16) overexpressed only genes already in each cell, so it was token-positive by construction. The pre-2026-09-22 targeted 50-gene panel overexpressed into every cell.

### 3h.2 Primary: token-positive cells for both operations (Option A)
- s_overexpress(g,d) = mean per-cell shift over the **same k token-positive cells** that delete uses.
- **Reason:** the s.5b concordance test asks whether deletion and overexpression move the classifier in opposite directions. That is only meaningful if both arms perturb the same cells. With all cells, the overexpress arm would be mostly insertion into cells that never expressed the gene (88 of 100 cells at k = 12). The clause in s.3 happens to encode this. The ruling would be the same without it.
- **Cost, stated:** the overexpress arm's n per (gene, donor) falls from 100 to k.

### 3h.3 Sensitivity B (unregistered, labelled)
- Overexpress over all 100 perturbed cells; delete unchanged. It is reported next to the primary and **never changes a status**.

### 3h.4 Reconstruction rule
- Replay Geneformer's own calls on each donor's ISP input: the start-state filter (`origin == tumor_primary`, a no-op because all 100 cells are tumour), then `Dataset.sort("length", reverse=True)`.
- The token-positive subset = the positions, in that order, of cells whose `input_ids` contain the gene's token.
- `datasets` 5.0.1, `pyarrow` 25.0.1 and `torch` 2.13.0 are identical on thinkstation1 and thinkstation2.
- Ties in length resolve by the sort's stability (pyarrow `sort_indices`).
- **Measured:** 12,184 of the 15,142 (gene, donor) pairs with k >= 1 (80.5%) contain a length shared by a token-positive and a token-negative cell.

### 3h.5 Count checks (must-match; refusal on any failure)
- (i) For every delete call, the replayed positive count equals the pickle's per-cell length.
- (ii) For every overexpress call, the reconstructed positive count equals the marker's `n_token_cells`, and the pickle length equals the donor's cell count.
- These test **cardinality only**. They cannot detect the right number of wrong cells.

### 3h.6 GPU identity check (tests identity; stop rule)
- **Pair selection rule, frozen before any pair runs.**
  - Candidates are (control gene, donor) pairs with 10 <= k <= 99.
  - A pair is assigned to the host that ran that gene in Phase 6.
  - "Mixed-tie positives" = token-positive cells whose length is shared with a token-negative cell.
  - Per host, four distinct pairs, chosen in this order. Every tie-break is by ascending (Ensembl ID, donor ID).
    1. **Tie-rich:** the most mixed-tie positives.
    2. **Tie-free positive control:** zero mixed-tie positives, largest k. If this pair also fails, the reconstruction model is wrong in general, not only at ties.
    3. **Low k:** 10 <= k <= 12 and >= 1 mixed-tie positive; the most mixed-tie positives.
    4. **High k:** 90 <= k <= 99 and >= 1 mixed-tie positive; the most mixed-tie positives.
  - An empty category is reported as empty, never substituted.
  - The selected pairs are written to `phase7_prep/identity_pairs.json` and committed **before any of them runs**. No pair is added afterwards.
- **Procedure.**
  - For each pair, build a token-positive-only copy of the donor's ISP input, in dataset order. This is the 09-22 construction, where the order is known.
  - Run `run_isp.py` for that gene and donor on its Phase 6 host, with the same fold model, goals, code path and package set.
  - Compare the overexpress output with the full Phase 6 pickle at the reconstructed positions.
  - Delete is compared the same way and reported.
- **Pass:** for every pair and all three states, max|Δ| <= 1e-3 **and** Spearman rho >= 0.999. This is the tolerance of the equivalence gates.
  - **Bit identity is not expected.** The forward batch is padded differently: k cells versus 100.
  - The tolerance is the criterion; 0.0 is not.
  - The known failure signature for misalignment is a comparator control with reversed order: Δ ≈ 0.14, rho ≈ 0.
- **Stop:** if any pair fails, Phase 7 does not run and the matter goes to god.
- **Scope:** up to eight control genes' shift values are read, as a technical check, not an outcome. Controls are not registered outcomes, and the analysis rules are frozen in `phase7_rules.json` and compared byte for byte at launch.

### 3h.7 Registered fallback C
- If 3h.6 fails, overexpress is re-run by construction on token-positive-only inputs for all 15,142 pairs. The end-of-run equivalence is then repeated.
- Before committing to C, one gene is timed, and the measured estimate goes to god. The unmeasured 3-4 h figure is not used.

### 3h.8 Withdrawn instruction, stated plainly
- An instruction given at 14:35Z, to **refuse on any positive/negative length tie**, is **withdrawn, not amended**.
- It was keyed to the wrong condition. A tie makes the subset ambiguous only if the sort is unstable, and 3h.6 tests exactly that.
- Refusing on ties would have refused 80.5% of the analysis.

### 3h.9 Launch gate additions
- `required_gates.json` gains two gates:
  - `ovx_count_checks`: `phase7_prep/ovx_index_checks.json` has `failures` = 0;
  - `ovx_identity_check`: `phase7_prep/identity_check_result.json` has `PASS` = true.
- The launcher passes the committed position index to the analysis.

### 3h.10 Late comparator control, disclosed
- After the end-of-run gate had passed, `probe_cmp.py` was shown to detect planted differences:
  - +1e-6 in one cell: detected, within tolerance;
  - +2e-3 in one cell: FAIL;
  - reversed order: FAIL, rho -0.16 to 0.02.
- This control was run late, not before the gate.
