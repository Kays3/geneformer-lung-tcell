# E2 registration: the balanced-donor Geneformer design on colorectal cancer (Pelka et al. 2021)

**Status: REGISTERED, awaiting the independent gatekeeper (Stanley). No model has been fine-tuned, no
embedding extracted and no perturbation run for E2. Nothing in this file may change after the first GPU
call except the dated deviation log in s.12.**

- Written by Phyllis, 2026-10-01 JST, on branch `analysis/e2-pelka-crc-20261001` of `Kays3/geneformer-lung-tcell`.
- Authorisation: the human's GO of 2026-10-01 about 00:05 JST, relayed by god (conversation
  `geneformer-e2-20261001`): E2 on Pelka 2021 with the 19-donor unsorted-only cohort, E2 GPU time (design
  estimate 18 to 52 GPU-hours) and the GSE178341 count download.
- Design source: `hive/reports/geneformer-alt-dataset-experiment-design-20261001.md` s.6 (gate-passed),
  and the E0 feasibility count `hive/reports/geneformer-e0-feasibility-20261001.md` (gate-passed).
- Everything computed so far is CPU-only and outcome-free: the cohort draw, tokenisation, the fold split,
  panel eligibility, matched controls, the null-gene draw and the run order. Their files are committed with
  this registration and pinned by hash in s.4.

## 1. Question and hypotheses

The balanced-donor LUAD study produced three results we want to test outside lung
(`balanced_donor_luad/RESULTS.md`, the IMRaD report `hive/reports/balanced-donor-null-imrad-report-20260930.md`):
R1, a fine-tuned 316M classifier separates tumour from adjacent-normal T cells of the same donor (pooled
held-out balanced accuracy 0.825, 43 of 43 donors above chance); R2, 13 of the 34 testable curated Panel B
genes reached a stable, dose-concordant status; R3, across random eligible genes detectable in LUAD tumour
T cells, deletion and overexpression shifts are negatively rank-correlated (rho -0.593 at N=100, p at the
permutation floor of 100,000 permutations; confirmed at N=200, rho -0.608).

E2 re-runs the same design, from the base Geneformer-V2-316M weights, on colorectal tumour and normal-colon
T cells of 19 donors. It asks whether R1 to R3 are properties of the model and the design, or of lung.

- **H2a (R1):** a donor-paired tumour/normal classifier fine-tuned by the same recipe separates colorectal
  tumour from normal-colon T cells. Test: the balanced-donor classifier gate, unchanged (s.6.1).
- **H2b (R3), primary:** deletion and overexpression shifts are negatively rank-correlated across random
  eligible genes detectable in colorectal tumour T cells. Test: the Amendment 5 null analysis, unchanged
  (s.6.2).
- **H2c (R2):** the LUAD Panel B pattern holds: the 13 LUAD reference genes keep the direction of their
  deletion shift, and the per-gene statuses are reported beside the LUAD ones (s.6.3).

## 2. Data and cohort

**Source.** GEO GSE178341 (Pelka et al. 2021, *Cell*), downloaded to thinkstation1 on 2026-09-30
15:06-15:07 UTC (2026-10-01 00:06-00:07 JST). Raw reads (dbGaP phs002407) are not used.

| File | sha256 |
|---|---|
| `GSE178341_crc10x_full_c295v4_submit.h5` (1,203,550,558 bytes) | `f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670` |
| `GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz` | `42ed840feb62aec786f092a18ea8c0cb88d07f919df207cc8f12dc2b08e96699` |
| `GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz` | `0452149b0bcff9660951de87913b82a4d72784c08d0c782e49395f756a4b9853` |

The metadata hashes equal those read in E0. The count file holds all 62 patients because GEO offers no
per-patient file; only cells of the 19 registered donors are read out of it (s.2.1) and nothing else in it
is used. The raw files stay on thinkstation1 (`~/workspace/pelka_crc_e2/data/raw`) and are not committed.

### 2.1 Selection (`scripts/build_cohort_pelka.py`)

- **Cells:** author annotation `clMidwayPr` in {TCD4, TCD8}; `PROCESSING_TYPE == "unsorted"`;
  `SPECIMEN_TYPE` T (origin `tumor_primary`) or N (origin `normal_adjacent`).
- **Donors:** at least 100 such cells in both tissues. The recomputed set must equal the registered 19, or
  the script refuses: **C107, C110, C111, C123, C125, C126, C129, C130, C132, C134, C135, C137, C140, C142,
  C143, C155, C157, C162, C170.** This is the E0 unsorted-only count; restricting both tissues to unsorted
  cells removes the tumour-only CD45+ sorting imbalance E0 found in 11 of the 25 all-processing donors.
- **Draw:** unchanged from the balanced-donor study: seed 20260924; per (donor, origin) the sorted cell IDs
  are permuted by `default_rng(int(sha256(f"{seed}|{donor}|{origin}")[:16], 16))`; rank < 100 is the
  analysis cohort (exactly 100 cells per donor per tissue, 3,800 cells), rank < 300 the training
  pool (9,039 cells). A unit test checks that the ranks equal those of `balanced_donor_luad/scripts/build_cohort.py`.
- **Genes:** version and suffix are stripped from the feature IDs (`ENSG00000243485.5_4` becomes
  `ENSG00000243485`); the 28 pseudoautosomal `_PAR_Y` copies collapse onto their X-chromosome IDs by
  summation (the only collisions). The 35 non-gene features in the file (4 cell hashtags and
  31 antibody-derived tags) are dropped. Counts must be non-negative integers; the collapse must preserve
  the total.
- **Doublets and QC:** the GEO object contains only cells that passed the authors' QC; there is no separate
  doublet column, and no further cell filter is applied.

Result (`cohort/COHORT_DEFINITION.json`; 18,055 eligible cells; 43,050 genes after the collapse):

| Donor | MMR | 10x chemistry | Unsorted CD4+CD8 T cells, tumour | Normal | Held out in fold |
|---|---|---|---|---|---|
| C107 | MMRp | 3' v2 | 664 | 127 | 2 |
| C110 | MMRd | 3' v2 | 279 | 759 | 1 |
| C111 | MMRd | 3' v2 | 226 | 219 | 0 |
| C123 | MMRd | 3' v2 | 1,577 | 257 | 3 |
| C125 | MMRp | 3' v2 | 235 | 149 | 0 |
| C126 | MMRp | 3' v2 | 172 | 159 | 4 |
| C129 | MMRp | 3' v2 | 719 | 411 | 1 |
| C130 | MMRd | 3' v2 | 600 | 161 | 2 |
| C132 | MMRd | 3' v2 | 2,220 | 731 | 1 |
| C134 | MMRp | 3' v2 | 137 | 248 | 2 |
| C135 | MMRp | 3' v2 | 659 | 104 | 4 |
| C137 | MMRd | 3' v2 | 254 | 247 | 0 |
| C140 | MMRp | 3' v2 | 325 | 124 | 0 |
| C142 | MMRd | 3' v2 | 450 | 127 | 3 |
| C143 | MMRd | 3' v2 | 746 | 191 | 4 |
| C155 | MMRp | 3' v3 | 950 | 161 | 3 |
| C157 | MMRp | 3' v3 | 205 | 157 | 3 |
| C162 | MMRp | 3' v2 and v3 (both tissues) | 641 | 1,457 | 1 |
| C170 | MMRd | 3' v3 | 691 | 516 | 2 |

MMRd 9, MMRp 10. Chemistry is constant within every donor except C162, whose tumour and normal samples
each include both 3' v2 and 3' v3 libraries (in the drawn pool: tumour 115 v2 / 185 v3, normal 133 v2 /
167 v3), so chemistry never differs between a donor's two tissues. CD4:CD8 in the pool: tumour 3,215 :
1,893, normal 2,818 : 1,113.

### 2.2 Tokenisation (`scripts/tokenize_cohort.py`)

A copy of the balanced-donor script with only the file prefix changed: Geneformer V2
`TranscriptomeTokenizer` (model input 4096, special tokens on), loom input path, and the same QC rows
(`qc/QC_TOKENIZATION.md`). All QC rows passed: 9,039 cells in and out, 3,800 analysis cells flagged, median 769 tokens per cell (LUAD 755), minimum 192; 6 cells reach the 4,096-token input limit and are truncated by the tokenizer (2 of them analysis cells). 45.5% of the 43,050 input genes are in the V2 dictionary (report only; the unmapped IDs are mostly non-coding). Dataset sha256 in `provenance/tokenized_dataset_sha256.txt` (arrow `9289e6cf…`).

## 3. Compute environment (thinkstation1 only)

| Item | Value |
|---|---|
| Host, GPU | thinkstation1, NVIDIA GB10 (the only host; no cross-host split, so no host effect enters any contrast) |
| Python environment | `~/workspace/geneformer-uv-starter/sclc_analysis/.venv` (Python 3.12, torch 2.13.0+cu130, datasets 5.0.1, pandas 2.3.3, which satisfies the pandas < 3 rule) |
| Geneformer | `~/workspace/geneformer-uv-starter/Geneformer` at `f45a6c7de57f` with the bf16 export patch (sha256 of `git diff` `d8e43471…0215281f`, as recorded in `balanced_donor_luad/provenance/ENVIRONMENTS.md`) |
| Base model | Geneformer-V2-316M, `model.safetensors` sha256 `965ceccea81953d362081ef3843560a0e4fef88d396c28017881f1e94b1246f3` |
| Token dictionary | `token_dictionary_gc104M.pkl`, sha256 `67c445f4385127adfc48dcc072320cd65d6822829bf27dd38070e6e787bc597f` |
| Precision | bf16 (dtype_cast), deterministic kernels off, as in the balanced-donor study |
| Package set | sorted `uv pip freeze` without the editable Geneformer line, recorded 2026-10-01 before any GPU call: 208 lines, sha256 `acf2c2d83cfa3cb6fe85cf22ea3a32d2021898b34f789e257a3b42af800a4c05` (`provenance/env_freeze_ts1_sorted.txt`). It differs from the LUAD-run freeze (`353d6596cbcab16c…`, `balanced_donor_luad/provenance/ENVIRONMENTS.md`) in one line only: pandas 3.0.5 (LUAD) versus 2.3.3 (now, Oscar's pandas < 3 downgrade after A6). Replacing that single line in the current freeze reproduces `353d6596cbcab16c…` exactly. datasets 5.0.1, pyarrow 25.0.1, torch 2.13.0, transformers 4.46.0 and numpy 2.5.2 are unchanged. |

Before each GPU step the driver re-checks the pinned files (`registration/pins_prep.sha256`,
`registration/pins_isp.sha256`) and refuses on any mismatch. The GPU is checked free before launch; no
other user's job is touched.

## 4. Code: what is reused unchanged and what is new

**Reused unchanged from `balanced_donor_luad/scripts/`** (the files on `main` at `6d4b2be`, pinned):
`finetune_fold.py`, `classifier_gate.py`, `goal_embeddings.py`, `pre_isp_gates.py` (no-op gate),
`run_isp.py`, `model_cache.py`, `inproc_map.py`, `noop_spotcheck.py`, `isp_order.py`,
`build_ovx_index.py`, `analyse.py`, `stats_core.py`, `null_analysis.py`, and
`sclc_validation/bf16_bench/{dtype_cast,eval_tie_fix}.py`,
`sclc_validation/perturbation_workflow/s100_isp/matched_controls.py`.

**New or copied with changes, in `pelka_crc_e2/scripts/`:**

| Script | Relation to the lung pipeline |
|---|---|
| `build_cohort_pelka.py` | new: the Pelka reader (s.2.1); the draw function is identical to `build_cohort.py` |
| `tokenize_cohort.py` | copy; file prefix only |
| `make_split.py` | copy; `N_EVAL = 2` instead of 4 (s.5.1) |
| `eligibility_controls.py` | copy; Panel B only, one host; rules, tolerances, seed and `d_min` unchanged |
| `select_null_e2.py` | combines `select_null_genes.py` and `select_estimable_prefix.py`; new seed and exclusions (s.5.3) |
| `make_run_order.py` | new: the frozen ISP order and its cost prediction |
| `build_design_e2.py` | copy of `build_design.py` without Panel A and the ambient table; not-run set derived, not hard-coded |
| `compare_panel_b.py` | new: the H2c primary reading (s.6.3) |
| `run_e2_gpu_prep.sh`, `run_e2_isp.sh`, `run_e2_analysis.sh` | drivers modelled on `run_phase4.sh` and `run_phase8_null.sh` |

`pelka_crc_e2/tests/test_e2.py` covers the new logic with synthetic data (16 tests pass on
thinkstation1; five deliberate code mutations were tried and four were caught; the fifth, dropping the
explicit zero check in `compare_panel_b.py`, is equivalent because no LUAD reference median is zero).

Pinned hashes at registration: `registration/pins_prep.sha256` (16 files: the tokenised dataset, base weights, token dictionary, folds, cohort table, and every script the GPU prep steps execute; file sha256 `3328afe4…`) and `registration/pins_isp.sha256` (31 files: the same shared inputs plus every script, rules file, gene list, frozen count table and LUAD reference file the ISP and analysis steps read; file sha256 `2cc4dc90…`). Paths are absolute on thinkstation1. Outputs of GPU steps that cannot exist before the run (fold models, `donor_manifest.json`) are hashed at ISP launch into `post_prep_inputs.sha256`.

## 5. Procedure

Each numbered step runs only after the previous one finished without a stop (s.9).

### 5.1 Folds (CPU, done, committed)

Five-fold donor cross-fitting, seed 20260924, unit = donor, both tissues of a donor always together. With
19 donors the folds hold out 4, 4, 4, 4 and 3 (fold 0: C111, C125, C137, C140; fold 1: C110, C129, C132, C162; fold 2: C107, C130, C134, C170; fold 3: C123, C142, C155, C157; fold 4: C126, C135, C143) donors. Within each fold, 2 of the training donors are set
aside as eval (reported, never used for selection). The balanced-donor study used 4; with 19 donors that
would have left 11 donors to train each fold, so 2 is used and 13 train. `cohort/folds.json`.

### 5.2 Panel eligibility and matched controls (CPU, done, committed)

Panel B is the balanced-donor file `balanced_donor_luad/registration/panel_B.json` unchanged: 39 listed,
36 perturbable (TRAC, TRBC1, TRBC2 are NOT_RUN by construction). Estimable (g, d): gene token in at least 10
of donor d's 100 tumour analysis cells. Eligible: estimable in at least `d_min = 11` donors, the balanced-donor
Panel B value; with an exact two-sided Wilcoxon test and Holm over the 36-gene family, 11 is the smallest
number of donors at which a gene can reach adjusted p <= 0.05 (2/2^11 x 36 = 0.035), so it is kept. Strata
and 20 matched controls per stratum by the registered procedure (detection within 0.5 log2, rank percentile
within 5 points, seed 20260924, panel genes and the ambient anchors excluded from the candidate pool). A
stratum with fewer than 20 candidates is NOT_ESTIMABLE_CONTROLS and its members are not run.

Result: 32 of the 36 perturbable genes are eligible; CCR7 (8 estimable donors), THEMIS (8), IKZF2 (4) and
LEF1 (3) are not, and become NOT_RUN rows. The 32 fall into 16 strata. Two strata have fewer than 20
candidate controls because their members are among the most highly detected genes in these cells: B13
(CD7, CD69) and B15 (CD2, CD3D). These four genes are NOT_ESTIMABLE_CONTROLS and are not run. That leaves
28 panel genes to run, with 241 unique control genes across the 14 remaining strata
(`controls/pre_gpu_eligibility.csv`, `pre_gpu_strata.json`, `pre_gpu_summary.json`). Two of the 13 LUAD
reference genes (CD3D, CD7) therefore have no controls in E2 and one (CCR7) is ineligible; 10 remain
testable for H2c.

Panel A (the 15 July LUAD screen hits) is not part of E2: it was a list of lung tumour-associated genes,
mostly epithelial, and 13 of 15 were not testable in LUAD T cells. The ambient classifier is therefore not
refitted either: ambient labels enter only Panel A statuses, and every Panel B row carries the fixed label
CURATED_T_CELL.

### 5.3 Null genes (CPU, done, committed)

Sampling frame: every gene in the pinned token dictionary except the four special tokens, the 39 Panel B
genes and every E2 control and stratum member. Seeded shuffle, seed 20261001 (the balanced-donor draw used
20260929; a new seed gives a fresh draw). The draw order is walked and the first 100 ESTIMABLE genes are
frozen, estimable meaning token in at least 10 of a donor's 100 tumour analysis cells in at least 10 donors
(the Amendment 5 numbers, unchanged; with 19 donors this asks for detection in more than half the cohort,
where LUAD asked for 10 of 43). Every walked gene's per-donor token-positive count is frozen with it: it is
the gate-identity reference `null_analysis.py` checks every ISP marker against.

Result: the 100th estimable gene sits at position 787 of the draw (687 walked genes were not estimable); frame 19,994 genes after excluding 280 (39 Panel B, 241 controls; the 32 stratum members are Panel B genes). Files: `null/null_genes_100_estimable.json`, `null/frozen_100_estimable.json`.

### 5.4 GPU steps, thinkstation1

1. **Fine-tunes** (`run_e2_gpu_prep.sh` -> `finetune_fold.py`): per fold, recipe unchanged: 1 epoch, lr 5e-5,
   batch 8, 6 frozen layers, seed 43, label `origin`, bf16, eval_tie_fix; training on the training donors'
   pool cells; held-out donors scored on their analysis cells (exactly 100 + 100).
2. **Classifier gate** (`classifier_gate.py`): s.6.1. A FAIL stops E2 here. The gate JSON is committed as
   `pelka_crc_e2/phase4_results/classifier_gate.json` before the ISP starts (`build_design_e2.py` reads it).
3. **Goal centroids and ISP inputs** (`goal_embeddings.py`): per donor, the exact-mean CLS embedding of
   the donor's own normal-tissue pool cells (the registered goal), its tumour centroid, and the fold's
   training-donor normal centroid (S1 sensitivity), each with the donor's own fold model; ISP input = the
   donor's 100 tumour analysis cells.
4. **No-op gate** (`pre_isp_gates.py --mode noop`): two forward passes, no edit, on the first donor with
   B2M: every per-cell shift must be exactly 0. A FAIL stops E2.
5. **ISP** (`run_e2_isp.sh` -> `run_isp.py`): every gene in `controls/isp_run_order.json`, in that order,
   deletion and overexpression, all 19 donors, each donor with its own fold model and its own goal. Order:
   the 100 null genes (H2b first), then the eligible Panel B genes, then the controls.
6. **No-op spot checks** (`noop_spotcheck.py`): every 20th gene of the null list and of the panel+control
   list, first donor, two independent runs; the shift difference must be exactly 0.

**GPU budget.** Predicted ISP: 12.3 GPU-hours (369 genes, of which 100 null, 28 panel and 241 controls, x 19
donors x 2 operations; the balanced-donor Phase 6 per-call fit applied to the actual token-positive counts;
`controls/isp_run_order.json`). Fine-tunes about 1.5 GPU-hours (13 training donors per fold, about 6,200
cells, at the measured 8.19 cells per second, plus scoring); goals under 0.5. The authorised ceiling is 52 GPU-hours in total. The ISP
driver enforces a hard ceiling equal to 52 hours minus the GPU time already used by steps 1 to 4 (measured,
recorded at launch), summed from the markers' own `seconds` field. If a stop is reached, the run halts, no
partial arm is analysed, and god is told.

## 6. Analysis (CPU, after the ISP completes; `run_e2_analysis.sh`)

### 6.1 H2a: classifier gate (unchanged)

PASS if pooled held-out balanced accuracy >= 0.60 AND the exact two-sided sign test over donors of
"per-donor balanced accuracy > 0.5" has p <= 0.05 with more donors above than below (donors at exactly 0.5
dropped). With 19 donors, the sign test reaches p <= 0.05 at 15 of 19 above chance (p = 0.019; 14 of 19
gives 0.064). Per-fold values and the nondeterminism fields are reported; the band inside
`classifier_gate.py` was measured on LUAD and is reported as is, never used to change PASS.

Readings: PASS, the tumour/normal T-cell signal transfers to colon under the same recipe. FAIL, E2 stops
before any perturbation and is reported as "classifier gate failed" with the numbers; H2b and H2c are then
not run.

### 6.2 H2b: null anti-correlation, primary (`null_analysis.py`, unchanged)

Per null gene and operation: the donor-level mean per-cell cosine shift toward the donor's own normal
centroid, over token-positive cells (overexpress positions from `build_ovx_index.py`); the gene's value is
the median over donors estimable for that gene (>= 10 token-positive cells), requiring >= 10 such donors.
No control adjustment (Amendment 4 s.4.4). Primary test: Spearman rho between the 100 deletion medians and
the 100 overexpression medians, one-sided permutation test (H1: rho < 0), 100,000 permutations, seed
20260929 inside the script. Stability: leave-one-gene-out (each with 2,000 permutations) and a 10,000-draw
gene bootstrap. Gate identity: every marker's token-positive count must equal the frozen count (s.5.3),
else the analysis halts. No-op: any failed or missing spot check forces `no_op_failed`.

Statuses (ISP-STD-1 A.6, as in the script): `positive` (rho < 0, p <= 0.05, every leave-one-out also
p <= 0.05 with rho < 0, bootstrap stable fraction >= 0.95); `control_draw_sensitive_open` (rho < 0 and
p <= 0.05 but not stable); `negative` (not met; reported with direction, p and n, never "no effect");
`opposite_direction` (rho > 0 with upper-tail p <= 0.05); `not_estimable`, `no_op_failed`,
`stopped_not_analysed`.

Validity check (unchanged): the same pipeline on every E2 control gene; it passes if rho < 0 (sign only;
magnitudes are not compared with any LUAD figure). If it fails, the primary is reported but marked not
trusted.

### 6.3 H2c: Panel B pattern

Per-gene statuses: `analyse.py` unchanged with the balanced-donor rules file
(`registration/phase7_rules.json`: control_min_cells 10, min_controls_per_donor 10, Holm over the full
36-gene family): control-adjusted donor values, exact two-sided Wilcoxon, Holm, dose concordance, the
T_CELL_SIGNAL_TOWARD / AWAY / DELETION_ONLY / OPEN / NOT_RUN / NOT_ESTIMABLE_CONTROLS statuses, the S3
fold downgrade, S1 (global goal) and S4 (cell-weighted) sensitivities. Host drift is false (one host).

Primary H2c reading (`compare_panel_b.py`): the LUAD reference set is the 13 Panel B genes with status
T_CELL_SIGNAL_TOWARD or T_CELL_SIGNAL_AWAY in `balanced_donor_luad/phase7_results/outcome_rows.json`
(sha256 `083b1ef0e83682c093486b58224a2081f8b41641b51430b45991ff79582ad471`): CCR7, CD247, CD27, CD3D, CD3G, CD7, CD8A, GZMA, ICOS, LAT, LCK (TOWARD) and
CTSW, ITK (AWAY). Among those tested in E2 (a control-adjusted deletion median exists), count the genes
whose E2 deletion median has the same sign as in LUAD (zero counts as not agreeing). Exact one-sided
binomial test against 1/2. With the expected n = 10 (s.5.2), `pattern_holds` needs at least 9 of 10 genes
with the LUAD sign (P(X >= 9) = 11/1024 = 0.0107; P(X >= 8) = 56/1024 = 0.0547). The bar is strict and is
stated here before any number exists.

| Reading | Rule |
|---|---|
| pattern_holds | n_tested >= 5 and p <= 0.05 |
| pattern_not_replicated | n_tested >= 5 and p > 0.05, reported with a/n and p |
| not_testable | n_tested < 5 |

Descriptive only: the same count on overexpression medians, the full LUAD-versus-E2 status table, Spearman
rho between LUAD and E2 deletion medians over genes tested in both, and per-gene medians within MMRd and
MMRp donors.

## 7. How the results will be read

| H2b in colon | Reading |
|---|---|
| positive | The anti-correlation is present in lung and colon. It is then more likely a property of Geneformer and this perturbation design than of lung T cells; model and design cannot be separated here. |
| negative or opposite_direction | The anti-correlation did not appear in colon. This cannot be attributed to tissue alone: study, dissociation, chemistry, annotation and a new fine-tune all differ. |
| control_draw_sensitive_open | Direction as in lung but not stable; reported as open. |

H2c `pattern_holds` says the direction of the LUAD reference genes carries over to colon T cells under
this design; it does not make any gene a regulator of T-cell state. H2a is a precondition, not a finding
about any gene.

## 8. Sensitivities and descriptive outputs (never change a primary reading)

S1 global goal and S4 cell-weighted (from `analyse.py`); MMRd versus MMRp medians (s.6.3); per-fold
balanced accuracies. The 25-donor all-processing cohort from E0 is **not run**; it remains a possible later
sensitivity for the processing imbalance and would need its own registration.

## 9. Stop conditions (each stops E2 and is reported to god)

- Any pinned input hash mismatch before a GPU step.
- Classifier gate FAIL (after step 2): no goals, no ISP.
- No-op gate FAIL (step 4).
- Projected or measured GPU use above 52 hours in total (the ISP ceiling).
- `build_ovx_index.py` count-check failure, a `null_analysis.py` gate-identity mismatch, or an
  `IncompleteRun` in `analyse.py`.
- Any failure this registration did not foresee: stop, record it, tell god, and change nothing until god
  answers.

## 10. What E2 cannot show

- A different result in colon cannot be attributed to tissue alone: Pelka is a different study, with its
  own dissociation protocol, 10x chemistry version, T-cell annotation and sorting history.
- E2 trains new fold models, so a colon result is not a test of the models that produced R1 to R3.
- 19 donors from one study give per-gene power well below the LUAD study's 43 donors from five studies;
  a gene that is OPEN in colon is not evidence against its LUAD status.
- The within-donor ambient difference between tumour and normal tissue remains, as in LUAD; no ambient
  classifier is fitted for colon.
- The E2 null genes are a more highly detected population than the LUAD null genes. The unchanged
  estimable rule (>= 10 token-positive cells in >= 10 donors) asks for detection in more than half of 19
  donors, against 10 of 43 in LUAD; E2 reached its 100th estimable gene at draw position 787, LUAD A5 at
  984. An E2 H2b result that differs from LUAD may partly reflect this difference in the gene
  population, not only tissue.
- Nothing here establishes that any gene has a causal role in T-cell state.

## 11. Differences from the design document (s.6 of the alt-dataset design)

- The design expected about 36 donors with normal tissue; E0 found 25 qualifying and the human chose the
  19 unsorted-only donors.
- Panel A and the colon ambient refit are dropped (s.5.2); the design's "36 panel genes" is Panel B.
- The control count is whatever the registered procedure yields (at most 20 per stratum), not a chosen
  120 or 318.
- The 5-fold design is kept; eval donors per fold 2 instead of 4 (s.5.1).
- One host instead of two; the Amendment 3h GPU identity check of the overexpress position replay is not
  re-run. The code is identical and the package set has one change (pandas 3.0.5 -> 2.3.3, s.3); the
  libraries 3h depends on (datasets, pyarrow, torch) and transformers are unchanged, so the identity check
  is not re-run. The two count checks in `build_ovx_index.py` still run on every call.

## 12. Deviation log (dated entries only; empty at registration)

(none)
