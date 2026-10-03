# Criteria inventory: the protocol's and ours

Written 2026-10-03 as the frozen list that the cross-evaluation in
[`EXECUTION_PLAN.md`](EXECUTION_PLAN.md) tests. Each criterion has an ID, a one-line rule and the
file that defines it. Nothing here is new; it collects rules from two sources:

- **P**, the Geneformer protocol: Zhang, Venkatesh & Theodoris, *Nature Protocols* (2026),
  doi:10.1038/s41596-026-01364-8. Step numbers refer to its Procedure.
- **O**, this repository's criteria, built between July and October 2026. Most are written into the
  balanced-donor LUAD registration
  ([`PHASE5_ISP_REGISTRATION.md`](../balanced_donor_luad/registration/PHASE5_ISP_REGISTRATION.md),
  below "LUAD reg"), the E2 registration
  ([`E2_REGISTRATION.md`](../pelka_crc_e2/registration/E2_REGISTRATION.md)) and the standard
  ISP-STD-1 v1.1 (`hive/standards/isp-outcome-criteria.md`). ISP-STD-1 is not in this repository; it is
  copied in and pinned by sha256 in Phase 0 of the plan.

## P: the protocol's criteria

| ID | Rule | Source |
|---|---|---|
| P1 | Train, validation and test sets hold different patients (`attr_to_split="individual"`), with confounders such as disease, age and sex balanced across them (`attr_to_balance`) | Experimental design; Step 19 |
| P2 | Hyperparameters are tuned on the target data (Ray Tune, or Optuna for multi-task) | Steps 19, 29 |
| P3 | Classification is reported on held-out patients as a confusion matrix, accuracy and macro F1 | Steps 20–21 |
| P4 | Each gene's ISP effect is `Shift_to_goal_end`, the cosine shift of the start state toward the goal, compared with the shifts from perturbing random genes by Wilcoxon test, BH-corrected; `Sig` = 1 if `Goal_end_FDR` < 0.05 | Step 26, `goal_state_shift` |
| P5 | Candidates should shift toward the goal and away from an alternative state (`Shift_to_alt_end`, `Alt_end_FDR`) | Steps 23, 26; Fig. 2b |
| P6 | `N_Detections`, the number of cells in which the gene was detected, is reported with each gene | Step 26 |
| P7 | Quantized and full-precision perturbation shifts agree (reported Pearson r > 0.99) | Steps 33–34; Fig. 2d |
| P8 | The embedding layer and cell filter are the same for the state embeddings and the perturbation (`emb_layer=0` for a fine-tuned model) | Steps 22–24 |
| P9 | Embeddings are inspected (UMAP) for separation of the states of interest | Steps 15–16, 22 |
| P10 | A top-ranked target is validated experimentally (GSN knockout in an iPSC model of DCM) | Anticipated results |

## O: this repository's criteria

| ID | Rule | Source |
|---|---|---|
| O1 | Classes are audited for study and donor confounding before training: no class may be carried by studies or donors that contribute no other class | `current_workflow/METHODS.md`, "Class-construction defects" (2026-09-24) |
| O2 | Donor-paired design where the biology allows (tumour and normal from the same donor), with equal analysis cells per donor (100 per tissue) so each donor weighs the same | LUAD reg s.2–3; E2 reg s.2 |
| O3 | Donor-disjoint k-fold cross-fitting; both tissues of a donor in one partition; eval donors never used for selection; every donor scored by a model that never saw it | LUAD reg s.2; E2 reg s.5 |
| O4 | Classifier gate: pooled held-out balanced accuracy ≥ 0.60 and an exact two-sided sign test over donors (p ≤ 0.05); a failure stops ISP and is a reportable finding | LUAD reg s.2, Amendment 3d |
| O5 | Nondeterminism band: results within ±0.045 of a threshold, measured from same-seed repeats, are reported as within run-to-run noise | LUAD reg Amendment 3e |
| O6 | No-op gate: a perturbation with no edit must give a shift of exactly 0, before the screen and in spot checks during it | LUAD reg s.3, s.8; E2 reg |
| O7 | Per-donor value from the raw per-cell output against the same donor's goal state, never from the cell-pooled `Shift_to_goal_end` | LUAD reg s.3 |
| O8 | Matched controls: each gene is adjusted against 20 control genes matched on detection and token rank; fewer than 20 gives `NOT_ESTIMABLE_CONTROLS` | LUAD reg s.3, Amendment 3c |
| O9 | Eligibility: a gene is tested only if estimable in at least d_min donors, with d_min set by the smallest attainable p under Holm | LUAD reg s.4 |
| O10 | Exact two-sided Wilcoxon over donors, Holm within a registered family | LUAD reg s.5a |
| O11 | Both deletion and overexpression; a positive status needs both arms significant with opposite signs (`COHERENT`); same-sign pairs are `DOSE_INCOHERENT` | LUAD reg s.5b–5c |
| O12 | Ambient RNA: a non-circular, leave-one-anchor-out classifier (AUC ≥ 0.8) flags genes that may come from ambient RNA; a flagged or undetermined gene can never be a plain positive | LUAD reg s.0a, s.1c |
| O13 | Random-gene null: across random eligible genes, the rank correlation of deletion and overexpression shifts, by permutation test, with leave-one-out and a bootstrap stability bar of 95% | LUAD reg Amendments 4–5; E2 reg s.6.2 |
| O14 | Fold heterogeneity: a positive whose deletion median flips sign in more than one fold is downgraded to `OPEN` | LUAD reg s.7 (S3) |
| O15 | Cross-host equivalence: one probe gene on both GPU hosts, per-cell max abs difference ≤ 1e-3 and Spearman rho ≥ 0.999, before and after the run | LUAD reg Amendment 3, s.3.5 |
| O16 | Precision canary: bf16 against fp32 on the same model, Spearman rho ≥ 0.95 and top-20 overlap ≥ 19 of 20; model size is a separate factor (104M against 316M: rho 0.489, 9 of 20) | `sclc_validation/bf16_bench/RESULTS_BF16_316M.md` |
| O17 | Every gene keeps a row with a status from a fixed vocabulary; "no effect" is never written | LUAD reg s.5c, s.4.8 |
| O18 | Pre-registration before any GPU output, a gatekeeper sign-off, sha256 pins, and dated append-only amendments | LUAD reg throughout; ISP-STD-1 A.5, E |
| O19 | Replication in an independent cohort, run unchanged | E2 reg s.1 |

## Where the two sets meet

Most protocol criteria have a counterpart in ours. The cross-evaluation compares the counterparts and also
asks what each unmatched criterion catches.

| Protocol | Ours | How they differ |
|---|---|---|
| P1 | O1, O2, O3 | Both split by patient. O1 also audits that the class labels are not carried by study or donor, and O2 weights donors equally |
| P3 | O4, O5 | Macro F1 pooled over cells, against balanced accuracy pooled and per donor with a sign test and a noise band |
| P4 | O7, O8, O9, O10 | Cell-pooled Wilcoxon against random genes, against a donor-level Wilcoxon on control-adjusted values |
| P5 | O11 | Toward the goal and away from an alternative state, against toward the goal under both deletion and overexpression with opposite signs |
| P7 | O15, O16 | Quantized against full precision, against cross-host, bf16 against fp32, and model size |
| P10 | O19 | Wet-lab validation of one target, against replication in another cohort; neither replaces the other |
| none | O6, O12, O13, O14, O17, O18 | No protocol counterpart |
| P2, P6, P8, P9 | none or partial | No counterpart in ours: we fix the recipe across cohorts instead of tuning (P2), and report detection through eligibility (P6) |
