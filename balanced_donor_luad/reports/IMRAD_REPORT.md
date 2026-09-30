# Balanced-Donor LUAD Tumor-vs-Normal T-Cell In Silico Perturbation Study

**Status: DRAFT — Results and Discussion are placeholders pending the independent gatekeeper's
post-run check of the null-study analysis outputs. Nothing in this document states or implies any
result of that pending check. No number below is a scientific outcome (a correlation, a p-value, or
an ISP-STD-1 status) unless it was already closed, published, and merged to `main` before this null
study began.**

*Repository: `Kays3/geneformer-lung-tcell`. This document lives on branch
`report/imrad-balanced-donor-null` and is built from the registration documents, commit history, and
run artifacts cited inline; every design number below is cited to a file and a commit so a reader can
verify it independently of this prose.*

---

## Abstract

*(Placeholder. Written after Results/Discussion are filled in.)*

---

## 1. Introduction

### 1.1 Background

Geneformer is a transformer model pretrained on single-cell transcriptomes; in silico perturbation
(ISP) uses it to estimate how removing (`delete`) or forcing (`overexpress`) a gene's rank in a cell's
input shifts that cell's embedding toward or away from a target biological state. An earlier screen
over the Lung Cancer Atlas (LuCA) core atlas ("the July analysis") ran an all-gene deletion screen
against a three-class classifier (LUAD tumor / LUSC tumor / normal) built from CELLxGENE's own disease
labels. On inspection (2026-09-24, read-only), that "normal" class drew from six non-cancer studies
that contributed zero tumor T cells — the class was study-confounded — and the "LUAD" class was
selected by disease label alone, a label that in the core atlas spans tumor, adjacent-normal,
metastasis, and non-cancer tissue of LUAD patients, so the class likely mixed sites. Twenty of the 43
donors used in the present study overlap the core atlas's LUAD donor pool the July screen drew from.

### 1.2 Motivation for the balanced-donor redesign

This study asks a narrower, better-controlled version of the July question: for each donor, does
Geneformer read a consistent shift between that donor's own tumor and adjacent-normal tissue, with
every other axis (study, chemistry, donor identity, cell count) held equal between the two arms of the
comparison? Pairing within donor cancels study and chemistry; a fixed per-donor, per-tissue cell cap
cancels the raw 28–31× cell-count imbalance between donors seen before capping (Provenance table row
P1).

### 1.3 Study questions

1. **Panel A (replication).** Do the 15 top LUAD→normal deletion hits from the July screen replicate
   under the donor-paired design, in both the deletion and overexpression directions?
2. **Panel B (T-cell-intrinsic signal).** Does the model read any T-cell-intrinsic signal on a
   curated panel of 36 canonical T-cell genes the July screen never ranked highly?
3. **Null study (false-positive-rate calibration).** Ranging over eligible random genes detectable in
   LUAD tumor T cells (never "genome-wide" — the population is bounded by detectability in this
   cohort, ISP-STD-1 B.9), does the same delete/overexpress raw-shift pipeline show a systematic
   negative concordance even for genes with no a priori relationship to the tumor/normal axis? A
   negative answer would suggest Panel A/B's own concordance findings reflect gene-specific biology,
   not a pipeline artifact; a positive answer would narrow what any concordance finding in this study
   licenses.

### 1.4 Governing standard

All work after 2026-09-28 (13:35 JST) is governed by ISP-STD-1 v1.1, the floor standard the human
adopted for every in-silico-perturbation study (`hive/standards/isp-outcome-criteria.md`). Its fixed
outcome vocabulary — `positive`, `negative`, `opposite_direction`, `not_estimable`,
`control_draw_sensitive_open`, `no_op_failed`, `not_run_by_design`, `stopped_not_analysed` — is used
throughout this report; "no effect" is never used in place of `negative`, and "genome-wide" is never
used in place of "eligible random genes detectable in LUAD tumor T cells." Work completed before
2026-09-28 (Phases 4–7 below) was not retroactively re-gated against ISP-STD-1, consistent with the
standard's own adoption note, but already satisfied most of its MUSTs by construction (matched
controls, a no-op gate, donor-level statistics, hash-pinned runner/weights/data, ambient flagging via
leave-one-anchor-out).

---

## 2. Methods

### 2.1 Data source and cohort (Phase 0–3)

**Source.** The single-cell Lung Cancer Atlas, extended atlas (Salcher et al., *Cancer Cell* 2022).
CELLxGENE dataset version `33165751-ae33-4a65-94ee-52d9bc38f97e`; 1,283,972 cells × 17,764 genes, raw
integer counts. Downloaded 2026-09-24 to a dedicated data root and verified against the published
object by S3 multipart ETag match and a full-file sha256 (`f0f7f43413bfe088281ac68487ad4bd73647d9057280da202a82fabce9937cb9`)
— Provenance table row P2.

**Selection rule** (mechanical, registered before any screening; `provenance/PHASE0_SELECTION_RULE.md`,
committed `3af3ca98`, 2026-09-24 18:19 JST, amended `8ea11784` same day to require a donor-level primary
analysis, ≥6 paired donors held out, and ≥12 eligible genes):

| Filter | Value |
|---|---|
| Disease | Lung adenocarcinoma (LUAD) |
| Cells | CD4 and CD8 T cells (`cell_type_major`), singlets |
| Chemistry | 10x only (3′ v2; one study also 5′ v1) |
| Tissue sites | `tumor_primary` vs `normal_adjacent` only — no metastases, effusions, blood, or nodes |
| Donor rule | ≥100 T cells in **both** tissues |

**Resulting cohort: 43 donors, both tissues each, 5 studies** (Leader_Merad 2021: 23; Kim_Lee 2020: 10;
He_Fan 2021: 5; Lambrechts_Thienpont 2018: 3; Laughney_Massague 2020: 2). 132,619 T cells were eligible
before capping; raw donor-level cell counts ranged 237–6,764 (tumor) and 345–10,384 (adjacent normal),
a 28–31× donor imbalance. A fixed-seed draw (seed 20260924) then fixes, per donor per tissue, exactly
100 analysis cells (8,600 total, max/min = 1.0) and up to 300 pool cells (25,700 total) for training.
Cells are tokenized for Geneformer V2 (median 755 tokens/cell; max observed 2,656 < the 4,096 context
length, so no truncation occurs). Committed on `analysis/balanced-donor-luad` @ `3af3ca98` (Provenance
row P3).

**Excluded, and why.** LUSC tumor tissue (only 9 qualifying donors on 10x within LuCA, below the
registered ≥12 floor; a wider public-data survey was scoped under a draft Amendment 3/L0–L2 process but
superseded when the human settled on the registered 2-class design, `e1a3ed94`, 2026-09-24 19:51 JST);
17 LUAD donors with only one qualifying tissue; and a structural selection bias toward T-cell-rich
samples, since the donor rule requires ≥100 T cells in both tissues.

### 2.2 Model and fine-tuning design (Phase 4)

**Model.** Geneformer-V2-316M, bf16, per the human's directive (no 104M arm — registered in Amendment 1
as a deliberate design choice, not an oversight; any Panel A non-replication is therefore reported
with the pre-registered caveat "non-replication on 316M; model-vs-design attribution not established,"
fixed in the wording before any output existed, `f5805cbc`, 2026-09-24 18:51 JST). Calibrated
throughput (measured before any GPU spend on the real cohort, not a result of this study): 8.19 cells/s
fine-tuning; ≈20 s per gene for ISP on 300 cells (both upper bounds, deterministic kernels off) —
committed as Amendment 1, `e2e10421`, 2026-09-24 18:47 JST.

**Design.** Five-fold, donor-stratified, study-proportional cross-fitting (`scripts/make_split.py`,
seed 20260924): every donor is held out in exactly one fold; both tissues of a donor always share a
partition; recipe is 1 epoch, lr 5e-5, batch 8, freeze 6 layers, seed 43, label `origin`.

**Gate (pre-registered, evaluated pooled across folds so every donor is scored by a model that never
saw it):** held-out balanced accuracy ≥ 0.60 (chance 0.33 for the original 3-class framing; the
donor-paired design here is 2-class, chance 0.50) **and** an exact sign test across all 43 donors
(per-donor balanced accuracy > 0.5) with p ≤ 0.05. The achieved value is a Results-section quantity,
withheld here as a placeholder — see §3.

### 2.3 ISP design and panels (Phase 5 registration, Amendments 1–3h)

**Panels**, fixed before any GPU output existed: **Panel A** — the 15 top LUAD→normal deletion hits
from the July screen (registered as a replication test, ≥10 ambient anchors required per gene).
**Panel B** — 36 canonical T-cell genes (TRAC/TRBC1/TRBC2 excluded: not present in the Geneformer V2
token dictionary, reason not otherwise inferred). Registered `552f0dfa`, 2026-09-24 18:33 JST.

**Host split.** ISP calls are split **by gene** (not by donor or by operation): each Ensembl ID is
assigned to thinkstation1 or thinkstation2 alternately, and *both* operations for that gene run on the
assigned host — registered as Amendment 3, `a5fe96a7`, 2026-09-24 19:24 JST. A cross-host equivalence
gate requires per-cell shifts on the same (gene, donor, operation) to be bit-identical within bf16
tolerance (max |Δ| ≤ 1e-3, Spearman ρ ≥ 0.999) before any host-split data is pooled.

**Matched controls.** Six shared control strata of 20 genes each (secondary use only), assembled via
`s100_isp/matched_controls.py` merged in from a sibling workstream, seed 20260924; ambient status is
kept distinct from contamination status throughout (e.g., EEF1G is ambient but not classed as
contamination) — Amendment 3c, `6ebf7fb1`, 2026-09-24 20:07 JST.

**Statistical design.** Donor-level (never cell-level) control-adjusted shift over ≈120 matched
controls per gene; exact Wilcoxon signed-rank test over donors, Holm correction per family; a gene's
delete and overexpress arms must both be significant with **opposite** signs for a `REPLICATED` /
T-cell-signal call — same-signed concordance is explicitly registered as *not* evidence of specificity,
since the existing 120/318-control null cloud itself shows the same-direction correlation. Stability
gate: leave-one-control-out plus a 10,000-draw bootstrap must each preserve significance in ≥95% of
draws for a call to stand, rather than being reported as `control_draw_sensitive_open`.

**A registered stop, disclosed as part of the record.** Phase 7 launch was **refused** on 2026-09-25
(~14:45 JST) when read-only inspection found that Geneformer's `InSilicoPerturber` silently skips its
own token filter and sorts cells by length internally, so the overexpress pickle's cell order does not
match the dataset order the registration had assumed — a genuine ambiguity in which cells the
overexpress arm's donor-level statistic was actually summarizing. A second blocker followed:
12,184 of 15,142 gene–donor pairs had tied token-positive lengths between the positive and negative
donor pool, making a tie-free identity check infeasible as first proposed. Both were resolved before
any GPU output was read: a deterministic order-reconstruction rule (verified bit-identical against a
tie-free positive control spanning k = 10 and k = 100 token-positive counts, 2,958 tie-free pairs
existed to draw from) plus a rerun-by-construction fallback for any pair the reconstruction could not
resolve. This sequence — registered stop, root-cause read before any fix, a verified reconstruction
rule, a disclosed fallback — is recorded here as part of the study's gate history (§2.8), not smoothed
over, since ISP-STD-1's reproducibility section treats a caught pipeline-integrity gap as a finding
about the pipeline, not a result about biology.

### 2.4 Perturbation execution (Phase 6)

**Scope.** 171 genes (Panel A 15 + Panel B 36 + ≈120 controls, later extended to 318 unique control
genes across 19 strata as the registration matured) × {delete, overexpress} × all 43 donors' 100
tumor analysis cells, using each donor's own fold model and a per-donor goal centroid (that donor's
held-out adjacent-normal pool mean).

**A finding, disclosed rather than silently corrected.** The runner's own docstring described the host
split as gene-disjoint ("Amendment 3.1: split by gene"), but the actual `"host"` field recorded inside
every completed marker shows each host ran a near-complete **independent copy of the full 181-gene
grid** — thinkstation1: 15,566 markers (181 × 43 × 2 exactly, all stamped `host=thinkstation1`);
thinkstation2: 14,792 (all stamped `host=thinkstation2`) — not a clean complementary split as the
docstring implied. This does not change any already-reported result (the cross-host bit-identity
equivalence gate had already verified both hosts agree on shared (gene, donor, operation) triples), but
it does change what the earlier GPU-hour close-out (ts1 14.25 + ts2 15.02 GPU-h) actually described:
two overlapping near-full runs, not two complementary halves. Lesson carried into every later stage of
this study: when a docstring or design document describes an intended data split, verify it against
the actual per-record provenance field before treating it as a cost-accounting boundary.

### 2.5 Null study design (Amendments 4–6)

**Motivation.** The delete/overexpress raw-shift pipeline already showed a negative correlation across
the study's own 120–318 matched control genes (ρ ≈ −0.62 in the stratum-adjusted RESULTS.md figure).
Because control-adjustment machinery does not apply to genes outside any registered stratum, the null
study asks the same direction-of-concordance question on genes drawn independently of any panel or
control assignment, using the **raw** (non-stratum-adjusted) shift — adjusting a gene newly drawn into
no stratum would require inventing a post-hoc matching rule, which the registration explicitly avoids.

**Gene frame and draw.** 19,902 genes = the Geneformer V2 dictionary (20,271) minus Panel A (15) minus
Panel B (36) minus every gene ever used as a matched control across the study's 19 strata (318 unique).
A single fixed-seed draw (seed 20260929) fixes the entire run order in advance, so any budget stop
always yields a clean prefix run / suffix `stopped_not_analysed` — never a gap in the middle. The draw
was later extended (byte-identical on its shared prefix at every extension) from 200 to 1,000 to 3,000
entries as the design's target gene count grew; eligibility (≥10 donors with ≥10 token-positive cells,
checked read-only, at zero GPU cost) is walked in draw order to find the Nth *estimable* gene, never
resampled.

**Design history, in order (every commit hash + JST time is the timestamp of the commit itself, not
of writing this report):**

- **Amendment 4** — a 200-**drawn**-gene design, `622efaf`, 2026-09-29 19:20 JST. Its first cost
  estimate (wall-clock GPU-h ÷ call count) was corrected the same session, `3d26d4a`, 19:25 JST, to use
  each call's own recorded `seconds` field instead (ISP-STD-1 E.3's "count GPU arms, not wall clock");
  a division error in the corrected estimate (149 vs. the correct 124.86 genes fitting a 10 GPU-h
  budget) was caught and fixed the same evening, `f945374`, 19:30 JST.
- **The estimable-count finding.** Among Amendment 4's 200 *drawn* genes, only 12 (6%) were estimable
  — catastrophically underpowered relative to the design's own power table. The human separately ruled
  (crossing with this finding) "allow 16 GPU-h, focus on 100 genes," read by both the orchestrator and
  the gatekeeper as **100 estimable genes**, not 100 drawn. Walking the same seed-20260929 draw, the
  100th estimable gene sits at position 984 (884 ineligible genes skipped at zero GPU cost).
- **Amendment 5** — the N=100-estimable redesign, registered `3c4a5a4`, 2026-09-29 19:49 JST; an
  out-of-sample/extrapolation check on the fitted cost model added as s.5.11, `e4c3a4e`, 19:52 JST;
  four run-package wiring defects (the driver still pointed at the old 200-gene file despite the
  registration text claiming otherwise; the no-op check was computed but never read; a bare `python3`
  instead of the pinned venv interpreter; missing hash-pins) fixed as s.5.12, `6fe5a84`, 20:01 JST; one
  further one-line fix to the no-op due-set computation (it must exclude `stopped_not_analysed` genes
  before taking every 20th, or a ceiling stop landing on a due position would wrongly force
  `no_op_failed`) as s.5.13, `36476f7`, 20:06 JST — this is the version the gatekeeper passed
  (registration and run package), and the version launched.
- **Amendment 6** — a confirmatory extension to N=200 estimable genes (positions 985–1791 of the same
  draw; the first 100 verified byte-identical to Amendment 5's own draw order, confirming the extension
  is a strict continuation, never a re-draw), registered `1b2f2fe`, 2026-09-30 02:27 JST. Framed
  explicitly as *extension-only*: the extension's GPU cost covers only the 100 new genes, reusing
  Amendment 5's already-computed markers for the first 100 — re-running all 200 from scratch would have
  cost materially more and required a separate escalation. Section 6.7, the decision rule for reading
  the N=100/N=200 pair (added verbatim per the gatekeeper's wording, since "both are reported, neither
  supersedes" left the headline undefined — an ISP-STD-1 A.1 multiplicity gap), committed `037683c`,
  2026-09-30 02:32 JST:
  1. The study's registered status is the Amendment 5 N=100 primary, α = 0.05.
  2. The N=200 set is a **nested extension**, never an independent replication.
  3. Primary positive **and** N=200 positive → "positive, confirmed at N=200."
  4. Primary positive **and** N=200 not positive → "positive (N=100), not confirmed at N=200."
  5. Primary not positive → the headline stays "negative (primary)"; a positive N=200 is reported
     separately as "N=200 extension positive, not a registered primary result" and never upgrades it.
  6. A `stopped_not_analysed` or `no_op_failed` outcome on either run is reported as that status for
     that run, and never overrides or is overridden by the other run's status.
  - A deviation note (not a numbered amendment, since no registered design/threshold/status changed)
    was appended after the gatekeeper's run-package re-gate found `apply_noop_gate` could not parse the
    real no-op result file's format (concatenated, indented JSON objects, not one object per line) —
    fixed, and disclosed as part of commit `c287178`, 2026-09-30 04:54 JST, together with two further
    fixes from the same gate: the combined-analysis script no longer trusts a supplied N=100 primary
    result file as given, instead recomputing it fresh from the same rows before using it (an
    authentication step, not a statistical change), and the two 100-gene ovx indexes are now checked
    for disjoint position keys before being merged.

**Cost model.** Fit by least squares directly from Phase 6's own recorded `seconds` + `n_token_cells`
fields (never from wall clock): delete ≈ 0.0688 + 0.0537 × cells; overexpress ≈ 4.309 + 0.0096 × cells
(overexpress carries a large near-fixed per-call cost, largely independent of detection). Cross-checked
out-of-sample against thinkstation2's independent Phase 6 data (−8.5% / −10.0% error, both well under
the pre-registered 50% refit trigger) and checked for extrapolation (the fit's training domain and the
100-selected-genes' domain are both bounded at [0, 100] by the cohort's own 100-cell-per-donor cap, so
no extrapolation occurs) — Amendment 5 s.5.11.

**Budgets and ceilings.** Amendment 5 (N=100 estimable): predicted 6.75 GPU-h raw / 10.12 GPU-h at a
flat 1.5× margin, against a 16 GPU-h ceiling (human-set) — 5.9 h of registered headroom. Amendment 6
(extension, 100 new genes only): predicted 6.80 GPU-h raw / 10.20 GPU-h margin-adjusted, registered as
the extension's own hard ceiling (36,720 s) — chosen specifically because it sits under the
orchestrator's 16 GPU-h escalation line, whereas a from-scratch 200-gene run's 20.3 GPU-h margin
estimate would not have.

**Power by simulation** (ISP-STD-1 C.3): an empirical-null critical value (40,000 simulations of
independent continuous data per N, one-sided α = 0.05 quantile) plus power estimation (20,000
simulations per effect size, Gaussian-copula construction mapping a target Spearman ρ to an underlying
Pearson r via r = 2 sin(ρ·π/6)), cross-checked analytically via the Fisher-z critical value
−1.645/√(N−1). At N=100: 43%/64%/81%/92% power for ρ = −0.15/−0.20/−0.25/−0.30. At N=200: 69%/88%/97%/100%
for the same effect sizes — the registered reason the confirmatory extension was judged worth its
additional GPU cost.

**Statistical test (both A5 and the A6 extension).** One-sided permutation test (100,000 permutations)
on the Spearman correlation between each gene's raw donor-level median delete shift and raw median
overexpress shift, H1: ρ < 0, over every estimable drawn gene; leave-one-out (2,000 permutations per
held-out gene) and a 10,000-draw gene-level bootstrap (2,000 permutations per draw) for stability;
positions due a no-op spot check (every 20th gene in run order, excluding any `stopped_not_analysed`
suffix) force `no_op_failed` on any failure or missing check, overriding the statistics; a gate-identity
check compares every post-GPU marker's token-positive cell count against its pre-GPU frozen value and
halts before any status is reported on any mismatch; an R3 validity check reruns the identical raw-shift
pipeline on the study's own 318 already-computed control genes as the closest available positive
control (a negative recovered ρ there is required before the null-study result is treated as trusted).

### 2.6 Reproducibility and provenance discipline

Every registration amendment is **append-only**: verified both by `git diff --numstat` (0 deletions on
the registration file) and by a byte-prefix sha256 check (the first N bytes of the new file must hash
identically to the entire pre-amendment file) before any push. Every script, weight file, token
dictionary, and data artifact a GPU run consumes is sha256-pinned in that run's launch script and
re-verified immediately before any GPU time is spent; a mismatch refuses the launch. A single-host
claim lock per run tree prevents two overlapping runs. Ceiling enforcement reads each call's own
recorded `seconds` field against the real process group (verified by `ps -o pid,pgid,cmd`, never
`pgrep -f`, which has repeatedly matched the launcher's own shell in earlier stages of this and sibling
studies) — never wall clock. Before the Amendment 6 extension's markers were treated as poolable with
Amendment 5's, a bridge check re-ran three already-completed (gene, donor) calls fresh into a scratch
directory and required bitwise-identical raw output against the original, both operations, before any
new gene was run (6/6 pairs matched). Before any analysis script reads Amendment 5's output tree, it
re-verifies a full sha256 manifest of that tree recorded at the gatekeeper's post-run clearance — and
that manifest file's own hash is itself pinned in the analysis code, so the manifest cannot be silently
substituted for one that always passes.

---

## Provenance table (partial — design and process facts only; extended with Results-section citations
once available)

| # | Claim | File | Commit / artifact |
|---|---|---|---|
| P1 | Raw donor cell-count imbalance 28–31× before capping; 1.0 after | `provenance/PHASE0_SELECTION_RULE.md` | `3af3ca98` |
| P2 | LuCA extended atlas download verified, sha256 `f0f7f434...` | download manifest, data root | Phase 2 close-out, 2026-09-24 09:06Z |
| P3 | 43-donor cohort, 8,600 analysis / 25,700 pool cells, seed 20260924 | `provenance/PHASE0_SELECTION_RULE.md`; cohort/tokenisation artifacts | `3af3ca98` |
| P4 | Model 316M bf16, no 104M arm (directive, not oversight) | `registration/PHASE5_ISP_REGISTRATION.md` Amendment 1/2 | `e2e10421`, `f5805cbc` |
| P5 | Calibrated throughput 8.19 cells/s; ≈20 s/gene ISP | `registration/PHASE5_ISP_REGISTRATION.md` Amendment 1 | `e2e10421` |
| P6 | 5-fold donor-stratified cross-fit design | `scripts/make_split.py`; registration Amendment 1 | `e2e10421` |
| P7 | Host split by gene; cross-host equivalence tolerance (max|Δ|≤1e-3, ρ≥0.999) | registration Amendment 3 | `a5fe96a7` |
| P8 | Panel A (15 genes), Panel B (36 genes) | registration s.1b | `552f0dfa` |
| P9 | Matched-control strata (6×20, later 19 strata / 318 unique) | registration Amendment 3c | `6ebf7fb1` |
| P10 | Phase 7 launch refusal (overexpress cell-order ambiguity) and its two-part fix | registration Amendment 3e/3f/3g/3h | `b9aafc6`, `10d0099`, `0f572cb`, `b2473ef` |
| P11 | Phase 6 duplicate-grid finding (both hosts ran a near-full independent copy) | Phase 6 output marker `"host"` field, both hosts' `phase6/` trees | read-only verification, 2026-09-29 |
| P12 | Null gene frame (19,902 genes) and draw mechanism (seed 20260929) | `scripts/select_null_genes.py`, `phase8_null/null_genes_*.json` | `622efaf` |
| P13 | Cost-model coefficients (delete/overexpress) and out-of-sample/extrapolation checks | registration Amendment 5 s.5.11 | `e4c3a4e` |
| P14 | Amendment 5 N=100-estimable design; 100th estimable gene at draw position 984 | registration Amendment 5 | `3c4a5a4` |
| P15 | Amendment 6 N=200 extension; s.6.7 pair-decision rule | registration Amendment 6 / s.6.7 | `1b2f2fe`, `037683c` |
| P16 | `apply_noop_gate` real-format parsing fix; deviation note | registration deviation note; `scripts/null_analysis.py` | `c287178` |
| P17 | Power-by-simulation table (N=100 and N=200) | registration Amendment 6 s.6.4 | `1b2f2fe` |
| P18 | ISP-STD-1 v1.1 adoption | `hive/standards/isp-outcome-criteria.md` | adopted 2026-09-28 13:35 JST |

*(Results-section rows — the achieved classifier gate value, Panel A/B outcome rows, and the null
study's A5/A6/combined statuses — are added once the gatekeeper's post-run check clears each output
file, citing that file's own path and sha256 rather than a number typed into this table.)*

---

## 3. Results

*PLACEHOLDER. Pending the independent gatekeeper's post-run check of:*
- *the closed Phase 4–7 analysis (classifier gate value, Panel A/B outcome table, GPU-hour totals with
  the duplicate-grid caveat applied) — already published via `RESULTS.md` and merged PR #38, restated
  here once cross-checked against this report's own provenance table;*
- *the null study's Amendment 5 (A5) N=100 primary result;*
- *the null study's Amendment 6 (A6) N=200 confirmatory-extension result and the s.6.7 pair headline.*

*No number from any of the above appears anywhere in this document until that check clears.*

---

## 4. Discussion

*PLACEHOLDER, pending §3.*

---

## 5. Limitations (design-level; independent of the pending result)

- **Ambient contamination.** Tumor-vs-normal comparison within a donor still differs in ambient RNA
  contamination from the tissue type itself; pairing removes study/chemistry/donor confounds but
  cannot remove this one. `REPLICATED_AMBIENT` and per-gene ambient flags (via leave-one-anchor-out) are
  the registered mechanism for disclosing this rather than absorbing it into a clean-looking headline.
- **Model change without a control arm.** The 316M model differs from the July screen's 104M model
  (measured ρ = 0.489 between them on a prior comparison); by directive, no 104M arm was run in this
  study, so any Panel A non-replication cannot be attributed to design vs. model — this wording was
  fixed in the registration (Amendment 2) before any output existed.
- **Donor overlap.** 20 of this study's 43 donors are among the core atlas's LUAD donor pool the July
  screen drew from; the design changed, not the underlying population of people.
- **Selection toward T-cell-rich samples.** The ≥100-T-cells-in-both-tissues donor rule structurally
  favors donors with more T-cell infiltration.
- **LUSC excluded.** Only 9 LUSC donors qualify on 10x within LuCA (below the registered ≥12 floor); a
  wider public-data survey was scoped but superseded by the human's choice to keep the registered
  2-class design.
- **Null study's extension is nested, not independent.** Amendment 6's N=200 set contains Amendment 5's
  N=100 primary genes; per the registered s.6.7 rule, it is read as a confirmatory extension of the
  same draw, never as an independent replication, and a positive extension result can never upgrade a
  non-positive primary headline.
- **Phase 6 GPU accounting reflects two near-complete overlapping runs**, not a clean complementary
  split, discovered from the output markers' own `"host"` field rather than from the runner's
  docstring — disclosed since it changes what the original ts1+ts2 GPU-hour close-out actually
  describes, though it does not change any already cross-host-verified result.

---

## 6. References

- ISP-STD-1 v1.1 — `hive/standards/isp-outcome-criteria.md` (adopted 2026-09-28).
- `balanced_donor_luad/registration/PHASE5_ISP_REGISTRATION.md` (all amendments, this repository,
  branch `analysis/balanced-donor-null-20260929`, current sha256 `465149ff...`).
- `balanced_donor_luad/provenance/PHASE0_SELECTION_RULE.md`.
- PR #38 (`Kays3/geneformer-lung-tcell`) — Phases 4–7, merged to `main` as `87349d41`.
- Salcher et al., *Cancer Cell* 2022 — the Lung Cancer Atlas, extended atlas.
- Geneformer V2 (vendored checkout, commit pinned in the run-launch hash-check dict).
