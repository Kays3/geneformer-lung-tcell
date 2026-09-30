# Balanced-Donor LUAD Tumor-vs-Normal T-Cell In Silico Perturbation Study

**Status: draft. Results and Discussion are placeholders pending the independent gatekeeper's
post-run check of the null-study analysis outputs. Nothing below states or implies any result of that
pending check, and no number in this document is a scientific outcome (a correlation, a p-value, or an
ISP-STD-1 status) unless it was already closed, published, and merged to `main` before this null study
began.**

Repository: `Kays3/geneformer-lung-tcell`. This document lives on branch
`report/imrad-balanced-donor-null` and is built from the registration documents, commit history, and
run artifacts cited inline, so a reader can check any design number against its source independently
of this prose.

---

## Abstract

*Placeholder, written after Results and Discussion are filled in.*

---

## 1. Introduction

Geneformer is a transformer model pretrained on single-cell transcriptomes. In silico perturbation
(ISP) uses it to estimate how forcing a gene's rank down (deletion) or up (overexpression) in a cell's
input token sequence shifts that cell's embedding toward or away from a target biological state. An
earlier screen over the Lung Cancer Atlas (LuCA) core atlas, referred to here as the July analysis, ran
an all-gene deletion screen against a three-class classifier (LUAD tumor, LUSC tumor, normal) built
from CELLxGENE's own disease labels. Read-only inspection on 2026-09-24 found that class definition to
be weaker than it first appears. The "normal" class drew from six non-cancer studies that contributed
no tumor T cells at all, so it was confounded with study identity rather than being a clean within-
population contrast. The "LUAD" class was selected by disease label alone, and in the core atlas that
label spans tumor, adjacent-normal, metastatic, and non-cancer tissue from LUAD patients, so the class
likely mixed several tissue sites without the mixing being visible in the label itself. Twenty of the
43 donors used in the present study overlap the donor pool the July screen drew its LUAD class from,
so the people are not new even though the design is.

The present study asks a narrower version of the July question. For each donor, does the model read a
consistent shift between that donor's own tumor tissue and their own adjacent-normal tissue, with every
other axis, study, chemistry, donor identity, and cell count, held equal between the two arms of the
comparison? Pairing within donor cancels study and chemistry by construction. A fixed per-donor,
per-tissue cell cap removes a raw 28 to 31-fold imbalance in cell counts across donors that existed
before capping (provenance table, row P1). Three questions follow from this design. The first, Panel A,
asks whether the 15 top LUAD-to-normal deletion hits from the July screen replicate under the
donor-paired design, in both the deletion and overexpression directions. The second, Panel B, asks
whether the model reads any T-cell-intrinsic signal on a curated panel of 36 canonical T-cell genes
that the July screen never ranked highly. The third is a null study: ranging over eligible random genes
detectable in LUAD tumor T cells, a population bounded by detectability in this cohort rather than
"genome-wide" in the usual sense (ISP-STD-1 B.9), does the same delete-overexpress raw-shift pipeline
show a systematic negative concordance even for genes with no a priori relationship to the tumor-normal
axis? A negative answer would suggest that any concordance seen in Panel A or B reflects gene-specific
biology rather than an artifact of the pipeline itself; a positive answer would narrow what such a
finding licenses.

All work after 2026-09-28 (13:35 JST) follows ISP-STD-1 v1.1, the floor standard the human adopted for
every in-silico-perturbation study going forward (`hive/standards/isp-outcome-criteria.md`). Its fixed
outcome vocabulary, `positive`, `negative`, `opposite_direction`, `not_estimable`,
`control_draw_sensitive_open`, `no_op_failed`, `not_run_by_design`, and `stopped_not_analysed`, is used
throughout this report. "No effect" is never substituted for `negative`, and "genome-wide" is never
substituted for "eligible random genes detectable in LUAD tumor T cells." Work completed before
2026-09-28 (Phases 4 through 7, described below) was not retroactively re-gated against ISP-STD-1,
consistent with the standard's own adoption note, though it already satisfied most of the standard's
requirements by construction: matched controls, a no-op gate, donor-level rather than cell-level
statistics, hash-pinned code and data, and ambient-signal flagging by leave-one-anchor-out.

---

## 2. Methods

### 2.1 Data source and cohort (Phases 0 through 3)

The data source is the single-cell Lung Cancer Atlas, extended atlas (Salcher et al., *Cancer Cell*
2022), CELLxGENE dataset version `33165751-ae33-4a65-94ee-52d9bc38f97e`, 1,283,972 cells across 17,764
genes with raw integer counts. It was downloaded on 2026-09-24 and verified against the published
object by S3 multipart ETag match and a full-file sha256 (`f0f7f43413bfe088281ac68487ad4bd73647d90572
80da202a82fabce9937cb9`; provenance row P2).

Donor and cell selection followed a mechanical rule fixed before any screening
(`provenance/PHASE0_SELECTION_RULE.md`, committed `3af3ca98` on 2026-09-24 at 18:19 JST, amended
`8ea11784` the same day to require a donor-level primary analysis, at least six paired donors held out,
and at least twelve eligible genes). The rule restricted disease status to lung adenocarcinoma, cell
type to CD4 and CD8 T cells (singlets, by the atlas's own `cell_type_major` label), chemistry to 10x
(3′ v2, with one contributing study also on 5′ v1), and tissue site to primary tumor versus
adjacent-normal only, excluding metastases, effusions, blood, and lymph nodes. A donor qualified only
if at least 100 T cells were available in both tissues.

This left 43 donors with both tissues represented, drawn from five studies: 23 from Leader_Merad 2021,
10 from Kim_Lee 2020, 5 from He_Fan 2021, 3 from Lambrechts_Thienpont 2018, and 2 from
Laughney_Massague 2020. Before any capping, 132,619 T cells met the eligibility rule, with donor-level
cell counts ranging from 237 to 6,764 in tumor tissue and 345 to 10,384 in adjacent-normal tissue, a
28 to 31-fold imbalance across donors. A fixed-seed draw (seed 20260924) then fixed exactly 100
analysis cells per donor per tissue (8,600 cells total, so every donor contributes an equal share) and
up to 300 pool cells per donor per tissue for training (25,700 cells total). Cells were tokenized for
Geneformer V2 with a median of 755 tokens per cell; the observed maximum of 2,656 tokens falls well
inside the model's 4,096-token context, so no truncation occurred. This cohort was committed on
`analysis/balanced-donor-luad` at `3af3ca98` (provenance row P3).

Three groups of cells were excluded, and the reasons are worth stating plainly rather than glossing
over. LUSC tumor tissue was excluded because only 9 donors qualify on 10x chemistry within LuCA,
short of the registered floor of 12; a broader public-data survey for a third class was scoped as a
possible follow-on (an early draft of Amendment 3) but was superseded once the human settled on the
registered two-class design (`e1a3ed94`, 2026-09-24 at 19:51 JST). Seventeen LUAD donors with only one
qualifying tissue were dropped, since the design requires both tissues from the same donor. And the
cohort as a whole is biased toward donors with more T-cell infiltration, because the donor rule
requires at least 100 T cells in both tissues; this selection effect cannot be removed by the pairing
design and is carried forward as a limitation (Section 5).

### 2.2 Model and fine-tuning design (Phase 4)

The model is Geneformer-V2-316M in bf16 precision, per the human's directive; no 104M arm was run in
parallel, a deliberate design choice registered in advance rather than an oversight (Amendment 1 and 2,
`e2e10421` and `f5805cbc`, both 2026-09-24). Because that comparison arm does not exist, any Panel A
non-replication is reported with a caveat fixed in the registration before any output existed: that a
non-replication on 316M does not by itself establish whether the cause is the redesigned cohort or the
change in model. Throughput was calibrated before any GPU time was spent on the real cohort, not as a
result of this study: roughly 8.19 cells per second for fine-tuning and about 20 seconds per gene for
ISP on 300 cells, both upper bounds measured with deterministic kernels disabled (Amendment 1,
`e2e10421`, 2026-09-24 at 18:47 JST).

Fine-tuning used a five-fold, donor-stratified, study-proportional cross-fitting design
(`scripts/make_split.py`, seed 20260924), so every donor is held out in exactly one fold and scored by
a model that never saw them, and both tissues from a donor always fall in the same partition. The
recipe was one epoch, learning rate 5e-5, batch size 8, the first six layers frozen, and seed 43. The
pre-registered gate required held-out balanced accuracy of at least 0.60, evaluated pooled across
folds, together with an exact sign test across all 43 donors (per-donor balanced accuracy above 0.5)
at p ≤ 0.05. The value the fine-tuned models actually achieved against this gate is a Results-section
quantity and is withheld here as a placeholder; see Section 3.

### 2.3 ISP design and panels (Phase 5 registration, Amendments 1 through 3h)

Two gene panels were fixed before any GPU output existed. Panel A consists of the 15 top LUAD-to-normal
deletion hits from the July screen, registered as a replication test and requiring at least 10 ambient
anchors per gene. Panel B consists of 36 canonical T-cell genes; three receptor-chain genes
(TRAC, TRBC1, TRBC2) were excluded because they are absent from the Geneformer V2 token dictionary, for
a reason that was not otherwise established and is not guessed at here. Both panels were registered at
`552f0dfa` on 2026-09-24 at 18:33 JST.

ISP calls were split by gene rather than by donor or by operation, registered as Amendment 3
(`a5fe96a7`, 2026-09-24 at 19:24 JST): each Ensembl ID was assigned alternately to thinkstation1 or
thinkstation2, and both operations for that gene ran on the assigned host. A cross-host equivalence
gate required that per-cell shifts computed on the same gene, donor, and operation on both hosts agree
within bf16 tolerance, a maximum absolute difference of 1e-3 and a Spearman correlation of at least
0.999, before any host-split data was pooled. Six shared control strata of 20 genes each, for secondary
use, were assembled with `matched_controls.py` (merged in from a sibling workstream, seed 20260924);
ambient status and contamination status were kept distinct throughout the study, so a gene such as
EEF1G could be flagged as ambient without being treated as evidence of contamination (Amendment 3c,
`6ebf7fb1`, 2026-09-24 at 20:07 JST).

The statistical design compares donor-level, control-adjusted shifts (never cell-level shifts) against
roughly 120 matched control genes per stratum, using an exact Wilcoxon signed-rank test over donors
with a Holm correction applied per family. A gene counts as replicated or as carrying a T-cell signal
only if both its deletion and overexpression arms are significant with opposite signs; same-signed
concordance between the two arms was registered in advance as insufficient evidence of specificity,
because the study's own pool of matched control genes already shows the same directional correlation.
A stability requirement, leave-one-control-out together with a 10,000-draw bootstrap, must each
preserve significance in at least 95% of draws for a call to stand; short of that, the result is
reported as `control_draw_sensitive_open` rather than as a clean positive.

One part of this history is worth recording in full rather than summarizing away, since it reflects
directly on what the eventual statistics can be trusted to mean. Phase 7 launch was refused on
2026-09-25, at roughly 14:45 JST, after read-only inspection found that Geneformer's own
`InSilicoPerturber` silently skips its stated token filter and instead sorts cells by length
internally, so the order of cells in the overexpression pickle did not match the order the registration
had assumed. This left it genuinely ambiguous which cells the overexpression arm's donor-level
statistic was summarizing. A second complication followed close behind: 12,184 of 15,142 gene-donor
pairs had tied token-positive lengths between the positive and negative donor pools, which made the
tie-free identity check first proposed for resolving the ambiguity infeasible on its own. Both issues
were resolved before any GPU output was read. A deterministic rule for reconstructing cell order was
written and verified bit-identical against a tie-free positive control spanning token-positive counts
from 10 to 100 (2,958 tie-free pairs existed to draw such a control from), and a rerun-by-construction
fallback was registered for any pair the reconstruction could not resolve. The sequence, a registered
stop, a root cause read before any fix was attempted, a verified reconstruction rule, and a disclosed
fallback, is recorded here as part of the study's gate history rather than smoothed into a footnote,
because ISP-STD-1's reproducibility requirements treat a caught pipeline-integrity gap as a finding
about the pipeline, distinct from any finding about biology.

### 2.4 Perturbation execution (Phase 6)

The perturbation run covered 171 genes (15 in Panel A, 36 in Panel B, and roughly 120 controls at the
time of launch, later extended to 318 unique control genes across 19 strata as the registration
matured), each perturbed by both deletion and overexpression, against all 43 donors' 100 tumor analysis
cells. Each donor was perturbed using their own fold's fine-tuned model and their own held-out
adjacent-normal pool as the goal centroid.

One finding from this stage is disclosed here rather than corrected quietly. The runner's own
docstring described the host split as gene-disjoint ("Amendment 3.1: split by gene"), but the `host`
field recorded inside every completed marker tells a different story: each host ran a near-complete
independent copy of the full 181-gene grid rather than a clean complementary half. Thinkstation1
produced 15,566 markers, exactly 181 genes times 43 donors times two operations, all stamped
`host=thinkstation1`; thinkstation2 produced 14,792, all stamped `host=thinkstation2`. This does not
change any result already reported, since the cross-host equivalence gate had independently verified
that both hosts agree on every shared gene-donor-operation triple, but it does change what the original
GPU-hour close-out for this phase actually describes: two overlapping near-full runs rather than two
complementary halves. The lesson carried forward into every later stage of this study is to check a
docstring's claimed data split against the provenance recorded in the output itself before treating
that claim as a cost-accounting boundary.

### 2.5 Null study design (Amendments 4 through 6)

The motivation for the null study came directly from the study's own matched-control data. The
delete-overexpress raw-shift pipeline already showed a negative correlation across the 120 to 318
matched control genes used elsewhere in the study (a stratum-adjusted correlation of roughly −0.62 is
reported in `RESULTS.md`). Because the stratum-adjustment machinery does not extend to genes outside
any registered control stratum, the null study asks the same concordance question on genes drawn
independently of any panel or control assignment, using the raw, non-stratum-adjusted shift.
Stratum-adjusting a newly drawn gene would require inventing a post-hoc matching rule, which the
registration explicitly avoids.

The gene frame for the draw is 19,902 genes: the Geneformer V2 dictionary of 20,271 genes, minus
Panel A's 15, minus Panel B's 36, minus the 318 unique genes ever used as a matched control across the
study's 19 strata. A single fixed-seed draw (seed 20260929) fixed the entire run order in advance, so
any budget stop always leaves a clean run prefix and a `stopped_not_analysed` suffix rather than a gap
in the middle. This draw was later extended, staying byte-identical on its shared prefix at every
extension, from 200 to 1,000 to 3,000 entries as the design's target gene count grew. Eligibility (at
least 10 donors with at least 10 token-positive cells) was checked read-only, at zero GPU cost, by
walking the draw order to find the Nth eligible gene; ineligible genes were skipped, never resampled.

The design history that follows is given with every commit hash and JST timestamp taken directly from
`git log`, not transcribed from memory, since accuracy here matters more than brevity. Amendment 4
registered an initial 200-gene draw (`622efaf`, 2026-09-29 at 19:20 JST). Its first cost estimate,
built from wall-clock GPU-hours divided by call count, was corrected the same evening to use each
call's own recorded `seconds` field instead, in line with ISP-STD-1's instruction to count GPU arms
rather than wall clock (`3d26d4a`, 19:25 JST); a simple division error in that corrected estimate, 149
genes fitting a 10 GPU-hour budget rather than the correct 124.86, was caught and fixed the same
evening (`f945374`, 19:30 JST).

A more consequential finding followed. Of Amendment 4's 200 drawn genes, only 12, or 6%, turned out to
be estimable, far short of what the design's own power table assumed. At close to the same time, and
independently, the human ruled that the null study could use up to 16 GPU-hours and should focus on
100 genes; both the orchestrator and the gatekeeper read this as 100 estimable genes rather than 100
drawn genes, since a literal 100-gene draw would yield only about five estimable genes and immediately
trip the orchestrator's own stop-and-ask rule for under 60 estimable genes. Walking the same
seed-20260929 draw, the hundredth estimable gene sits at draw position 984, with 884 ineligible genes
skipped along the way at no GPU cost.

Amendment 5 registered this N=100-estimable redesign at `3c4a5a4` (2026-09-29, 19:49 JST). An
out-of-sample and extrapolation check on the fitted cost model was added as section 5.11
(`e4c3a4e`, 19:52 JST). Four wiring defects in the run package, most seriously that the launch script
still pointed at the old 200-gene file even though the registration text claimed otherwise, plus a
no-op check that was computed but never actually read, a bare `python3` invocation instead of the
pinned virtual environment's interpreter, and missing hash pins, were fixed as section 5.12
(`6fe5a84`, 20:01 JST). One further one-line fix followed: the no-op due-gene set had to exclude any
gene already marked `stopped_not_analysed` before taking every twentieth gene, since otherwise a
ceiling stop landing exactly on a due position would wrongly force `no_op_failed` on an outcome that
had nothing wrong with it (section 5.13, `36476f7`, 20:06 JST). This was the version the gatekeeper
passed, for both the registration and the run package, and the version that was launched.

Amendment 6 registered a confirmatory extension to 200 estimable genes, covering draw positions 985
through 1791, with the first 100 verified byte-identical to Amendment 5's own draw order to confirm the
extension is a strict continuation of the same draw rather than a fresh one (`1b2f2fe`, 2026-09-30 at
02:27 JST). The extension was framed deliberately as extension-only: its GPU cost covers only the 100
newly drawn genes, reusing Amendment 5's already-computed markers for the first 100 rather than
re-running them, since a from-scratch 200-gene run would have cost materially more and required a
separate escalation past the orchestrator's GPU-hour threshold. Section 6.7 of the registration states
the rule for reading the resulting pair of results, added in the gatekeeper's own wording after he
pointed out that the earlier text, "both are reported, neither supersedes," left the eventual headline
undefined, a multiplicity gap under ISP-STD-1's pre-registration requirement (committed `037683c`,
2026-09-30 at 02:32 JST). The rule has six parts. The study's registered status is the Amendment 5
N=100 primary result at alpha 0.05. The N=200 set is a nested extension of that primary result, never
an independent replication. If the primary result is positive and the N=200 result is also positive,
the finding is reported as positive, confirmed at N=200. If the primary result is positive but the
N=200 result is not, the finding is reported as positive at N=100, not confirmed at N=200. If the
primary result is not positive, the headline remains negative at the primary level regardless of what
the extension shows; a positive extension result is reported separately as an N=200 extension finding
that is not a registered primary result, and it never upgrades the headline. Finally, a
`stopped_not_analysed` or `no_op_failed` outcome on either run is reported as that status for that run
specifically, and never overrides or is overridden by the other run's status.

A further change followed from the gatekeeper's run-package review rather than from any new
registration decision, and is recorded as a deviation note rather than a numbered amendment, since it
altered no registered design, threshold, or status definition. The review found that `apply_noop_gate`
could not parse the real no-op result file, because the file format is concatenated, indented JSON
objects rather than one object per line; this was fixed, along with two related findings from the same
review, that the combined-analysis script had been trusting a supplied N=100 primary result file at
face value rather than recomputing it from the same rows before use, and that the two 100-gene indexes
built for the primary and extension gene sets were not checked for overlapping keys before being
merged. All three fixes are disclosed together in commit `c287178` (2026-09-30, 04:54 JST).

The cost model behind these budgets was fit by ordinary least squares directly from Phase 6's own
recorded `seconds` and `n_token_cells` fields, never from wall clock: roughly 0.0688 plus 0.0537 times
cell count for deletion, and roughly 4.309 plus 0.0096 times cell count for overexpression, the latter
carrying a large near-fixed cost per call that is largely independent of how many cells are detected.
The fit was checked out of sample against thinkstation2's independently recorded Phase 6 data, with
errors of −8.5% for deletion and −10.0% for overexpression, both comfortably under the 50% threshold
that would have triggered a refit, and checked for extrapolation risk, finding none, since both the
fit's training domain and the 100 selected genes' own domain are bounded at 100 cells by the cohort's
own per-donor cap (Amendment 5, section 5.11). Under this model, Amendment 5's 100 estimable genes were
predicted to cost 6.75 GPU-hours raw, or 10.12 GPU-hours with a flat 1.5-fold margin, against a
human-set ceiling of 16 GPU-hours, leaving 5.9 hours of registered headroom. Amendment 6's extension,
covering only the 100 new genes, was predicted to cost 6.80 GPU-hours raw, or 10.20 GPU-hours with the
same margin, and that margin-adjusted figure was registered as the extension's own hard ceiling
(36,720 seconds), chosen specifically because it falls under the orchestrator's 16-hour escalation
line, whereas a from-scratch 200-gene run's margin estimate of roughly 20.3 GPU-hours would not have.

Power was estimated by simulation, following ISP-STD-1 section C.3: an empirical null critical value
from 40,000 simulations of independent continuous data at each sample size, combined with power
estimates from 20,000 simulations per effect size using a Gaussian-copula construction that maps a
target Spearman correlation to an underlying Pearson correlation, and cross-checked analytically
against the Fisher-z critical value. At N=100, power ranged from 43% at a true correlation of −0.15 to
92% at −0.30. At N=200, the same range of effect sizes gave power from 69% to essentially 100%. This
gain in power was the registered reason the confirmatory extension was judged worth its additional GPU
cost.

Both the Amendment 5 primary run and the Amendment 6 extension use the same statistical test: a
one-sided permutation test, with 100,000 permutations, on the Spearman correlation between each gene's
raw donor-level median deletion shift and its raw median overexpression shift, with the alternative
hypothesis that the correlation is negative, computed over every estimable drawn gene. Stability was
assessed by leave-one-out, recomputing the full permutation test with 2,000 permutations for each
held-out gene, and by a 10,000-draw gene-level bootstrap, each draw again using 2,000 permutations. A
no-op spot check runs on every twentieth gene in run order, excluding any genes already
`stopped_not_analysed`; any failed or missing check forces the overall result to `no_op_failed`,
overriding whatever the statistics alone would say. A gate-identity check compares every post-GPU
marker's token-positive cell count against its pre-GPU frozen value and halts before any status is
reported if the two disagree. A validity check reruns the identical raw-shift pipeline on the study's
own 318 already-computed control genes as the closest available positive control; recovering a negative
correlation there is a precondition for treating the null-study result as trustworthy.

### 2.6 Reproducibility and provenance discipline

Every amendment to the registration document is append-only, verified both by `git diff --numstat`
showing zero deletions to the registration file and by a byte-prefix sha256 check, in which the first N
bytes of the amended file must hash identically to the entire pre-amendment file, before any push. Every
script, weight file, token dictionary, and data artifact that a GPU run consumes is sha256-pinned in
that run's launch script and re-verified immediately before any GPU time is spent, and a mismatch
refuses the launch outright. A single-host claim lock on each run's output tree prevents two runs from
overlapping. Ceiling enforcement reads each call's own recorded `seconds` field against the real
process group, verified with `ps -o pid,pgid,cmd` rather than `pgrep -f`, which has repeatedly matched
a launcher's own shell process in earlier stages of this and related studies, and never from wall
clock. Before Amendment 6's markers were treated as poolable with Amendment 5's, a bridge check
re-ran three already-completed gene-donor pairs fresh into a scratch directory and required the raw
output to be bitwise identical to the original for both operations before any new gene was run; all
six checks passed. Before any analysis script reads Amendment 5's output tree, it re-verifies a full
sha256 manifest of that tree recorded at the gatekeeper's post-run clearance, and the manifest file's
own hash is itself pinned inside the analysis code, so the manifest cannot be silently swapped for one
that always passes.

---

## Provenance table

This table is partial. It covers design and process facts established before this report was written;
rows drawn from the pending analysis outputs will be added, citing each output file's own path and
hash, once the gatekeeper's check clears.

| # | Claim | File | Commit or artifact |
|---|---|---|---|
| P1 | Raw donor cell-count imbalance of 28 to 31-fold before capping, 1.0 after | `provenance/PHASE0_SELECTION_RULE.md` | `3af3ca98` |
| P2 | LuCA extended atlas download verified, sha256 `f0f7f434...` | download manifest, data root | Phase 2 close-out, 2026-09-24 09:06Z |
| P3 | 43-donor cohort, 8,600 analysis and 25,700 pool cells, seed 20260924 | `provenance/PHASE0_SELECTION_RULE.md`; cohort and tokenization artifacts | `3af3ca98` |
| P4 | Model is 316M bf16 with no 104M arm, by directive | `registration/PHASE5_ISP_REGISTRATION.md`, Amendments 1 and 2 | `e2e10421`, `f5805cbc` |
| P5 | Calibrated throughput: 8.19 cells/s fine-tuning, about 20 s/gene ISP | Amendment 1 | `e2e10421` |
| P6 | Five-fold donor-stratified cross-fit design | `scripts/make_split.py`; Amendment 1 | `e2e10421` |
| P7 | Host split by gene; cross-host tolerance max\|Δ\|≤1e-3, ρ≥0.999 | Amendment 3 | `a5fe96a7` |
| P8 | Panel A (15 genes) and Panel B (36 genes) | registration section 1b | `552f0dfa` |
| P9 | Matched-control strata, six of 20 genes, later 19 strata / 318 unique | Amendment 3c | `6ebf7fb1` |
| P10 | Phase 7 launch refusal (overexpression cell-order ambiguity) and its fix | Amendments 3e, 3f, 3g, 3h | `b9aafc6`, `10d0099`, `0f572cb`, `b2473ef` |
| P11 | Phase 6 duplicate-grid finding: both hosts ran a near-full independent copy | output marker `host` field, both hosts' `phase6/` trees | verified read-only, 2026-09-29 |
| P12 | Null gene frame of 19,902 genes and draw mechanism, seed 20260929 | `scripts/select_null_genes.py`, `phase8_null/null_genes_*.json` | `622efaf` |
| P13 | Cost-model coefficients and out-of-sample/extrapolation checks | Amendment 5, section 5.11 | `e4c3a4e` |
| P14 | Amendment 5 N=100-estimable design; hundredth estimable gene at draw position 984 | Amendment 5 | `3c4a5a4` |
| P15 | Amendment 6 N=200 extension and the section 6.7 pair-decision rule | Amendment 6, section 6.7 | `1b2f2fe`, `037683c` |
| P16 | `apply_noop_gate` real-format parsing fix and its deviation note | registration deviation note; `scripts/null_analysis.py` | `c287178` |
| P17 | Power-by-simulation table, N=100 and N=200 | Amendment 6, section 6.4 | `1b2f2fe` |
| P18 | ISP-STD-1 v1.1 adoption | `hive/standards/isp-outcome-criteria.md` | adopted 2026-09-28, 13:35 JST |

---

## 3. Results

Placeholder, pending the independent gatekeeper's post-run check of three outputs: the closed Phase 4
through 7 analysis (the achieved classifier-gate value, the Panel A and B outcome table, and GPU-hour
totals with the duplicate-grid caveat applied, already published in `RESULTS.md` and merged as PR #38,
restated here once cross-checked against this report's own provenance table); the Amendment 5 N=100
primary result; and the Amendment 6 N=200 confirmatory-extension result together with the section 6.7
pair headline it produces. No number from any of these three outputs appears anywhere in this document
until that check clears.

---

## 4. Discussion

Placeholder, pending Section 3.

---

## 5. Limitations

Several limitations are inherent to the design and do not depend on what the pending results turn out
to show. Ambient RNA contamination differs between tumor and adjacent-normal tissue within the same
donor in ways that pairing cannot remove, since pairing controls for study, chemistry, and donor
identity but not for tissue-of-origin contamination; the study's `REPLICATED_AMBIENT` label and its
per-gene ambient flag, computed by leave-one-anchor-out, exist specifically to disclose this rather
than let it disappear into a clean-looking headline. The 316M model differs from the 104M model used in
the July screen (a prior comparison measured a correlation of 0.489 between the two), and because no
104M arm was run here by directive, any non-replication on Panel A cannot be attributed cleanly to the
redesigned cohort as opposed to the change in model; this caveat was written into the registration
before any output existed, not added after the fact. Twenty of the 43 donors in this study overlap the
donor pool the July screen's LUAD class drew from, so the two studies are not independent at the level
of the underlying population. The cohort is structurally biased toward donors with more T-cell
infiltration, since donors need at least 100 T cells in both tissues to qualify. LUSC tissue was
excluded because only nine donors qualify on 10x chemistry within LuCA, short of the registered floor
of twelve; a wider public-data survey for a third class was scoped but not pursued once the human
settled on the two-class design. The null study's Amendment 6 extension is nested inside Amendment 5's
primary gene set rather than independent of it, and the registered section 6.7 rule treats it strictly
as a confirmatory extension of the same draw: a positive extension result can never upgrade a
non-positive primary headline. And the Phase 6 GPU-hour accounting reflects two near-complete
overlapping runs on the two hosts rather than a clean complementary split, a fact discovered from the
output markers' own host field rather than from the runner's docstring; this changes what the original
combined GPU-hour figure actually describes, though it does not change any result that the cross-host
equivalence check had already verified.

---

## 6. References

- ISP-STD-1 v1.1, `hive/standards/isp-outcome-criteria.md`, adopted 2026-09-28.
- `balanced_donor_luad/registration/PHASE5_ISP_REGISTRATION.md`, all amendments, this repository,
  branch `analysis/balanced-donor-null-20260929`, current sha256 `465149ff...`.
- `balanced_donor_luad/provenance/PHASE0_SELECTION_RULE.md`.
- PR #38, `Kays3/geneformer-lung-tcell`, Phases 4 through 7, merged to `main` as `87349d41`.
- Salcher et al., *Cancer Cell* 2022, the Lung Cancer Atlas, extended atlas.
- Geneformer V2, vendored checkout, commit pinned in the run-launch hash-check dictionary.
