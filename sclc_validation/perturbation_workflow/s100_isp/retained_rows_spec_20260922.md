# Retained-row analysis contract — S100 ISP

Build `s100_isp_retained_rows.csv` before any result filtering. Its grain is
one planned `gene × source × goal × perturbation` row. For Module A, the fixed
two sources (SCLC, LUAD), two goals per source, 12 core S100 genes, and two
operations yield **96 core rows**. Add control rows at the same grain after
the six fixed 20-control strata have been locked.

The table is the left join of the planned Cartesian design, per-gene completion
markers, paired-cell eligibility manifest, no-op manifest, and ISP stats. A
missing stats row must remain a row; it is never filtered away.

Required columns:

| Column | Meaning |
|---|---|
| `panel_id`, `run_id`, `runner_sha256`, `panel_sha256` | Exact provenance. |
| `gene`, `ensembl_id`, `role`, `stratum` | Pre-registered panel membership. |
| `source_state`, `goal_state`, `alt_state`, `perturbation_type` | Planned contrast and arm. |
| `status` | `eligible_completed`, `not_estimable_cell_count`, `not_estimable_donor_count`, `not_estimable_control_stratum`, `run_failed`, or `no_op_failed`. |
| `status_reason` | Machine-readable completion-marker or preflight reason. |
| `n_token_positive_cells`, `n_eligible_donors`, `donor_cell_counts` | Gate evidence before the 100-cell-per-donor cap. |
| `paired_cell_manifest_sha256`, `paired_cell_count` | Proof delete/OE used the same deterministic cells. |
| `no_op_status`, `no_op_score` | Required no-op diagnostic result. |
| `raw_completion_marker`, `raw_stats_path`, `raw_stats_sha256`, `stats_row_present` | Output provenance and explicit missingness. |
| `donor_balanced_shift`, `donor_sign_fraction`, `matched_control_percentile_q` | Nullable analysis outputs; null is retained when a row is not estimable. |

`status` is determined before the stats join. A per-gene completion marker
applies to both goal rows for that source/operation. The final report must give
counts for every status, list all 96 core rows, and state that unestimable rows
were retained rather than treated as zero effects.

## Amendment — 2026-09-23 (ruled by Michael, s100-isp-execution-20260922)

Three gaps surfaced while building the analysis layer against this spec.
Amended visibly here, per the same standard as the preflight hash
correction — never silently edit the enum or the rule above.

**1. `status` gets a seventh value, `not_run`, outside the six listed above.**
Every value in the `status` column's enum describes a *completed* outcome;
none describes "eligible, but the run has not happened yet" — a state that
necessarily exists before Module A executes and while it is partially
complete. `not_run` labels exactly that. **Ruling: this changes no
threshold, no statistic, and no interpretation — it labels a state that
exists prior to any result, which is exactly what pre-registration
doctrine is meant to allow, not what it exists to prevent.** Implemented in
`retained_rows.py` as `STATUS_NOT_RUN`.

**2. `donor_balanced_shift` must be computed from the raw per-cell pickle
plus the paired-cell manifest's donor column, never from
`InSilicoPerturberStats`' own `Shift_to_goal_end` column.** That library
column is a plain cell-weighted mean over all cells in the arm, not the
donor-weighted mean this spec's grain and the design doc both require ("the
donor, not the cell, is the unit of replication"). Concretely, on the arm
this decides: LUAD `S100A2` has four eligible donors at 40, 4, 4, and 25
cells. A cell-weighted mean gives them 54.8%, 5.5%, 5.5%, and 34.2% of the
result — the 40-cell donor at **2.19x** its intended (donor-balanced, 25%
each) weight, each 4-cell donor at 0.22x. `S100A2` is the only
ambient-flagged gene surviving the LUAD arm at all. This is not a rounding
difference, and it is not a hypothetical: the design doc names this
failure mode explicitly as one that has "already happened once" in the
existing lung analysis. `raw_stats_path`/`raw_stats_sha256`/
`stats_row_present` remain in the table as output provenance, but they are
provenance only — they are never the source of `donor_balanced_shift`.

**3. The exact percentile formula behind `matched_control_percentile_q`
is fixed as documented in `ambient_stats.py`'s `midrank_percentile()`
(value's rank among the stratum's N controls plus itself, ties averaged,
mapped to `(rank - 0.5) / (N + 1) * 100`), and a row's `status` is never
overwritten to `not_estimable_control_stratum` when its own stratum's
control set fails the 20-count gate — `status` continues to describe only
that gene's own run outcome, and `matched_control_percentile_q` simply
stays null. Ruling on both, and why they are one decision, not two:**
`Q` is used only in rank-based/relative contexts everywhere it appears in
the design (the Spearman correlation, its permutation test, the LOO and
bootstrap recomputations, and the "flagged genes have a higher median `Q`
than clean controls" comparison) — never against an absolute threshold.
Every gene is scored against the **same fixed N = 20** for its stratum, so
any monotone rank-to-percentile map applied identically to all ten
non-anchor genes preserves both the cross-gene ordering and the median
comparison; the registered rho is numerically unchanged by which monotone
formula is used. **That inertness holds only because N is fixed at 20 for
every gene, with no exceptions.** The companion ruling — refusing to
compute `Q` on a short control set rather than "rescuing" a gene by
computing it on however many controls survived — is precisely what
guarantees N stays fixed at 20 everywhere `Q` is used. If a future change
ever computed `Q` for one gene on fewer than 20 controls, the percentile
formula would stop being a common monotone map across genes and could
start moving the registered rho; do not make that change without revisiting
this amendment.

## Amendment 2 -- 2026-09-23 (human ruling, s100-isp-execution-20260922)

Five items ruled together, unblocking the analysis layer's remaining stub
and setting the shape of the real run. Amended visibly here, same
discipline as Amendment 1 above -- append only, the original text and
Amendment 1 are unedited.

**1. The SCLC source arm is dropped.** In SCLC, none of the ambient-flagged
genes survive eligibility, leaving a two-cluster ambient-risk variable
where a pure two-group difference passes the primary `rho >= 0.70` gate
47% of the time -- the arm could return a passing number that carries no
information about ambient risk, whether or not it happened to pass.
Dropping a source is not dropping a goal state: `SCLC` remains a valid
**goal**. Module A's real invocation is source = LUAD only, goals
LUAD -> SCLC and LUAD -> Normal (both registered). Core row count is now
12 genes x 1 source x 2 goals x 2 operations = **48**, not the 96 in this
document's opening section, which described the two-source design as
originally registered and is superseded for Module A's actual invocation
by this amendment (never edited above -- `retained_rows.py` remains fully
generic over `sources` and its own test suite still exercises both source
branches to prove that generic code path is correct, independent of which
sources a given call site chooses to invoke). `retained_rows.py`'s
`MODULE_A_SOURCES` constant is updated to `("luad",)` accordingly.

**2. Median token rank (matched_controls.py's open question 1) is ruled:**
it is not in `ambient_risk_all_genes.csv` and detection fraction is not
substituted for it -- they are independent matching axes. It is computed
from the tokenized held-out dataset, `ALLGENE_ROOT/data/heldout_test.dataset`
(the same source the eligibility counts come from): for each gene, the
median, over cells in which that gene is token-positive, of the gene's
position in that cell's rank-ordered `input_ids`. **Convention: 0-based**
(Python's own list-index convention --
`matched_controls.median_token_rank_and_detection()`). This does not
affect matching, which is relative, but is frozen here so it is never a
live argument. The frozen per-gene table (detect_frac + median_token_rank,
LUAD-only) must be hashed from its own bytes at the moment it is written
(`matched_controls.load_and_freeze_luad_gene_stats()`, via
`provenance_utils.sha256_file`), same principle as every other input this
layer reads.

**3. Matching scope (open question 2) is ruled: per-source-state, not
global.** With SCLC dropped this means, in practice: the matched-control
table is built from LUAD held-out cells only. Both matching quantities are
source-dependent -- a gene's detection fraction and median token rank both
move with the source state -- so a global table would match on an average
that describes neither state, quietly rather than loudly wrong.
`ambient_risk_all_genes.csv`'s `detect_frac` is a global (cross-source)
figure and is therefore NOT used for the detection-fraction matching axis
either, for the same reason, not because detection fraction changed
meaning. `matched_controls.build_matched_control_table()` is unblocked: it
no longer raises `NotImplementedError`; its module docstring's old
"STATUS: BLOCKED" section (correct while the questions were open) is
replaced with "STATUS: UNBLOCKED 2026-09-23" pointing back to this
amendment, since it described a stub that no longer exists.
`synthetic_control_table()` is untouched and remains in use by
`ambient_stats.py`'s own tests.

**4. The primary test's real `n` is eight, not ten -- ruled to still pass,
with a correction to the supporting math first circulated for it.** Of the
twelve panel genes, ten are non-anchor; in LUAD only eight of those ten
are eligible (`S100P` and `S100A16` both fail the >=50-cell gate in LUAD --
confirmed against the live preflight manifest,
`s100_preflight_eligibility.csv`, built earlier under card
`isp-runner-paired-arms-20260922` gate 3(a)). The exact null distribution
at n=8 was recomputed independently by direct enumeration of all
`8! = 40,320` label permutations (matching
`ambient_stats.spearman_exact_permutation`'s own method, not an asymptotic
approximation), because the numbers first circulated for this ruling used
a one-sided convention while the design doc registers the primary test as
**two-sided** ("exact two-sided permutation p-value... the gate is
rho >= 0.70 and p <= 0.05," design doc, "Controls and analysis"). Corrected
exact two-sided table:

| n | P(rho >= 0.70), one-sided | min rho for exact two-sided p <= 0.05 |
|---|---|---|
| 10 | 0.0134 | 0.6485 |
| 8 | 0.0288 | 0.7381 |

**This changes the conclusion at n=8.** At n=10 the two-sided threshold
(0.6485) sits comfortably below the 0.70 rho-gate, so rho >= 0.70 was
correctly the binding constraint there. At n=8 it does not: the minimum
rho for exact two-sided p <= 0.05 is ~0.7381, HIGHER than the 0.70
rho-gate. An observed rho in `[0.70, 0.7381)` at n=8 clears the rho-gate
and still fails the combined gate on the p-value alone -- the two
conditions are not redundant at this n, and "rho >= 0.70 remains the
binding constraint" (the claim first offered for this ruling) is not
correct for the real n=8 run. **Ruling: the registered gate
(`rho >= 0.70 AND p <= 0.05`, both independently checked) is unchanged and
is not being loosened or tightened for n=8** -- this correction is to the
supporting arithmetic offered as reassurance, not to the pre-registered
gate itself, and the gate as implemented already handles this correctly by
construction (it checks both conditions, never assumes one from the
other). What changes is that the p-value condition can no longer be
assumed automatic once rho clears 0.70 at n=8, and
`ambient_stats.primary_test()`'s own report must always give the exact
two-sided p and never assume it from rho alone. `primary_test()` now also
reports `n` explicitly alongside `rho`/`p_exact` -- never quote an n=10
p-value for an n=8 test.

**5. Leave-one-gene-out is added to the primary test, registered now,
before the real number exists.** Of the eight LUAD-eligible non-anchor
genes, exactly one -- `S100A2` -- is ambient-flagged (confirmed against
the same live preflight manifest as item 4). If the other seven cluster
together in ambient risk, the observed rho is decided by where `S100A2`
lands -- a single comparison wearing the clothes of a rank correlation.
**Rule, registered before the number exists:** recompute the primary
Spearman rho with `S100A2` removed (rho only, via the same cheap
rank-Pearson used by the leave-one-control-out/bootstrap loops -- not a
second exact-permutation p, which is not the registered check here). If
the full eight-gene primary rho passes the gate but the leave-`S100A2`-out
rho (n=7) falls below **0.60** (the same floor as leave-one-control-out,
not a new number chosen for this situation), the result is declared
**"carried by a single gene"** and may **not** be stated as an
ambient-risk association across the panel -- only as an `S100A2`-specific
finding. Reported always, not only on request. Implemented as
`ambient_stats.primary_test_with_single_gene_check()`.

**CPU precheck (mandatory before any GPU work, per this amendment's own
order of work): ambient-risk values of the eight LUAD-eligible non-anchor
genes**, pulled directly from the existing global
`ambient_risk_all_genes.csv` (this table's cross-source `detect_frac` is
not used for the LUAD-only control matching per item 3 above, but
`ambient_risk` itself is the pre-registered primary-test x-variable, not a
matching axis, and is unaffected by that ruling):

| gene | ambient_risk | ambient_pct |
|---|---|---|
| S100A4 | 0.0000011 | 0.94 |
| S100A6 | 0.000034 | 1.52 |
| S100A11 | 0.000144 | 1.95 |
| S100A10 | 0.000312 | 2.23 |
| S100PBP | 0.744492 | 44.73 |
| S100A13 | 0.795706 | 49.90 |
| S100B | 0.859066 | 57.48 |
| S100A2 | 0.955296 | 79.68 |

**Spread: strongly bimodal, not smoothly distributed.** Four values sit
within `3.1e-4` of zero (the clean-control cluster); the other four sit
between `0.744` and `0.955`; the gap between the two clusters
(`0.744 - 0.000312 = 0.744`) is larger than the full range within either
cluster. This is a materially different shape from "eight roughly evenly
spread values" and is worth naming explicitly: it is not the same failure
mode that dropped the SCLC arm (there, no ambient-flagged gene survived
eligibility at all; here, `S100A2` -- the one ambient-flagged survivor --
sits at the extreme end of a real, non-degenerate four-gene high cluster,
and the within-cluster ordering among `{S100PBP, S100A13, S100B, S100A2}`
and among the four clean controls is genuine information, not an
artifact). But it does mean the eight-point rank correlation is closer to
"two clusters plus within-cluster ordering" than to a smooth rank spread,
which is exactly why item 5's leave-`S100A2`-out check matters here and
was registered before this number was computed, not after.

## Amendment 3 -- 2026-09-23 (human self-correction, s100-isp-execution-20260922)

Corrects two numbers in Amendment 2 above, same day, before any GPU work
started. Amendment 2's text is left exactly as written -- both errors are
corrected here, visibly, not edited into the original.

**1. Amendment 2, item 4's "min rho for exact two-sided p <= 0.05" table
was itself computed one-sided and is WRONG. It is superseded by this
corrected, independently-verified table (enumerated exactly, no floating
point, cross-checked by two independent parties by direct permutation
enumeration):**

| n | attainable rho just below the boundary | its exact two-sided p | smallest attainable rho that clears p <= 0.05 | its exact two-sided p |
|---|---|---|---|---|
| 10 | 0.636364 | 0.05443 (fail) | **0.648485** | 0.04898 |
| 8 | 0.714286 | 0.05759 (fail) | **0.738095** | 0.04583 |

The values 0.6485 (n=10) and 0.7381 (n=8) in Amendment 2 were correct as
written -- this amendment does not change them, it confirms them by an
independent method (exact integer `d^2`/`Fraction` enumeration, no
floating-point comparisons) and corrects the *reasoning* Amendment 2 used
to justify them, which was still using a one-sided framing in prose even
though its own numbers were two-sided. **The registered gate
(`rho >= 0.70 AND p <= 0.05`) does not change.** At n=8, `rho >= 0.70`
alone does NOT guarantee `p <= 0.05` -- the effective floor for the
combined gate is 0.7381, not 0.70. Both conditions must be evaluated
independently and explicitly; do not assume the p-condition from the
rho-condition at n=8.

**2. Amendment 2, item 5's leave-one-gene-out floor (0.60, borrowed from
leave-one-control-out) is WRONG and is replaced.** Leave-one-control-out
keeps n fixed (it drops a control, not a gene) -- its 0.60 floor has no
bearing on leave-one-*gene*-out, which shrinks the primary test's n from 8
to 7. This was a substitution error of the same kind Amendment 2 itself
warned against for token rank (using the wrong quantity because it was
convenient, not because it was correct).

A first replacement of 0.60 with 0.75 was also wrong, by one attainable
step: **Spearman's rho at small n is discrete -- for n=7, `d^2` is always
an even integer, so attainable rho values near the boundary are exactly
0.678571, 0.714286, 0.750000, 0.785714, and nothing between.** 0.750000's
own exact two-sided p is 0.06627, which FAILS `p <= 0.05` -- it is not the
critical value, it is one attainable step short of it. **The correct
value is 0.785714 (11/14, `d^2 = 12`), whose exact two-sided p is 0.04801
and clears the gate.** General lesson, worth keeping for the stability
code broadly: **find a critical value by scanning attainable statistic
values for the first whose own exact p clears the threshold -- never by
indexing into a sorted null at an approximate quantile position**, which
can land inside a gap between attainable values and silently pick the
wrong side of the boundary.

**Corrected rule, replacing Amendment 2 item 5, evaluated two-sided like
the primary test (not the cheap rho-only recomputation the leave-one-
control-out/bootstrap loops use -- this runs once, at n=7, so the full
exact permutation test costs nothing to compute properly):**

- **leave-`S100A2`-out rho >= 0.785714 AND its own exact two-sided
  p <= 0.05: "survives"** -- the ambient-risk association is not solely
  attributable to `S100A2`; a panel-wide statement is permitted (still
  subject to the full eight-gene test's own gate).
- **0 < leave-`S100A2`-out rho, but its exact p > 0.05: "gene_sensitive_open"**
  -- no panel-wide claim, and no denial either. Matches the design's
  existing "control-draw-sensitive / open" vocabulary rather than
  inventing a new one.
- **leave-`S100A2`-out rho <= 0: "carried_by_single_gene"** -- the
  association does not survive `S100A2`'s removal at all; may NOT be
  stated as a panel-wide ambient-risk finding, only as an `S100A2`-specific
  one.

Nothing here is borrowed from an unrelated gate; both the 0.7381 (item 4)
and 0.785714 (item 5) values are derived from the design's own registered
two-sided convention at the n the real analysis actually has (8 and 7
respectively). Implemented as `ambient_stats.leave_one_gene_out_check()`
(three-way `outcome` field plus a `carried_by_single_gene` boolean for
callers that only need the one flag) and
`ambient_stats.primary_test_with_single_gene_check()` (merges
`primary_test()` and `leave_one_gene_out_check()`). Reported always, not
only on request, same as Amendment 2 registered.

## Amendment 4 -- 2026-09-23 (human ruling, GPU held, s100-isp-execution-20260922)

Registered before GPU work starts and before any real result exists.
Amendments 1-3 stand unedited.

**GPU IS HELD** pending the human's review of the finding below. Nothing
in this amendment authorizes or resumes GPU work; the CPU precheck in
Amendment 2 did its job and surfaced a structural problem worse than the
one that dropped the SCLC arm.

**1. The bimodal ambient-risk spread reported in Amendment 2's CPU
precheck is NOT itself the problem -- the group split is.** Spearman's
rho depends only on ranks; it cannot see that the gap between the two
clusters (0.744) is larger than either cluster's own range, so Amendment
2's framing of that gap as the concerning feature was imprecise. **What
matters is that the 8 genes split into a 4-gene ambient-high cluster and a
4-gene ambient-low cluster at all.** Verified independently by direct
enumeration of all 4! x 4! = 576 within-cluster orderings (a pure
between-cluster Q difference, contributing zero within-cluster rank
information):

```
mean rho over the 576 arrangements:                 0.7619
P(rho >= 0.7381, the corrected n=8 two-sided floor): 365/576 = 0.6337
P(rho >= 0.70):                                       415/576 = 0.7205
```

**A pure 4-vs-4 group difference with no within-cluster rank information
clears the corrected primary gate 63.4% of the time -- worse than the 47%
that got the SCLC source arm dropped.** The 8-point primary rho has, in
the worst case, roughly one degree of freedom (which cluster a gene falls
in), not eight.

**Where this differs from the SCLC case, and why it does not simply repeat
the same ruling:** in SCLC, no ambient-flagged gene survived eligibility
at all, so any clustering there was incidental to eligibility, not to
ambient risk. Here, the high cluster **is** the four ambient-flagged/
ambient-related genes and the low cluster **is** the four clean controls --
the split itself is evidence of *something* tracking ambient risk. What it
is NOT is evidence of a *monotone* association across all eight genes,
which is what the registered primary test's rho actually claims. This
distinction does not resolve the problem; it changes what the honest
conclusion can say.

**2. Data-correction to Amendment 2's CPU precheck table: the four values
are confirmed correct by exact ensembl_id lookup; a proposed correction to
one of them was itself a gene-symbol mismatch.** Amendment 2 listed
`S100PBP` at ambient_risk 0.744492. A cross-check flagged this value as
absent from `ambient_risk_all_genes.csv` and proposed 0.002612 (from a row
matched on the symbol `PEBP1`) as the true value, which would move the
split to 5-vs-3. **Resolved by ensembl_id, not symbol, per the standing
"do not substitute a similar-looking quantity" discipline this whole
amendment chain has enforced elsewhere:**

```
ENSG00000116497 -> gene=S100PBP  ambient_risk=0.744492   (the panel's registered ensembl_id for S100PBP)
ENSG00000089220 -> gene=PEBP1    ambient_risk=0.002612   (an unrelated gene, "Phosphatidylethanolamine-binding protein 1")
```

`ambient_risk_all_genes.csv` has zero duplicate ensembl_ids and zero
duplicate gene symbols -- there is no ambiguity once looked up by the
panel's own registered ensembl_id (`s100_gene_panel_20260922.json`).
**Amendment 2's original value (0.744492) and 4-vs-4 split both stand.**
`PEBP1` is a different gene that happens to share the substring "PBP"
with `S100PBP`'s symbol -- a symbol-based lookup risk, not an ensembl-id
one. No code in this analysis layer looks up ambient risk by symbol
substring; this was caught during manual cross-checking, not in code, and
is recorded here so the same substring confusion is not repeated.

**3. One pre-existing worry is tested and refuted, recorded rather than
silently dropped:** whether the four low-ambient-risk genes (as low as
1.1e-6) might sit below `ambient_risk_all_genes.csv`'s detection floor,
making their mutual ranking arbitrary the same way `S100A7`/`S100A12`
were excluded from the panel. They do not -- `S100A7` and `S100A12` are
**absent from the table entirely** ("below the detection floor" meant "not
measured," not "measured as small"), while all four low-cluster genes
here have real, present, distinct measurements. This concern is tested and
closed, not merely unraised.

**4. Three diagnostics registered now, before any result exists, all
reported always (not only on request):**

- **The exact 4-vs-4 group-separation test** (`ambient_stats.
  exact_group_separation_test()`): exact two-sided Mann-Whitney U on `Q`
  between the ambient-high and ambient-low groups -- the test the data's
  actual structure supports. With 4-vs-4 and no ties, the minimum
  attainable two-sided p (complete separation) is 2 / C(8,4) = 2/70 =
  0.02857 (verified against `scipy.stats.mannwhitneyu(method="exact")`
  directly). A significant result here says "Q is higher in the
  ambient-high group," which is a real and reportable finding -- it is
  a different claim from "Q rises monotonically with ambient risk."
- **Within-cluster Spearman, computed separately for each cluster**
  (`ambient_stats.within_cluster_spearman()`): measures directly whether
  `Q` carries rank information beyond group membership. Reported for both
  clusters unconditionally.
- **Pre-declared interpretation rule, to be applied when the real result
  exists:** if both within-cluster rhos are near zero, the result must be
  reported as **a two-group contrast between ambient-high and
  ambient-low genes after detection/rank matching** -- NOT as a monotone
  ambient-risk association -- regardless of what the primary 8-point rho
  says. **No numeric "near zero" cutoff is fixed by this amendment**; that
  judgment is applied when the real within-cluster rhos exist, by whoever
  writes the final interpretation, using the raw numbers this diagnostic
  reports -- inventing a threshold now, before either cluster's real rho
  is known, would repeat exactly the mistake this amendment chain has
  spent the day correcting (assigning a number before it can be checked
  against anything real).

**Order of work, unchanged: nothing above authorizes GPU.** This CPU-only
diagnostic layer is complete and tested; the decision on whether/how to
proceed with Module A rests with the human.

## Amendment 5 -- 2026-09-23 (human ruling, s100-isp-execution-20260922)

Closes two open items from Amendment 4: confirms the `S100PBP` correction
was the human's own error (not mine), and replaces Amendment 4's
open "no numeric near-zero cutoff is fixed" note with a structural
finding that makes fixing one unnecessary. Amendments 1-4 stand unedited.

**1. `S100PBP` vs `PEBP1`, resolved: the human's proposed correction was
the error, confirmed by the human.** The gene-symbol set used to look up
the "high cluster" values included `PBP` and `PEBP1` but never
`S100PBP` itself -- eighteen rows in `ambient_risk_all_genes.csv` contain
"PBP" as a substring (`GTPBP2`, `TOPBP1`, `HSPBP1`, `TAPBP`, among
others), so a substring search over symbols was never safe on this table.
A second, independent signal was available and unused: `PEBP1`'s
`detect_frac` (0.436) is wildly out of line with the other three
high-cluster genes' (0.012-0.067), which alone should have flagged the
wrong row. **Amendment 4's original values and the 4-vs-4 split stand,
unchanged, confirmed by both parties independently.** Standing rule,
now confirmed twice on this exact panel file: **resolve gene identity by
`ensembl_id`, registered in `s100_gene_panel_20260922.json`, never by
symbol or symbol substring.**

**2. Amendment 4 declined to fix a numeric "near zero" cutoff for
within-cluster Spearman rho, pending the real numbers. That caution was
under-stated: no cutoff is possible in principle, not just premature.**
At the real within-cluster group size (n=4), the exact null distribution
has only `4! = 24` permutations; its minimum attainable two-sided p,
reached ONLY by a perfect correlation (`rho = +-1.0`), is `2/24 = 0.08333`
(verified directly by enumeration -- the eleven attainable rho values at
n=4 are `0, +-0.2, +-0.4, +-0.6, +-0.8, +-1.0`, and every one of them
except `+-1.0` has an even larger p). **No within-cluster rho, however
extreme, can ever be distinguished from chance at p <= 0.05 at this n.** A
numeric interpretation threshold on this quantity would therefore be
theatre, not rigor -- it cannot rescue the group-separation result (nothing
here can ever reach significance) and it cannot condemn it either (a rho
of 0 is not more "real" than a rho of 0.8 at this n; both merely fail to
clear a floor nothing can clear).

**Ruling, replacing Amendment 4's open item: `within_cluster_spearman()`
is DESCRIPTIVE ONLY, permanently, not pending a threshold.** It reports
rho, its own exact two-sided p, and a fixed disclaimer string, for a
reader to see the shape of the data -- it must never gate an
interpretation. **The claim that can pass or fail rests entirely on
`exact_group_separation_test()`** (the 4-vs-8-vs-C(8,4)=70 Mann-Whitney U,
genuinely reachable at p <= 0.05 -- minimum attainable p 2/70 = 0.02857 at
complete separation). The headline n=8 primary rho is reported with its
own two-sided p AND the pre-declared note (Amendment 4, item 1) that it is
dominated by the group split: a pure 4-vs-4 separation with zero
within-cluster rank information clears the corrected primary gate 63.4%
of the time. Implemented: `ambient_stats.WITHIN_CLUSTER_N4_MIN_ATTAINABLE_P
= 2/24`; `within_cluster_spearman()` now returns, per group,
`{rho, p_exact, n, min_attainable_p_note}` rather than a bare rho.

**Every number in this analysis-layer chain (Amendments 1-5) is now
derived from an exact enumeration, a registered design-doc line, or a
panel-file lookup by ensembl_id -- none is negotiable after the fact, and
none was invented ahead of a real result.** GPU remains held; nothing in
this amendment authorizes it.
