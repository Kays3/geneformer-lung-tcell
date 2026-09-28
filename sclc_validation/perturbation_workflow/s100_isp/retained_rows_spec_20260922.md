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

## Amendment 6 -- 2026-09-28 (human approval of the LUAD run, s100-isp-execution-20260922)

Registered before any Module A output exists: no S100 LUAD perturbation
has been run, no `Q`, `E`, or shift value exists for any panel gene.
Amendments 1-5 stand unedited. The human approved the LUAD arm on
2026-09-28 on the condition that the primary analysis is re-declared as
the exact two-group test. Budget: tiny no-op inference plus 2.98 GPU-hours,
a hard stop. Module B is not approved. The SCLC source arm stays dropped.

**1. The primary test is now the exact two-group test.** For each LUAD
contrast separately (`LUAD -> normal` and `LUAD -> SCLC`), compare `Q`
between two fixed groups with an exact two-sided Mann-Whitney U test
(`ambient_stats.exact_group_separation_test()`). The groups are fixed by
the ambient-risk table in Amendment 2 and identified by the `ensembl_id`
registered in `s100_gene_panel_20260922.json` (sha256 `06ee2a51...`):

| group | genes (ensembl_id) | ambient_risk |
|---|---|---|
| ambient-high | S100A2 (ENSG00000196754), S100B (ENSG00000160307), S100A13 (ENSG00000189171), S100PBP (ENSG00000116497) | 0.744 - 0.955 |
| ambient-low | S100A4 (ENSG00000196154), S100A6 (ENSG00000197956), S100A10 (ENSG00000197747), S100A11 (ENSG00000163191) | 1.1e-6 - 3.1e-4 |

The p-value is an explicit label permutation: all C(8,4) = 70 ways of
assigning the eight observed `Q` values to a 4-vs-4 split are enumerated,
U is computed on midranks, and p is the fraction of splits at least as far
from the null mean (8) as the observed U. This is exact with or without
ties.

**2. The exact p floor.** With no ties, the minimum attainable two-sided p
is 2/70 = 0.02857, reached only by complete separation. One inversion
already gives 4/70 = 0.05714, which fails p <= 0.05. In practice this test
passes only when every ambient-high gene's `Q` sits above every
ambient-low gene's `Q`.

**3. What counts as a positive, per contrast:** exact p <= 0.05 **and** the
ambient-high group above the ambient-low group. A separation significant
in the opposite direction (p = 2/70 with the clean genes above) is reported
as "significant in the direction opposite the prediction", never as a
positive. If any of the eight genes has no estimable `Q` in a contrast,
that contrast's primary result is `not_estimable`; it is not re-run on
fewer genes. The two contrasts are each reported. Carried over unchanged
from the design ("all six contrasts within a module are reported ...
multiple contrasts are descriptive unless a later, separately registered
family-level correction is added"): no family-wise correction is applied,
and none can be added after results exist. Note for whoever reads this
later: a Bonferroni correction across the two contrasts (alpha 0.025 each)
would make the test unreachable, since 2/70 = 0.02857 > 0.025.

**4. The monotone statistic is demoted to descriptive.** The eight-point
Spearman rho of `Q` against ambient risk (`primary_test()`) and its
leave-one-gene-out check (`leave_one_gene_out_check()`) are still computed
and always reported with their exact two-sided p, but they no longer gate
any conclusion. Reason, from Amendment 4: a pure 4-vs-4 group difference
with no within-group rank information clears the old rho gate 365/576 =
63.4% of the time. `within_cluster_spearman()` stays descriptive, as ruled
in Amendment 5. The claim this run can support is "`Q` is higher in the four
ambient-high genes than in the four clean genes after detection and
token-rank matching", not "`Q` rises with ambient risk".

**5. Knock-on re-declarations, proposed by Kevin, for Stanley and the
human to confirm or amend before any output exists.** The design ties
three further rules to the rho gate. Each is re-stated for the new primary
so no gate is left undefined:

- **Control-stability gate.** Same resampling as registered (the 120
  leave-one-control-out recomputations; 10,000 within-stratum bootstrap
  redraws, seed `20260922`), now recomputing `Q` and the group-separation
  test each time. The primary positive stands only if (a) the unresampled
  test is positive, (b) every leave-one-control-out recomputation is still
  positive, and (c) at least 95% of bootstrap redraws are positive.
  Otherwise the contrast is "control-draw-sensitive / open".
- **Raw diagnostic.** The same exact two-group test run on raw `E` instead
  of `Q`, so the "detection/token-frequency pattern" row compares like with
  like. The raw `E` Spearman stays descriptive.
- **Interpretation table.** Wherever the design says "the detection-adjusted
  `Q` rho gate passes/fails", read "the stable group-separation primary is
  positive/not positive". All other conditions in each row are unchanged
  (including `S100A2` being non-dose-responsive or donor-inconsistent for
  the ambient-compatible row).

**6. Correction to Amendment 4, item 4, found while registering this
amendment.** Amendment 4 said the group-separation test was "verified
against `scipy.stats.mannwhitneyu(method="exact")` directly". That is true
only without ties. scipy's exact mode uses the no-ties null even when ties
are present: for high `Q` = (0.90, 0.95, 0.97, 1.0) and low `Q` = (0.1, 0.2,
0.3, 1.0), scipy gives p = 0.343 where the true 70-split permutation p is
16/70 = 0.229. `Q` is a percentile against 20 controls, so ties across
genes can happen (two genes that both beat every control). The function
now enumerates the 70 splits itself; tests cover complete separation
(2/70), one inversion (4/70), the reversed direction (not positive), and
the tie case (16/70). The error was Kevin's.

**Order of work from here:** re-pin the runner hash on thinkstation1 after
PR #32; tiny no-op inference; CPU completion (gate 3); manifests to Stanley.
No GPU until Stanley signs. Nothing in this amendment starts GPU work.

## Amendment 7 -- 2026-09-28 (conditions from Michael and Stanley; gate-3 record)

Appended before any `Q`, `E` or perturbation shift exists for any gene.
The only S100 outputs on thinkstation1 are the tiny no-op (S100A2, LUAD,
shift exactly 0.0) and the CPU preflight tables. Amendments 1-6 stand
unedited; Stanley signed off Amendment 6 (sign-off 1, PASS).

**1. Ties rule, stated in the spec (Michael's condition).** Amendment 6
item 1 already registers it; restated here so it cannot be missed. The
eight `Q` values of a contrast are ranked together; tied values get the
average of the ranks they span (midranks). U for the ambient-high group is
the sum of its midranks minus 10. The p-value is the exact permutation p:
the fraction of all C(8,4) = 70 assignments of the eight observed values to
a 4-vs-4 split whose U is at least as far from 8 as the observed U. No
normal approximation, no tie correction formula, no library "exact" mode
that assumes no ties.

**2. Leave-one-gene-out under the new primary (Stanley's finding 2).** It
stays on rho and is descriptive only, as Amendment 6 item 4 says. No
leave-one-gene-out version of the group test is registered, because it
could never pass: with one gene removed the design is 4 vs 3, and the
minimum attainable exact two-sided p is 2/C(7,3) = 2/35 = 0.0571 > 0.05
(verified by enumeration). Any such number, if reported, carries that
floor next to it.

**3. Control drawing implements the design's "common eligible ...
non-S100" wording (commit `3733aa89`).** A control must pass the run's own
LUAD eligibility gate (>= 50 token-positive cells from >= 3 donors) and
must not be any S100-family gene (32 ensembl ids by symbol prefix in the
gc104M name dictionary), not only the twelve panel genes. Matching
tolerances, the all-members rule, the seed and the 20-control floor are
unchanged.

**4. Gate-3 record: as registered, the primary is not estimable.** Built on
2026-09-28 from LUAD-only statistics (Amendment 2 ruling):

| stratum | candidates | status |
|---|---:|---|
| clean_high (S100A4, S100A6, S100A10, S100A11) | 0 | not estimable: members' rank percentiles span 70.2-90.6 (20.3 points); a control within 5 points of every member needs a span <= 10 |
| low_p_a16 (S100P, S100A16) | 27 | controls drawn, but both members fail LUAD eligibility |
| low_a2_b (S100A2, S100B) | 43 | eligible |
| singleton_a13 (S100A13) | 141 | eligible |
| singleton_a8 (S100A8) | 69 | eligible |
| low_a9_pbp (S100A9, S100PBP) | 0 | not estimable: detection differs by 1.34 log2 units (max 1.0), rank by 48.8 points |

The zero counts do not come from the item-3 filters: with both filters
off they are still 0 and 0. Under Amendment 6 item 3, missing `Q` for any
of the eight primary genes makes the contrast `not_estimable`; five of the
eight (all four ambient-low genes and S100PBP) have no stratum, so both
contrasts' primary result is `not_estimable` as registered. The design
forbids widening or splitting a stratum after the preflight, and nothing
has been widened or split.

**No ruling is made here.** Whether to stop, or to re-register matching
(for example, each primary gene matched on its own with the same
tolerances), is the human's decision. That decision, if it changes the
matching rule, must be a further dated amendment registered before any GPU
run. No GPU beyond the tiny no-op has been spent.

## Amendment 8 -- 2026-09-28 (human ruling: per-gene matching, 3.4 GPU-h; s100-isp-execution-20260922)

**This amendment overrides design lines 171-172 ("it is not widened,
split, or given private controls after the preflight has been viewed"),
after the gate-3 preflight was viewed, because two of the six registered
strata are infeasible by arithmetic on the LUAD-only statistics (Amendment
7, item 4); no effect, `E`, `Q` or perturbation-shift data existed for any
gene when it was written.** The human ruled at about 12:37 JST on
2026-09-28 (relayed 12:38) and reconfirmed at 12:55 JST after being told
this reverses lines 171-172. Amendments 1-7 stand unedited.

**Root cause (Stanley's question).** The six strata were grouped on pooled
detection fraction only, before any token rank existed (22 Sep preflight:
"source-specific token ranks require the held-out dataset, absent
locally"). On the pooled table S100A9 and S100PBP differ by 0.07 log2
units and clean_high spans 0.40; on LUAD-only statistics they differ by
1.34 and clean_high's ranks span 20.3 percentile points. Amendment 2 moved
both matching axes to LUAD-only and nobody, Kevin included, re-checked that
the strata were still satisfiable.

**1. Matching rule.** Each of the eight primary genes gets its own 20
controls. A candidate qualifies for gene g if it is within 0.5 log2
detection units and 5 rank-percentile points of g (the registered
tolerances), passes the run's LUAD eligibility gate (>= 50 token-positive
cells from >= 3 donors), and is not an S100-family gene (32 ensembl ids)
or a panel gene. Statistics: the frozen LUAD tables from gate 3
(`luad_gene_stats.csv` sha256 `cbad7c7b...`, `luad_isp_eligibility_all_genes.csv`
sha256 `552ee81d...`). Drawn once, without replacement, by
`matched_controls.build_matched_control_table()` with one set per gene; the
generator is seeded from (20260922, gene symbol); qualifying candidates are
sorted by ensembl_id before the draw. Never re-sampled. A control may serve
more than one gene; it is run once and its `E` is reused.

| gene | group | candidates | drawn |
|---|---|---:|---:|
| S100A2 | ambient-high | 152 | 20 |
| S100B | ambient-high | 129 | 20 |
| S100A13 | ambient-high | 141 | 20 |
| S100PBP | ambient-high | 298 | 20 |
| S100A4 | ambient-low | 51 | 20 |
| S100A6 | ambient-low | 34 | 20 |
| S100A10 | ambient-low | 71 | 20 |
| S100A11 | ambient-low | 45 | 20 |

160 draws, 142 distinct controls, 18 shared by more than one gene. All 142
confirmed eligible by the runner's own `paired_eligible_dataset()`; paired
cell lists built. Committed before any GPU run: control sets
`matched_control_sets_pergene.json` sha256 `607ac6c9...`; run panel
`s100_luad_run_panel_pergene_20260928.json` (154 genes: 12 registered +
142 controls) sha256 `8a0b3668...`; eligibility manifest sha256
`b8d989f3...`. The registered panel file is unchanged (`06ee2a51...`).

**1a. Realised control overlap, S100A6's small pool (a limitation,
reported with the result).** Candidate pools overlap heavily (S100A6 and
S100A10 share 32 of S100A6's 34 candidates; S100A2 and S100A13 69%; S100B
and S100A13 59%). The drawn sets share 18 controls, all within a group,
none across: S100A2 x S100B 2, S100A2 x S100A13 4, S100B x S100A13 1,
S100A6 x S100A10 6, S100A6 x S100A11 2, S100A10 x S100A11 3 (of 20 each;
`pairwise_control_overlap.csv`, sha256 `91ea3afd...`). Q values within a
group are therefore correlated, which weakens the exchangeability the
permutation null assumes; the p-value is reported with that caveat.
S100A6 draws 20 of only 34 candidates, so its reference set is close to
exhaustive rather than a sample.

**1b. What a positive can claim.** Only: "sparse ambient-high genes beat
their own detection-matched controls by more than dense clean genes beat
theirs." A claim that ambient-high genes simply have larger effects is not
supported, for two separate reasons. (1) The two groups sit at opposite
ends of detection: Spearman -0.952 between LUAD detection fraction and
ambient risk across these eight genes (ambient-high genes detected in
1.1-4.6% of LUAD cells, clean genes in 63-88%; the design's own figure is
-0.9758 on the pooled ten-gene table). (2) Genes detected in fewer cells
tend to show larger nominal perturbation effects. Together, raw effects
would favour the ambient-high group from detection alone, which is why `Q`
is matched to each gene's own detection and raw `E` stays a diagnostic.

**2. Anchors and ineligible genes.** S100A8 and S100A9 run as panel genes
without their own control sets, so they have `E` and dose-response results
but no `Q`; the design's secondary circular `Q` analysis for the anchors
is `not_estimable`. S100P and S100A16 fail LUAD eligibility and are not
run; their rows are retained with status and reason. The gate-3 strata
table (Amendment 7) is not used.

**2a. No-op (human ruling, 12:55 JST).** The design's "unperturbed/no-op
inference on the same cell lists" runs on the 10 eligible S100 panel genes
(S100A2, S100B, S100A13, S100PBP, S100A4, S100A6, S100A10, S100A11, S100A8,
S100A9), on the same paired cell lists as their delete and overexpress
arms. Matched-control rows are not given a no-op arm; their `no_op_status`
is `not_run_by_design` in the retained-row output
(`retained_rows.NO_OP_NOT_RUN_BY_DESIGN`), never `not_run` or
`no_op_failed`.

**3. Stability gate, per gene (Amendment 6 item 5, confirmed by the
human).** Leave-one-control-out: 160 recomputations, each dropping one
control from one gene's set and recomputing that gene's `Q` and the
primary test. Bootstrap: 10,000 redraws, seed 20260922, of 20 controls with
replacement within every gene's set simultaneously (a shared control is
resampled independently in each set). The primary positive stands only if
the unresampled test is positive, all 160 leave-one-control-out
recomputations are positive, and at least 95% of bootstrap redraws are
positive. The design's line-188 note about shared controls correlating `Q`
within a stratum no longer applies as written: sharing is now partial (18
controls); the exact permutation p still conditions on the observed `Q`
vector, and the stability gate is what tests sensitivity to the draw.

**4. Budget and run.** 3.4 GPU-hours, a hard stop, not a target (human
ruling). Projection 3.15 GPU-h: 3.05 h for 152 genes (10 eligible panel
genes + 142 controls) x 2 operations (delete, overexpress), LUAD only,
40,144 capped cells per operation, plus 0.10 h for the no-op arm on the
10 panel genes (2,615 capped cells; a no-op arm is two forward passes per
cell, the same work as a perturbation arm);
per arm = 6.7 s + 0.1113 s/cell, fitted from the 316M bf16 LUAD
overexpress arms of run `316m_bf16` (50 genes at 300 cells, median 39.7 s,
range 39.0-41.5). Sensitivity: the one paired 316M calibration point has
delete 12% slower than overexpress (46.1 s vs 41.2 s), which gives about
3.33 h. During the run, cumulative GPU time is read from the completion
markers; if elapsed plus the remaining projection (rescaled by the observed
rate) exceeds 3.4 h, the run is stopped and reported. No gene, control,
source or direction is dropped to fit, and a partial run is not analysed.
Runner: main at `de8cda2` (sha256 `974535b7...`), `--sources luad`,
verified on thinkstation1 by fast-forward pull immediately before GPU.

**5. Provenance.** Build script `preflight_s100_luad_pergene_20260928.py`
(sha256 `5fe5715e...`) is committed with this amendment; it reuses the
frozen gate-3 tables and was run on thinkstation1 with every root set
explicitly. Its log and outputs are in
`bf16_bench/runs/s100_luad_preflight_pergene_20260928/` on thinkstation1.


## Amendment 9 -- 2026-09-28 (budget stop and resume at 4.2 GPU-h; s100-isp-execution-20260922)

Appended after the first launch was stopped by the registered budget rule and
before the resume. **No output has been analysed.** Effect data now exists:
the 7 completed delete arms wrote raw pickles and the no-op arm wrote 2 stats
tables. Stanley confirmed by access time that none of these 9 files had been
opened as of 04:16Z (relatime filesystem; access time equals write time for
all 9). Amendments 1-8 stand unedited.

**1. The stop.** Launched 04:02:14Z (13:02:14 JST) on runner `974535b7` at
`de8cda2`. The budget watcher stopped the run at 04:14:20Z (13:14:20 JST).
`s100_luad_20260928.BUDGET_STOP` reads "projected 3.473 h > 3.4 at elapsed
0.202 h". **GPU used: 0.202 h.** Completed: no-op on the 10 eligible panel
genes (plus not-estimable records for S100P and S100A16); delete on S100A4,
S100A6, S100A10, S100A11, S100A2, S100B and S100A13 (plus the two
not-estimable records). No overexpress arm and no control arm ran. The
S100A8 delete arm was in progress and left no files. For the first 3 min 6 s
(04:02:14-04:05:20Z) the run had no working budget stop, because watcher v1
was given the SSH shell's process group (Kevin's error). Watcher v2, on the
run's real group 1482351, made the stop.

**2. Cause: the projection, not the rule.** Amendment 8 projected each arm as
6.7 s + 0.1113 s/cell. Every 316M calibration run had exactly 300 cells, so
the fixed part could not be fitted; 6.7 s came from the tiny no-op, whose
setup is far smaller. Observed: each delete arm spends 18.5-20.3 s in library
setup before the first forward pass (new perturber, model load, dataset
filter and sort, repeated for every gene), whatever its cell count.

**3. Revised projection method.** Per delete or overexpress arm: 37.88 s +
0.0394 s/cell. This is the least-squares fit on the 7 completed delete arms
(residual sd 4.5 s, intercept standard error 3.5 s), reproduced independently
by Kevin and Stanley. Overexpress is projected at the delete rate (the
conservative choice; the old 316M overexpress arms ran about 17% faster).
The total, including the 0.202 h spent, is 3.84 h (overexpress 17% faster)
to 4.19 h (same speed). With the fixed cost 2 standard errors either side and
overexpress at delete speed, it is 3.62-4.76 h. Stanley's wider range is
3.14-4.95 h. **A 4.2 h budget sits at the central estimate, so a second
budget stop is possible.** The fit rests on 3 small-cell genes, and 54 of
the 152 run genes have 150 or fewer capped cells.

**4. Human ruling, 13:16 JST.** Budget raised to **4.2 GPU-hours total, a
hard stop, with the 0.202 h already spent counted toward it.** Resume the
same run. No runner change.

**5. Resume.** Same launcher, runner (`974535b7`), run panel (`8a0b3668`),
inputs (re-hashed unchanged after the stop), `--sources luad`, run tags and
seed. Completed arms are kept and not rerun: the runner skips any gene with
a completion marker, and they came from the identical runner, config and
inputs. Before the resume, the first launch's `run_config.json` for both run
tags, `started_utc`, log and `BUDGET_STOP` were copied with a
`first_launch` or `_1` suffix. The original `BUDGET_STOP` is moved to
`BUDGET_STOP_1`, so a new `BUDGET_STOP` can only come from the resumed run.
Watcher v3 (`watch_s100_budget_v3.py`, sha256 `ae87c5f5...`) counts the
0.202 h already spent plus wall clock since the resume. It projects the
remaining delete and overexpress arms with the section 3 method, scaled by
the observed rate when that is slower. It stops the run's process group at
4.2 h, and refuses to start unless its process group contains the run.

**6. Unchanged.** If the run is stopped again, nothing from it is analysed,
exactly as for this stop.

## Amendment 10 -- 2026-09-28 (human ruling 13:20 JST: no GPU-hour cap today; s100-isp-execution-20260922)

Appended before the resume; no new output exists since Amendment 9, and none
has been analysed. Amendments 1-9 stand unedited.

**1. Budget.** The human ruled at 13:20 JST that runs on 2026-09-28 have
**no GPU-hour cap**. This supersedes Amendment 9's 4.2 GPU-h (ruled 13:16
JST). It also supersedes a 5.0 GPU-h figure the human set a few minutes
earlier, relayed by Michael and superseded within minutes. The projection
method and range in Amendment 9 section 3 stay on the record as the
expected cost: central 3.84-4.19 GPU-h including the 0.202 h already
spent, ranges 3.62-4.76 h (Kevin, intercept +/-2 standard errors) and
3.14-4.95 h (Stanley). Kevin's first fit on the same 7 delete arms (34 s +
0.047 s/cell) gives about 2 h per operation, the same as the OLS fit.

**2. Monitor instead of a kill switch.** `watch_s100_monitor_v4.py` (sha256
`99edb6f3...`), derived from watcher v3. It logs elapsed GPU time (0.202 h
plus wall clock since the resume) and the Amendment 9 projection every 60 s.
It never signals the run. If no arm completes for 30 minutes it writes
`s100_luad_20260928.STALL_ALERT`, and Kevin relays that to Michael. The run
is not killed on a stall; that decision stays with Michael and the human.

**3. Scope unchanged.** The same 154-gene LUAD run: the same launcher
(`0a7813b2`), runner (`974535b7`), run panel (`8a0b3668`), `--sources luad`,
and no-op on the 10 eligible panel genes only. Module B (CAR-T) is not
approved. No arm, gene, source or direction is added because the budget is
now uncapped.

## Amendment 11 -- 2026-09-28 (monitor fix from Stanley's review; s100-isp-execution-20260922)

Appended before the resume; no new output exists. Amendments 1-10 stand
unedited. Stanley found that watcher v3 would have stopped the run during the
runner's CPU-only stats phase, after every GPU arm had finished, and so
discarded a complete run. The same pattern would have raised a false stall
alert in monitor v4. Monitor v4 now logs "all GPU arms complete" once every
delete/overexpress arm has a completion marker, and turns its stall alarm off
for the stats phase that follows. It still never signals the run. The version
used for the resume is `watch_s100_monitor_v4.py`, sha256 `efac8bdd...`, which
replaces the `99edb6f3...` named in Amendment 10. Watcher v3 (4.2 h kill) is
not used: the cap was removed at 13:20 JST.
