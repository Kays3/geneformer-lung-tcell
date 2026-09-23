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
