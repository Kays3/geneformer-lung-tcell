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
