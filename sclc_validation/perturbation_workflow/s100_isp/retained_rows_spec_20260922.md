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
