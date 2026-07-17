# Report source and chart notes

The HTML report is built from `reports/artifact.json` with the packaged Data
Analytics portable-report builder. It is self-contained and read-only.

## Evidence map

- Target-cohort chart and table: `tables/tcell_cancer_vs_normal_summary.csv`
- Disease-specific T-cell table: `tables/stage2_tcell_cohort_summary.csv`
- Classifier chart and table: `tables/classifier_metrics.csv`
- Embedding separation table: `tables/embedding_master_table.csv`
- Donor diagnostic table: `tables/stage1_donor_leakage_summary.csv`,
  `tables/stage2_donor_leakage_summary.csv`, and
  `tables/stage2_chance_adjusted_donor_enrichment.csv`
- Visual evidence archive: `figures/embeddings/` and `figures/classifiers/`

## Quality notes

- Stage 2 was sampled to 5,000 cells per disease class, so class-level totals
  are balanced but cell-type composition is not.
- The cancer T-cell target cohort contains substantially more cells than the
  normal cohort; donor-level aggregation and resampling are required.
- Donor purity is much higher than its chance expectation in Stage 2, so
  donor-aware splits and donor-level uncertainty are mandatory.
- Saved classifier metrics and confusion matrices are treated as verified
  outputs. Notebook figures restored from prior executions were not recomputed
  during packaging.
- Step 7 placeholder and random-fallback outputs are excluded from report claims.
