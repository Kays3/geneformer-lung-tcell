# Curated tables

These compact tables were extracted from the completed Geneformer workflow.
They intentionally exclude expression matrices, tokenized datasets, embeddings,
checkpoints, and model weights.

## Primary review tables

- `tcell_cancer_vs_normal_summary.csv`: target cohort size by state and T-cell class.
- `stage2_tcell_cohort_summary.csv`: T-cell counts and donor coverage by disease.
- `stage2_celltype_by_disease.csv`: full Stage 2 cell-type/disease cross-tabulation.
- `classifier_metrics.csv`: held-out accuracy and macro-F1 from saved metric dictionaries.
- `stage_1_confusion_matrix.csv` and `stage_2_confusion_matrix.csv`: held-out confusion matrices.
- `embedding_master_table.csv`: training and embedding summary for both stages.
- `stage1_donor_leakage_summary.csv`, `stage2_donor_leakage_summary.csv`, and
  `stage2_chance_adjusted_donor_enrichment.csv`: donor-structure diagnostics.
- `stage2_donor_by_disease.csv`: per-donor cell counts for the balanced Stage 2 subset.

## Interpretation boundary

Rows are descriptive evidence from the sampled analysis datasets. They do not
yet contain in-silico perturbation effects or causal tumor-microenvironment
estimates. Cancer combines lung adenocarcinoma and squamous cell lung carcinoma;
COPD is retained as an inflammatory comparator and is not grouped with normal.
