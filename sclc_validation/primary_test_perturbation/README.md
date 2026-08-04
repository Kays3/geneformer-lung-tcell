# Primary SCLC test-set perturbation analysis

This directory is the confirmatory analysis layer for the SCLC/LUAD/normal
Geneformer V2 experiment. It is deliberately separate from the exploratory
targeted-panel and whole-cohort analyses.

## Primary estimand

The primary results use only donor-held-out **test cells**. Disease-state
reference embeddings are calculated from training cells only. For each test
cell, Geneformer's native deletion and overexpression interventions are scored
as changes in cosine similarity to the training-derived disease centroids.

The inferential replicate is the donor, not the cell. Cell-level Geneformer
statistics are retained for auditability, but a candidate is not considered
robust without adequate detections and donor-level sign consistency.

## Build sequence

The all-gene computation is run by the separate two-node launcher. After all
test-set shards have completed and remote overexpression outputs have been
synchronized:

```bash
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh sync
./sclc_validation/primary_test_perturbation/run_primary_analysis.sh
```

The report builder is resumable and can be run directly:

```bash
python sclc_validation/primary_test_perturbation/scripts/build_primary_report.py
python sclc_validation/primary_test_perturbation/scripts/build_notebook.py
```

Override the compute artifact location with `SCLC_PERTURBATION_ROOT` and the
Python executable with `PYTHON_BIN`.

## Output contract

- `reports/primary_test_perturbation_report.md`: human-readable decision report
- `reports/primary_test_perturbation_report.html`: portable HTML report
- `tables/primary_arm_summary.csv`: delete/overexpression comparison status
- `tables/primary_concordant_hits.csv`: validated two-arm candidates
- `tables/coverage_audit.csv`: missing files, columns, and detection coverage
- `primary_test_perturbation.ipynb`: reproducible review notebook
- `analysis_manifest.json`: paths, filters, and completeness state

Missing computation artifacts are reported as `PENDING`; they are never
silently converted to empty results.

Performance benchmark conclusions for the DGX Spark GB10 environment are
recorded in [`PERFORMANCE_NOTE.md`](PERFORMANCE_NOTE.md). In this setup,
changing `forward_batch_size` and `nproc` did not materially improve runtime;
future optimization should focus on shard scheduling and additional GPUs.

The latest dated execution summary and post-run analysis plan are in
[`STATUS_2026-08-04.md`](STATUS_2026-08-04.md).
