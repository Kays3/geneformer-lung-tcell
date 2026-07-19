# Perturbation evaluation workspace

## Purpose

This workspace is the starting point for deciding whether the completed
perturbation screen can advance from a hypothesis-generating ranking to a
donor-stable, T-cell-intrinsic biological result.

The primary decision question is:

> Do the leading goal shifts reproduce across donors and T-cell subsets after
> controlling for ambient RNA, doublets, coverage, and centroid construction?

## Directory map

- [Analysis evaluation](analysis/README.md): statistical calibration,
  donor-level robustness, null models, rank stability, and representation
  sensitivity.
- [Biological evaluation](biology/README.md): T-cell specificity,
  contamination sensitivity, pathway interpretation, external replication,
  and experimental validation.
- [Evaluation tracker](evaluation_tracker.csv): work items, provisional gates,
  required inputs, and output locations.
- [Results workspace](results/README.md): conventions for future evaluation
  artifacts. Do not overwrite the completed primary analysis outputs.

## Recommended execution order

1. Build donor-level gene-effect tables and leave-one-donor-out rankings.
2. Quantify epithelial/alveolar burden and doublet risk, then repeat rankings
   after decontamination and exclusion sensitivity analyses.
3. Calibrate the null with donor-level label permutations and matched random
   genes.
4. Stratify stable signals by T-cell subtype and compare them with
   source-versus-goal differential expression.
5. Consolidate pathway leading-edge modules and test their donor stability.
6. Replicate the final shortlist in an independent cohort or orthogonal
   perturbation system.

## Promotion rule

A candidate should not be called a biological target unless it has:

- a consistent direction across donors;
- a meaningful, uncertainty-bounded effect size;
- stability after ambient-RNA and doublet sensitivity checks;
- expression and interpretation compatible with T-cell biology;
- robustness to centroid and model choices; and
- replication in an independent dataset or orthogonal experiment.

Thresholds in the tracker are provisional planning gates. They should be
finalized before looking at the corresponding validation result.
