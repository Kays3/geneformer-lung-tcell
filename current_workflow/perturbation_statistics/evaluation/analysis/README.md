# Statistical and technical evaluation

## Objective

Determine whether ranked goal shifts are calibrated, donor-stable, and robust
to analytical choices. Cell-level significance alone is insufficient because
cells from the same donor are not independent biological replicates.

## Priority 1: donor stability

For every comparison and gene:

- estimate the shift separately for each donor;
- report donor count, cells per donor, median donor effect, effect range, and
  sign concordance;
- bootstrap donors for confidence intervals;
- calculate leave-one-donor-out rank and effect stability; and
- flag results whose significance or direction depends on one donor.

Suggested provisional gates:

- the same effect direction in at least 70% of evaluable donors;
- a bootstrap interval that excludes zero for the donor-level aggregate; and
- top-50 rank retention in at least 70% of leave-one-donor-out runs.

LUSC has only five held-out donors, so its intervals and stability estimates
must be reported explicitly rather than compared mechanically with the larger
LUAD and normal donor groups.

## Priority 2: null calibration and multiplicity

- Permute disease labels at the donor level, not the cell level.
- Add expression-frequency- and detection-count-matched random-gene controls.
- Inspect null p-value uniformity and empirical false-discovery rates.
- Compare within-table Benjamini-Hochberg FDR with correction across all six
  directional tables.
- Use known irrelevant or synthetic control genes where technically possible.

Key outputs should include empirical type-I error, observed-versus-null hit
counts, p-value calibration plots, and false-discovery estimates at each
proposed threshold.

## Priority 3: effect-size robustness

- Report standardized effects and confidence intervals alongside FDR.
- Repeat rankings at `N_Detections` thresholds of 10, 25, 50, and 100.
- Measure top-20, top-50, and top-100 overlap across thresholds.
- Inspect outlier concentration for unusually large shifts, including
  LUSC → LUAD `MMP12`.
- Separate statistically precise but biologically tiny shifts from effects
  large enough to change the embedding materially.

## Priority 4: centroid and model sensitivity

- Recompute goal centroids while leaving out the evaluated donor.
- Compare mean centroids with robust alternatives such as medoids or trimmed
  centroids.
- Repeat the leading rankings using another model checkpoint or embedding
  layer when available.
- Quantify agreement using effect correlation, rank correlation, top-k overlap,
  and direction concordance.

## Priority 5: reciprocal and directional structure

- Compare A → B with B → A effects for the same genes.
- Classify genes as reciprocal, source-specific, goal-specific, or unstable.
- Test whether deletion removes a source-state feature versus inducing a
  coherent goal-state program.
- Compare goal shifts with source-versus-goal differential expression without
  assuming that the two quantities should have the same sign.

## Minimum analysis deliverables

- `donor_gene_effects.csv`
- `leave_one_donor_out_stability.csv`
- `null_calibration_summary.csv`
- `coverage_sensitivity.csv`
- `centroid_model_sensitivity.csv`
- `reciprocal_shift_summary.csv`
- figures showing donor effects, rank stability, p-value calibration, and
  sensitivity comparisons

Write future outputs under `../results/analysis/` and preserve scripts or
notebooks needed to reproduce them.
