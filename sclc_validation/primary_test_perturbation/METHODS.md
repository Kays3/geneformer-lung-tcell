# Primary analysis methods

## Population and split

The primary cohort is the HTAN/CELLxGENE SCLC/LUAD/normal T-cell object. Cells
are assigned to train, evaluation, and test at the donor level, with the
cross-disease patient guard documented in
`../audit/SCLC_DATA_FEASIBILITY_AUDIT.md`. No test donor contributes to model
fine-tuning or reference-centroid calculation.

## Model and intervention

The Geneformer V2 104M fine-tuned cell classifier is frozen before perturbation.
The native `InSilicoPerturber` operations are used:

- `delete`: remove one detected gene token from the rank-value encoding;
- `overexpress`: move one gene token to the front of the encoding.

These are representation-level interventions, not experimental knockout or
expression measurements. They support prioritization and model-mechanism
hypotheses, not causal biological claims by themselves.

## Primary statistic

For source state `a`, goal state `b`, and perturbation type `t`:

```text
shift(a -> b, t) = cosine(perturbed_cell_t, centroid_b_train)
                    - cosine(original_cell, centroid_b_train)
```

The primary comparison requires `Goal_end_FDR < 0.05`, positive goal shift for
the stated direction, and a predeclared minimum detection threshold. Delete and
overexpression are stored and analysed separately because Geneformer's stats
reader consumes all raw pickle files in a directory.

## Reliability checks

1. All six source-to-goal comparisons must have non-empty tables for each arm.
2. Gene and Ensembl identifiers must be unique within each comparison table.
3. Detection counts and FDR values must be present and numerically valid.
4. Concordance requires both arms to pass FDR and have opposite-signed shifts.
5. Candidate rankings must be checked for donor-level sign consistency using raw
   shard outputs; cell counts alone are insufficient.
6. Ribosomal, stress, hemoglobin, epithelial, stromal, and ambient-RNA markers
   receive a sensitivity flag and are not treated as T-cell regulators without
   independent evidence.

The full-cohort and evaluation-cell analyses are secondary sensitivity analyses
and must not be pooled into the primary estimate.
