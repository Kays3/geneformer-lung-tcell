# Primary donor-held-out SCLC perturbation report

**Status:** PENDING: computation artifacts are incomplete
**Generated:** 2026-08-01T09:56:37.934752+00:00

## Decision boundary

Only test-cell perturbations with training-only reference centroids are primary evidence. Training, evaluation, and pooled whole-cohort perturbations are excluded from this report's primary estimates.

## Completeness

- Required arm/comparison tables: 0/12
- Primary concordant gene-comparison hits: 0
- FDR threshold: 0.05
- Minimum detections: 25

## Arm summary

```text
        arm     comparison  status  genes_tested  qualified_primary  fdr_significant  positive_shift
     delete   sclc_to_luad missing             0                  0                0               0
     delete sclc_to_normal missing             0                  0                0               0
     delete   luad_to_sclc missing             0                  0                0               0
     delete luad_to_normal missing             0                  0                0               0
     delete normal_to_sclc missing             0                  0                0               0
     delete normal_to_luad missing             0                  0                0               0
overexpress   sclc_to_luad missing             0                  0                0               0
overexpress sclc_to_normal missing             0                  0                0               0
overexpress   luad_to_sclc missing             0                  0                0               0
overexpress luad_to_normal missing             0                  0                0               0
overexpress normal_to_sclc missing             0                  0                0               0
overexpress normal_to_luad missing             0                  0                0               0
```

## Coverage audit

```text
        arm     comparison                                                                                                                                                 path  exists  status  rows
     delete   sclc_to_luad             /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_sclc_to_luad.csv   False missing     0
     delete sclc_to_normal           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_sclc_to_normal.csv   False missing     0
     delete   luad_to_sclc             /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_luad_to_sclc.csv   False missing     0
     delete luad_to_normal           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_luad_to_normal.csv   False missing     0
     delete normal_to_sclc           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_normal_to_sclc.csv   False missing     0
     delete normal_to_luad           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_normal_to_luad.csv   False missing     0
overexpress   sclc_to_luad   /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_sclc_to_luad.csv   False missing     0
overexpress sclc_to_normal /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_sclc_to_normal.csv   False missing     0
overexpress   luad_to_sclc   /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_luad_to_sclc.csv   False missing     0
overexpress luad_to_normal /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_luad_to_normal.csv   False missing     0
overexpress normal_to_sclc /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_normal_to_sclc.csv   False missing     0
overexpress normal_to_luad /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_normal_to_luad.csv   False missing     0
```

## Interpretation

A concordant hit requires significant, adequately detected, opposite-signed deletion and overexpression shifts. This is a model-level prioritization criterion, not evidence of experimental causality. Donor-level consistency, ambient-RNA/doublet sensitivity, and independent validation remain required before biological claims.

The report is intentionally incomplete. Do not interpret missing comparisons as null effects.