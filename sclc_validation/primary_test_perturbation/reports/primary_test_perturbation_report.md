# Primary donor-held-out SCLC perturbation report

**Status:** COMPLETE
**Generated:** 2026-08-06T14:31:52.786090+00:00

## Decision boundary

Only test-cell perturbations with training-only reference centroids are primary evidence. Training, evaluation, and pooled whole-cohort perturbations are excluded from this report's primary estimates.

## Completeness

- Required arm/comparison tables: 12/12
- Primary concordant gene-comparison hits: 3586
- FDR threshold: 0.05
- Minimum detections: 25

## Arm summary

```text
        arm     comparison status  genes_tested  qualified_primary  fdr_significant  positive_shift
     delete   sclc_to_luad     ok         14415                405              954            7779
     delete sclc_to_normal     ok         14415                310              700            6537
     delete   luad_to_sclc     ok         15185                955             2881            5717
     delete luad_to_normal     ok         15185               1339             2803            9993
     delete normal_to_sclc     ok         11768                207              462            6274
     delete normal_to_luad     ok         11768                130              287            5326
overexpress   sclc_to_luad     ok         14414               2650             3763           10637
overexpress sclc_to_normal     ok         14414                846             2131            5168
overexpress   luad_to_sclc     ok         15183               1537             5090            4697
overexpress luad_to_normal     ok         15183               3106             4758           10240
overexpress normal_to_sclc     ok         11754                464             1172            3470
overexpress normal_to_luad     ok         11754                720             1198            8423
```

## Coverage audit

```text
        arm     comparison                                                                                                                                                 path  exists status  rows missing_columns  duplicate_ensembl_ids  n_qualified  max_detections
     delete   sclc_to_luad             /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_sclc_to_luad.csv    True     ok 14415                                      0          405            2423
     delete sclc_to_normal           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_sclc_to_normal.csv    True     ok 14415                                      0          310            2423
     delete   luad_to_sclc             /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_luad_to_sclc.csv    True     ok 15185                                      0          955            6384
     delete luad_to_normal           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_luad_to_normal.csv    True     ok 15185                                      0         1339            6384
     delete normal_to_sclc           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_normal_to_sclc.csv    True     ok 11768                                      0          207             566
     delete normal_to_luad           /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/delete/heldout_allgene_delete_normal_to_luad.csv    True     ok 11768                                      0          130             566
overexpress   sclc_to_luad   /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_sclc_to_luad.csv    True     ok 14414                                      0         2650            2420
overexpress sclc_to_normal /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_sclc_to_normal.csv    True     ok 14414                                      0          846            2420
overexpress   luad_to_sclc   /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_luad_to_sclc.csv    True     ok 15183                                      0         1537            6384
overexpress luad_to_normal /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_luad_to_normal.csv    True     ok 15183                                      0         3106            6384
overexpress normal_to_sclc /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_normal_to_sclc.csv    True     ok 11754                                      0          464             566
overexpress normal_to_luad /home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/stats/overexpress/heldout_allgene_overexpress_normal_to_luad.csv    True     ok 11754                                      0          720             566
```

## Interpretation

A concordant hit requires significant, adequately detected, opposite-signed deletion and overexpression shifts. This is a model-level prioritization criterion, not evidence of experimental causality. Donor-level consistency, ambient-RNA/doublet sensitivity, and independent validation remain required before biological claims.

