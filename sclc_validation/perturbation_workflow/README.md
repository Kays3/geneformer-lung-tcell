# SCLC / LUAD / normal T-cell fine-tuning and perturbation

**Compute environment:** NVIDIA GB10 (DGX Spark), Geneformer V2 104M
**Disease states:** small cell lung carcinoma (SCLC), lung adenocarcinoma (LUAD), normal

This is the SCLC-inclusive counterpart to
[`current_workflow/`](../../current_workflow/README.md)'s LUAD/LUSC/normal T-cell
classifier and all-gene perturbation screen. It follows the same donor-aware
design (donor-disjoint splitting, training-cells-only reference centroids,
shard-checkpointed execution) but on the HTAN/CELLxGENE SCLC cohort identified
by [`sclc_validation/audit/`](../audit/README.md), and adds an in-silico
overexpression arm alongside deletion.

Run status and final numbers are recorded in `METHODS.md` and `RESULTS.md`
(the latter written once the perturbation screen completes).

## Data source and cohort

HTAN/CELLxGENE T-cell object only (46,140 T cells, 42 donors: 19 SCLC, 22 LUAD,
4 normal) -- no other dataset is blended in. See `METHODS.md` for the full
cohort characterization, including the donor-disjoint-across-disease guard for
three patients (RU675, RU682, RU684) who each contribute both LUAD tumor and
normal-tissue T cells.

## Repository map

```text
METHODS.md    cohort, split design, fine-tuning and perturbation methodology
RESULTS.md    classifier metrics and perturbation findings (written after the run)
```

Large tokenized datasets, embeddings, checkpoints, and raw perturbation
pickles remain outside Git, under `KD/sclc_luad_normal_htan_finetune/` and
`KD/sclc_luad_normal_htan_heldout_allgene_perturbation/` on the compute host.
