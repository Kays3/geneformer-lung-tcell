# Methods and reproducibility

## 1. Source atlas and cell selection

The final cohort was sampled from the original integrated NSCLC atlas:

```text
/home/petadimensionlab/workspace/Geneformer/KD/data/nsclc/nsclc_integrated.h5ad
```

Selection criteria were:

- disease in `normal`, `lung adenocarcinoma`, or
  `squamous cell lung carcinoma`;
- `cell_type_major` in `T cell CD4` or `T cell CD8`;
- singlet cells;
- non-missing donor identifier;
- valid count information;
- no COPD in the final LUAD/LUSC/normal cohort;
- no cell oversampling.

The feasibility audit showed that the original target of 10,000 clean cells per
disease per split was impossible without oversampling because strict LUSC T-cell
availability was limiting. The final design therefore maximized clean LUSC use
at approximately 7,000 cells and balanced LUAD and normal to the same total.

### Final selected cohort

| Disease | CD4 T cells | CD8 T cells | Total |
|---|---:|---:|---:|
| LUAD | 3,882 | 3,118 | 7,000 |
| LUSC | 2,879 | 4,121 | 7,000 |
| Normal | 3,395 | 3,605 | 7,000 |
| **Total** | **10,156** | **10,844** | **21,000** |

The resulting matrix contained 21,000 cells and 17,764 genes.

## 2. Donor-disjoint splitting

Cells were assigned to `train`, `eval`, and `test` at the donor level. Every
donor occurs in exactly one split. The leakage audit passed before fine-tuning
and was repeated before perturbation.

| Split | LUAD cells/donors | LUSC cells/donors | Normal cells/donors | Total cells |
|---|---:|---:|---:|---:|
| Train | 4,503 / 58 | 5,458 / 17 | 4,321 / 37 | 14,282 |
| Eval | 1,086 / 15 | 982 / 4 | 1,271 / 9 | 3,339 |
| Test | 1,411 / 19 | 560 / 5 | 1,408 / 12 | 3,379 |

Cell counts are not exactly balanced within each split because donor isolation
was prioritized over cell-level balance. This is especially visible for LUSC,
where a small number of donors contribute many cells.

## 3. Tokenization

A slim h5ad copy was created for Geneformer. Raw integer counts came from
`layers['count']`; normalized expression was not substituted for tokenization.
The tokenized dataset retained these metadata fields:

- `cell_id`
- `individual`
- `celltype`
- `disease`
- `split`
- `length`

Local tokenized dataset:

```text
KD/tcell_luad_lusc_normal_luscmax7000_finetune/data/
balanced_lusc_max_7000_per_disease_tcells.dataset
```

## 4. Fine-tuning

| Parameter | Value |
|---|---|
| Base model | Geneformer-V2-104M |
| Model type | Cell classifier |
| Prediction label | `disease` |
| Classes | LUAD, LUSC, normal |
| Epochs | 1 |
| Learning rate | `5e-5` |
| Training batch size | 8 |
| Forward/evaluation batch size | 16 |
| Frozen transformer layers | 6 |
| Random seed | 43 |
| Oversampling | None |

The model was fitted only on training donors, evaluated during development on
eval donors, and tested on the untouched donor-held-out test dataset.

Local runner:

```text
KD/tcell_luad_lusc_normal_luscmax7000_finetune/scripts/run_finetune.py
```

Local selected model:

```text
KD/tcell_luad_lusc_normal_luscmax7000_finetune/runs/
260717_geneformer_cellClassifier_tcell_luad_lusc_normal_luscmax7000/ksplit1
```

## 5. Pre-perturbation evaluation

The saved model was evaluated on all 3,379 test cells. No test oversampling or
duplicate-cell balancing was applied. Accuracy, macro F1, per-class metrics,
and a confusion matrix were saved before any perturbation analysis.

## 6. Held-out all-gene deletion

### Reference states

Exact mean disease-state CLS embeddings were calculated from **training cells
only**. The held-out test cells were not used to construct the LUAD, LUSC, or
normal reference centroids.

### Perturbation unit

For each held-out cell, every non-special gene token present in that cell's
ranked Geneformer sequence is deleted once. The shortened sequence is passed
through the fine-tuned model and its CLS embedding is compared with all three
disease references.

For reference state \(s\), the recorded effect is:

```text
shift_s = cosine(perturbed_cell, reference_s)
          - cosine(original_cell, reference_s)
```

A positive value indicates movement toward state `s`; a negative value
indicates movement away from it. Because deletion removes a ranked token and
changes sequence context, this represents a complete in silico deletion, not a
quantitative partial knockdown.

### Efficient three-source design

Each cell-gene deletion is computed once while its shift is scored against all
three references. Three source-state screens therefore recover six directions:

| Source screen | Directional comparisons recovered |
|---|---|
| LUAD | LUAD to normal; LUAD to LUSC |
| LUSC | LUSC to normal; LUSC to LUAD |
| Normal | normal to LUAD; normal to LUSC |

This avoids repeating the same forward pass in six independent pairwise runs.

### Exact held-out workload

| Source | Test cells | Test donors | Mean tokens | Valid deletions |
|---|---:|---:|---:|---:|
| LUAD | 1,411 | 19 | 817.1 | 1,150,097 |
| LUSC | 560 | 5 | 624.0 | 348,313 |
| Normal | 1,408 | 12 | 1,024.3 | 1,439,366 |
| **Total** | **3,379** | **36** | — | **2,937,776** |

CLS and EOS special tokens are excluded from the deletion counts.

### Execution and recovery

- Source order: LUSC, LUAD, then normal.
- Shard size: 25 cells.
- Total shards: 137.
- Forward batch size: 16.
- Data workers: 4.
- Perturbation type: deletion.
- Genes: all eligible genes present in each cell.
- Embedding mode: V2 CLS, layer offset 0.
- Statistics mode: `goal_state_shift` with FDR correction.

A shard receives a completion marker only after all of its cells succeed. An
interrupted run can therefore restart without repeating completed shards.

Local commands:

```bash
cd /home/petadimensionlab/workspace/Geneformer

.venv/bin/python \
  KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/scripts/run_heldout_allgene.py \
  prepare

.venv/bin/python \
  KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/scripts/run_heldout_allgene.py \
  state-embeddings

.venv/bin/python \
  KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/scripts/run_heldout_allgene.py \
  smoke-test

.venv/bin/python \
  KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/scripts/run_heldout_allgene.py \
  all
```

## 7. Compute monitoring

The full run executes on an NVIDIA GB10 with approximately 119 GiB visible
unified memory. Hourly monitoring records GPU utilization, temperature, power,
CUDA-process memory, system-memory availability, cells written, completed
shards, and raw output counts. The report is committed automatically to
`hokudai_spark1` until the run finishes.

## 8. Planned aggregation

Final reporting should require:

1. minimum deletion coverage per gene;
2. effect size toward the goal state and away from the source state;
3. false-discovery-rate correction;
4. consistency across held-out donors rather than cell-count weighting alone;
5. sensitivity analysis excluding ribosomal and other rank-dominant genes;
6. pathway-level interpretation of robust genes;
7. targeted reruns of top candidates if further validation is needed.
