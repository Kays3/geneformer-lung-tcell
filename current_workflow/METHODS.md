# Methods and reproducibility

## 1. Source atlas and cell selection

The final cohort was sampled from the original integrated NSCLC atlas:

```text
/home/thinkstation2/workspace/KD/data/nsclc/nsclc_integrated.h5ad
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
cd /home/thinkstation2/workspace

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
unified memory. Monitoring records GPU utilization, temperature, power,
CUDA-process memory, system-memory availability, cells written, completed
shards, and raw output counts. The live report is generated by
`current_workflow/monitoring/report_generation_job.sh` on a 15-minute cadence
until the run finishes.

## 8. Planned aggregation

Final reporting should require:

1. minimum deletion coverage per gene;
2. effect size toward the goal state and away from the source state;
3. false-discovery-rate correction;
4. consistency across held-out donors rather than cell-count weighting alone;
5. sensitivity analysis excluding ribosomal and other rank-dominant genes;
6. pathway-level interpretation of robust genes;
7. targeted reruns of top candidates if further validation is needed.

## Class-construction defects (established 2026-09-24)

Three defects in how the LUAD/LUSC/normal classes in §1 were built, found
while auditing a different classifier for a study confound and then checked
back against this one. Each is established by direct verification against
the actual data this experiment used, not by inference from the selection
rule alone. This is a record of what was found and how; it is not a claim
about what it does to the numbers in §5–§8 (see "What this does not cover"
below).

### 1. Study confound (normal vs. tumor)

No study in the source atlas contributes both a normal-labelled and a
tumor-labelled (LUAD or LUSC) T cell. Verified directly on the whole source
atlas (`disease` in `normal`/`lung adenocarcinoma`/`squamous cell lung
carcinoma`, no cell-type filter): 265 donors across 18 studies, disease label
is constant within every donor (0/265 donors carry more than one disease
label) and study is constant within every donor (0/265 donors span more than
one study), and of the 18 studies, 0 contribute cells to both a normal label
and a tumor label. The `normal` class and the tumor classes are therefore
drawn from entirely disjoint sets of studies. The `LUAD`-vs-`LUSC` boundary
is **not** subject to this defect: every one of LUSC's contributing studies
is also a LUAD study.

### 2. Class definition mixes tissue origin, not just disease label

§1's selection criteria filter on `disease`, `cell_type_major`, donor/count
validity, and singlet status — there is no `origin` term. Confirmed against
the actual 21,000 cells that trained the classifier, not the selection rule
in the abstract: the pre-tokenization subset (`KD/tcell_luad_lusc_normal_10k_from_atlas/outputs/balanced_lusc_max_7000_per_disease_tcells.h5ad`)
and the tokenized training dataset (§3) were joined on `cell_id`; for both
the LUAD and LUSC classes, all 7,000 training `cell_id` values matched
exactly, 7000/7000, zero mismatches — an exact-identity join, not a
resemblance check.

**LUAD (7,000 cells), by `origin`:**

| origin | cells | % |
|---|---:|---:|
| tumor_primary | 2,950 | 42.1% |
| tumor_metastasis | 1,920 | 27.4% |
| normal_adjacent | 1,190 | 17.0% |
| normal | 940 | 13.4% |

58% of the LUAD class is not primary tumor tissue.

**LUSC (7,000 cells), by `origin`:**

| origin | cells | % |
|---|---:|---:|
| tumor_primary | 4,990 | 71.3% |
| normal_adjacent | 2,010 | 28.7% |

28.7% of the LUSC class is not primary tumor tissue — a narrower mix than
LUAD (two origins, not four), but not clean.

**Normal (7,000 cells), by `origin`:**

| origin | cells | % |
|---|---:|---:|
| normal | 6,638 | 94.8% |
| normal_adjacent | 362 | 5.2% |

The normal class is not 100% `origin == normal` either.

**Bounded range, not a point estimate.** The atlas carries a second origin
field, `origin_fine`, which mostly agrees with `origin` but not always: for
LUAD, `origin_fine` reclassifies 243 of the `normal_adjacent`-per-`origin`
cells as plain `normal`, and for the normal class it reclassifies 345 cells
the same way. LUSC shows no such disagreement. The two fields are not
reconciled here — which one is correct is not established, only that they
disagree. Reported as a range: the LUAD class's `normal`-origin content is
**13.4% under `origin`, 16.9% under `origin_fine`**; the normal class's
`origin == normal` purity is **94.8% under `origin`, 99.8% under
`origin_fine`**.

### 3. Asymmetric pollution between the tumor classes

`tumor_metastasis` and `normal` origin appear in the LUAD class (27.4% and
13.4% respectively, §2 above) and are entirely absent from the LUSC class
(0.0% both). Within the LUAD-vs-LUSC contrast specifically, this means
metastatic tissue origin is a perfect predictor of the LUAD label, covering
27.4% of that class, independent of any adenocarcinoma-versus-squamous
biology.

### How defects 1 and 2 interlock

Disease label is donor-constant (§1), so the 940 `origin == normal` cells
inside the LUAD class (§2) are LUAD **patients'** normal-origin tissue, not
healthy donors' tissue. The normal class is 94.8–99.8% `origin == normal`
(§2). So for that ~13–17% slice of the LUAD class, the contrast the
classifier is actually drawing on is: normal-origin lung tissue from a
lung-cancer patient, versus normal-origin lung tissue from a person without
one — and because disease status and study identity are perfectly
confounded (§1), the only thing observable in the data that predicts which
side of that contrast a cell falls on is which study it came from.

### What this does not cover

- **The SCLC/LUAD/normal classifier** (`sclc_validation/perturbation_workflow/`)
  is built from a different atlas (a single-site HTAN/CELLxGENE object, HTA8)
  and is not covered by any finding above. It carries its own separately
  documented tissue-of-origin skew.
- **Downstream impact is not quantified.** Nothing above measures how much,
  if at all, these defects changed the fine-tuning result in §4–§5 or the
  perturbation results in §6–§8. This section records defects in class
  construction; it does not establish that the published numbers are wrong.
- **The LUAD-vs-LUSC boundary is not study-confounded** (§1) — every LUSC
  study is also a LUAD study. That is true and separate from defect 3: a
  boundary can be clean on study identity and still carry a perfect
  tissue-origin shortcut.

Verification for all of the above: obs-level metadata only, read via
`h5py`/`anndata` in backed mode against the source atlas and the two
artifacts named in §2 and §3; expression data (`X`) was never loaded for
this audit. All figures independently reproducible from files already
referenced in this document.
