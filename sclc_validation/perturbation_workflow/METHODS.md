# Methods and reproducibility

## 1. Source cohort

HTAN/CELLxGENE T-cell object, the primary cohort identified by the SCLC data
feasibility audit (`sclc_validation/audit/SCLC_DATA_FEASIBILITY_AUDIT.md`):

```text
source_metadata/HTAN_MSK_T_cells_6fde3ad9.h5ad
(CELLxGENE dataset 6fde3ad9-c2dc-4bea-bcb1-100192dd5877, HTAN MSK collection)
```

46,140 T cells, 24,540 genes, raw integer counts in `.raw.X` (Ensembl gene
IDs), 42 donors across three disease states:

| Disease | Cells | Donors |
|---|---:|---:|
| Small cell lung carcinoma (SCLC) | 11,791 | 19 |
| Lung adenocarcinoma (LUAD) | 29,829 | 22 |
| Normal | 4,520 | 4 |

No other dataset is blended in. An earlier plan to hybridize the thin
normal class with the LuCA NSCLC atlas was explicitly dropped in favor of
using HTAN alone, accepting the normal class's low donor count as a named
limitation (see "Limitations" below).

### Donor identity resolution

- `PleuralEffusion` is a `donor_id` string used for one SCLC biospecimen
  (1,816 cells) instead of an "RU####" clinical code. `HTAN_Participant_ID`
  confirms it maps to exactly one participant (`HTA8_2001`) -- a single
  genuine donor, not multiple identities merged under one label.
- Three `donor_id` values (`RU675`, `RU682`, `RU684`) each map to **two**
  `HTAN_Participant_ID` values -- one for a tumor (LUAD) biospecimen, one for
  a normal-tissue biospecimen, from the same clinical patient. This confirms
  these are genuine paired tumor/normal donors, not a `donor_id` collision
  artifact. `donor_id` (which correctly groups paired biospecimens from the
  same patient) is therefore used as the donor-disjoint splitting key, not
  the finer-grained `HTAN_Participant_ID`.

## 2. Donor-disjoint splitting

Cells are assigned to `train`, `eval`, and `test` at the donor level via a
greedy, cell-count-balanced assignment (target ~60/20/20 by cells within each
disease, seed 43).

**Cross-disease donor guard**: `RU675`, `RU682`, and `RU684` contribute cells
to both the LUAD and normal disease rows. Each is pinned to a single split
across *both* rows before the remaining single-disease donors are assigned,
so a patient's tumor cells and that same patient's normal cells cannot land
in different splits. This extends the standard donor-leakage check (`assert
no donor appears in more than one split`) to also assert no donor appears
under more than one split *across* disease labels.

### Final split

| Disease | Split | Cells | Donors |
|---|---|---:|---:|
| LUAD | train | 17,831 | 14 |
| LUAD | eval | 5,611 | 4 |
| LUAD | test | 6,387 | 4 |
| Normal | train | 2,334 | 2 |
| Normal | eval | 1,620 | 1 |
| Normal | test | 566 | 1 |
| SCLC | train | 7,037 | 13 |
| SCLC | eval | 2,330 | 3 |
| SCLC | test | 2,424 | 3 |

Both the donor-leakage check and the cross-disease guard passed with zero
leaked donors.

## 3. Tokenization

Raw counts (`.raw.X`) tokenized with Geneformer's `TranscriptomeTokenizer`
(V2, gc104M gene median / token dictionaries), retaining `cell_id`,
`individual`, `celltype`, `disease`, `split`, `length`.

```text
KD/sclc_luad_normal_htan_finetune/data/sclc_luad_normal_htan_tcells.dataset
```

## 4. Fine-tuning

| Parameter | Value |
|---|---|
| Base model | Geneformer-V2-104M |
| Model type | Cell classifier |
| Prediction label | `disease` |
| Classes | SCLC, LUAD, normal |
| Epochs | 1 |
| Learning rate | `5e-5` |
| Training batch size | 8 |
| Forward/evaluation batch size | 16 |
| Frozen transformer layers | 6 |
| Random seed | 43 |
| Oversampling | None |

Fitted on train donors, validated during development on eval donors, tested
on the untouched donor-held-out test set.

```text
KD/sclc_luad_normal_htan_finetune/scripts/run_finetune.py
```

## 5. Classification and metrics report

Run on the donor-held-out test set only. Reports overall accuracy and macro
F1, per-class precision/recall/F1 with both cell-count and donor-count
support, a confusion matrix (table + chart), and a context comparison
against the prior LUAD/LUSC/normal classifier (0.7834 accuracy, 0.7577 macro
F1 -- different classes and cohort, calibration reference only). Normal-class
metrics are flagged explicitly given 1 test donor.

This report does not gate the pipeline -- the run proceeds directly into the
perturbation screen regardless of classifier quality, per an explicit
"run end-to-end unattended, no pause" instruction. Interpret perturbation
findings in light of the classifier metrics.

```text
KD/sclc_luad_normal_htan_finetune/scripts/report_classification_metrics.py
```

### One test cell is not scored, by upstream design

The recorded confusion matrix totals **9,376** cells while the tokenized test
split holds **9,377**. The missing cell is real and reproducible, and it is not
a fault in this pipeline.

Geneformer's `classifier_predict` (`geneformer/evaluation_utils.py`) truncates
the evaluation set when the final batch would contain exactly one example:

```python
# ensure there is at least 2 examples in each batch to avoid incorrect tensor dims
evalset_len = len(evalset)
max_divisible = find_largest_div(evalset_len, forward_batch_size)
if len(evalset) - max_divisible == 1:
    evalset_len = max_divisible
```

The guard exists because the loop calls `torch.squeeze(outputs.logits)`, which
on a one-example batch collapses the batch dimension from `(1, 3)` to `(3,)` and
corrupts the subsequent `torch.cat`. Upstream avoids the case rather than
handling it.

With this cohort it fires because `9377 % 16 == 1` at the configured
`forward_batch_size=16` — and it would fire at 8 or 32 as well, since 9,377 is
one more than a multiple of all three. The dropped cell is the final row of
`sclc_luad_normal_htan_labeled_test.dataset`:

| | |
|---|---|
| `cell_id` | `RU1138_230816753564460` |
| donor | RU1138 |
| celltype | CD4-positive helper T cell |
| label | lung adenocarcinoma |

which matches the per-class arithmetic exactly: LUAD 6,386 recorded against
6,387 in the data, with SCLC and normal unaffected.

**Impact: none at any reported precision.** Scoring all 9,377 cells gives
accuracy 0.919377 and macro F1 0.903318, against the recorded 0.919369 and
0.903315. The headline 91.9% / 0.903 is unchanged.

**What to do with this.** The recorded confusion matrix is a valid 9,376-cell
result, not a 9,377-cell one. Anyone recomputing metrics directly from the
tokenized dataset will land one cell off and should not read that as an error.
`tools/verify_gpu_env.py` reproduces the full-dataset numbers and reports the
coverage difference explicitly rather than failing on it.

## 6. Held-out all-gene deletion + overexpression screen

### Reference states

Mean disease-state CLS embeddings calculated from **training cells only**.
Held-out test cells never contribute to the SCLC, LUAD, or normal reference
centroids.

### Perturbation unit and types

For each held-out cell, every non-special gene token present in that cell's
full-context ranked Geneformer sequence is perturbed individually, once per
perturbation type:

- **Delete**: gene removed from the rank-value encoding (as in the prior
  LUAD/LUSC/normal workflow).
- **Overexpress**: gene moved to the front of the rank-value encoding
  (Geneformer's native `perturb_type="overexpress"`), simulating maximal
  relative expression.

Both types use the full gene universe present in each cell -- no
highly-variable-gene restriction. (An earlier plan to restrict to the top
3000 variable genes was reverted: Geneformer's `InSilicoPerturber` has no
built-in "sweep individually, restricted to a list" mode -- passing a gene
list perturbs the whole list as one combined event, not one gene at a time.
A true restricted sweep would have required pre-filtering each cell's
tokenized context down to the target genes before perturbing, changing the
"original"/"perturbed" embeddings from full-transcriptome context to
reduced context. Running all genes avoids that deviation and matches the
prior workflow's methodology.)

For reference state \(s\) and perturbation type \(t\):

```text
shift_s = cosine(perturbed_cell, reference_s) - cosine(original_cell, reference_s)
```

Positive = movement toward state `s`; negative = movement away. A gene whose
**deletion** moves cells away from a state and whose **overexpression**
moves them toward it (or vice versa) is a substantially stronger candidate
than either signal alone -- this concordance is the primary evidence tier in
`RESULTS.md`.

### Efficient three-source design (per perturbation type)

Each cell-gene perturbation is computed once per type while its shift is
scored against all three references. Three source-state screens recover six
directional comparisons, per type (12 directional comparisons total across
delete + overexpress):

| Source screen | Directional comparisons recovered |
|---|---|
| SCLC | SCLC to LUAD; SCLC to normal |
| LUAD | LUAD to SCLC; LUAD to normal |
| Normal | normal to SCLC; normal to LUAD |

### Execution and recovery

- Perturbation types: delete, overexpress (run and stored separately --
  `InSilicoPerturberStats` ingests every `*_raw.pickle` file in a directory
  regardless of type, so mixing them would silently corrupt the stats).
- Shard size: 25 cells. Forward batch size: 16. Data workers: 4.
- Embedding mode: V2 CLS, layer offset 0.
- Statistics mode: `goal_state_shift` with FDR correction.
- A shard receives a completion marker only after all of its cells succeed
  for that perturbation type; an interrupted run resumes without repeating
  completed shards.

```text
KD/sclc_luad_normal_htan_heldout_allgene_perturbation/scripts/run_heldout_allgene.py
```

## 7. Execution

Runs as a single unattended background job on the remote GB10 box: cohort
prep -> tokenize/fine-tune -> classification & metrics report -> perturbation
prepare -> state-embeddings -> smoke-test (both types) -> full perturb sweep
(both types, all sources) -> stats. No pause between stages.

```text
KD/sclc_luad_normal_htan_finetune/scripts/run_pipeline.sh
```

## 8. Planned aggregation

Same aggregation discipline as the prior workflow, extended for the second
perturbation type:

1. minimum perturbation coverage per gene, per type;
2. effect size toward the goal state and away from the source state, per type;
3. false-discovery-rate correction, per type;
4. **deletion vs. overexpression concordance per gene** -- the primary
   evidence tier;
5. consistency across held-out donors rather than cell-count weighting alone;
6. sensitivity analysis excluding ribosomal and other rank-dominant genes;
7. pathway-level interpretation of robust genes;
8. targeted reruns of top candidates if further validation is needed.

Results are scoped first to the audit's pre-registered 21-gene panel
(`PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX, LAYN, NKG7, GNLY, PRF1, GZMB, GZMH,
IFNG, TCF7, SLAMF6, IL7R, CCR7, ASCL1, NEUROD1, POU2F3, YAP1`), with
broader all-gene findings presented as exploratory.

## Limitations

- **Thin normal class**: only 4 donors (2 train / 1 eval / 1 test). Normal
  eval/test metrics and perturbation reference embeddings rest on 1 donor
  each -- they cannot distinguish a real population-level signal from that
  single donor's idiosyncrasies. This is a named limitation, not silently
  treated as equivalent in statistical power to the SCLC/LUAD arms.
- **No pipeline checkpoint**: the run proceeds from fine-tuning straight into
  the (expensive) perturbation screen regardless of classifier quality, per
  explicit instruction. If classifier metrics turn out to be poor, the
  perturbation results should be treated as exploratory pending a
  re-fine-tune.
