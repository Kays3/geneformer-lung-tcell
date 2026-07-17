# Results and interpretation

## 1. Preliminary fine-tuning experiments

### LUAD/LUSC/normal minimum-1,000 experiment

The first perturbation-oriented model excluded COPD and oversampled only the
train/eval pool so every disease had at least 1,000 train/eval cells. Held-out
test donors were not oversampled.

| Metric | Result |
|---|---:|
| Test accuracy | 0.8381 |
| Test macro F1 | 0.8634 |

The natural test set contained 385 LUAD, 183 LUSC, and only 56 normal cells.
Consequently, these favorable metrics were treated as an implementation result,
not as sufficient evidence for the final perturbation model.

### LUAD/normal/COPD minimum-1,000 experiment

Donors were separated before oversampling, and the donor leakage audit passed.
Cells were then oversampled within each of train, eval, and test to at least
1,000 observations per disease.

| Metric | Result |
|---|---:|
| Oversampled test accuracy | 0.7747 |
| Oversampled test macro F1 | 0.7496 |

| True / predicted | Normal | COPD | LUAD |
|---|---:|---:|---:|
| Normal | 392 | 124 | 484 |
| COPD | 65 | 935 | 0 |
| LUAD | 3 | 0 | 997 |

Although donors remained isolated, duplicate test cells mean these are not
independent cell-level performance estimates. COPD was not included in the
final perturbation analysis.

## 2. Final no-oversampling classifier

### Overall performance

| Split | Accuracy | Macro F1 |
|---|---:|---:|
| Eval | 0.7772 | 0.7660 |
| Test | 0.7834 | 0.7577 |

### Held-out test confusion matrix

Rows are true disease labels and columns are predictions.

| True / predicted | LUAD | LUSC | Normal |
|---|---:|---:|---:|
| LUAD | 1,323 | 69 | 19 |
| LUSC | 249 | 311 | 0 |
| Normal | 390 | 5 | 1,013 |

### Per-class test performance

| Class | Precision | Recall | F1 | Support |
|---|---:|---:|---:|---:|
| LUAD | 0.6743 | 0.9376 | 0.7845 | 1,411 |
| LUSC | 0.8078 | 0.5554 | 0.6582 | 560 |
| Normal | 0.9816 | 0.7195 | 0.8303 | 1,408 |

### Interpretation before perturbation

- LUAD has high sensitivity but lower precision because 390 normal and 249
  LUSC cells were predicted as LUAD.
- Normal predictions are highly precise, but 28.1% of true normal cells were
  classified as LUAD.
- LUSC is the most difficult class by recall: 249 of 560 true LUSC cells were
  classified as LUAD.
- No LUSC test cells were classified as normal, and only five normal cells were
  classified as LUSC. The principal ambiguity is therefore LUAD versus the
  other states, especially LUAD versus LUSC.

This confusion pattern is important for perturbation interpretation. A shift
from LUSC toward LUAD may be easier for the model to produce than a shift from
LUSC directly toward normal. Candidate rankings should therefore consider both
target movement and movement away from the source.

## 3. Perturbation smoke test

One held-out LUSC cell from donor `Wu_Zhou_2021_P7` was used to verify the full
all-gene deletion pathway. The cell contained 1,101 eligible genes. Every gene
produced shifts against all three disease references, and the two expected raw
checkpoint files were written.

Selected single-cell effects are shown only as technical examples:

| Deleted gene | Shift toward LUAD | Shift toward normal | Shift toward LUSC |
|---|---:|---:|---:|
| RPL21 | +0.0357 | +0.0065 | -0.0414 |
| CBX3 | +0.0196 | +0.0035 | -0.0236 |
| TYROBP | +0.0189 | +0.0035 | -0.0210 |
| BTG1 | +0.0142 | +0.0050 | -0.0246 |
| RPS26 | -0.1094 | -0.0306 | +0.0965 |

RPL21 deletion moved this particular cell toward LUAD and away from LUSC,
whereas RPS26 deletion had the opposite effect. These are not biological hit
claims. The prominence of ribosomal genes in a single cell shows why final
ranking must use coverage, donor consistency, and sensitivity filters.

## 4. Full held-out deletion status

The production screen passed preparation, donor audit, state-reference
calculation, and smoke testing. It is processing all 3,379 held-out cells and
2,937,776 eligible cell-gene deletions in resumable 25-cell shards.

Current status, GPU utilization, memory, and hourly history are maintained in
[the live report](../heldout_allgene_perturbation/GPU_PROGRESS_REPORT.md).

No final perturbation ranking should be reported until all shards for a source
state have completed and the aggregate statistics have run.

## 5. Expected final outputs

For each of the six directional comparisons, the final tables are expected to
contain:

- Ensembl ID and gene name;
- mean shift toward the requested goal state;
- shift toward the alternate state;
- rank-sum significance against the empirical deletion distribution;
- FDR-adjusted significance;
- number of cells in which the gene was detected and deleted;
- downstream donor-level consistency summaries.

The six comparisons are LUAD to normal, LUAD to LUSC, LUSC to normal, LUSC to
LUAD, normal to LUAD, and normal to LUSC.

## 6. Main limitations

1. The held-out LUSC set contains only five donors and 560 cells.
2. Cell-level aggregation can over-weight donors contributing many cells;
   donor-level summaries are required.
3. Deleting a Geneformer token is an embedding perturbation, not a calibrated
   experimental knockout or partial knockdown.
4. Rank-abundant and ribosomal genes may create broad sequence-context effects.
5. The classifier's LUAD/LUSC confusion can influence apparent transition
   directions.
6. Statistical significance does not establish causal or therapeutic validity;
   robust hits require orthogonal biological validation.
