# Prior NSCLC workflow archive

This directory preserves the earlier whole-cohort Geneformer workflow. It is
kept for provenance and comparison; the active T-cell fine-tuning and
perturbation experiment is in [`current_workflow/`](../../current_workflow/).

## Contents

```text
notebooks/   Step1-Step7 analysis notebooks
figures/     classifier and embedding figures
tables/      curated metrics and donor diagnostics
geneformer_nsclc_progress_summary.html
artifact.json
source-notes.md
```

## Notebook order

1. `Step1-non-small-cell-download-explore-cell-type.ipynb`
2. `Step2-Geneformer-27k-cell-prep.ipynb`
3. `Ste3-Geneformer-embeddign-27k.ipynb`
4. `Step4-Geneformer-27k.ipynb`
5. `Step5-Stage1-embeddings.ipynb`
6. `Step6-Stage2-embeddings.ipynb`
7. `Step7-synthesis-stage1-stage2.ipynb`

The original Step3 filename is preserved. Step7 contains draft synthesis logic;
placeholder metrics and random fallback matrices are not treated as evidence.

Earlier held-out results were 0.9221 accuracy/0.9211 macro F1 for the nine-class
cell-type model and 0.7360 accuracy/0.7423 macro F1 for the four-class disease
model. These tasks differ from today's three-state T-cell classifier.
