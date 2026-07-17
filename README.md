# Geneformer NSCLC workflow

This branch contains the reproducible notebook workflow developed in the local
`Geneformer/KD` workspace for non-small-cell lung cancer (NSCLC) analysis.
Selected low-resolution notebook outputs and curated evidence are included for
review. Large datasets, checkpoints, model weights, and full generated runs
remain excluded from Git.

## Workflow

Run the notebooks from the repository root in numerical order:

1. `Step1-non-small-cell-download-explore-cell-type.ipynb` — download and explore the NSCLC data.
2. `Step2-Geneformer-27k-cell-prep.ipynb` — prepare the 27,000-cell stage-one dataset.
3. `Ste3-Geneformer-embeddign-27k.ipynb` — tokenize and train the primary cell-type model.
4. `Step4-Geneformer-27k.ipynb` — prepare and train the disease model.
5. `Step5-Stage1-embeddings.ipynb` — extract and analyze stage-one embeddings.
6. `Step6-Stage2-embeddings.ipynb` — extract and analyze stage-two embeddings.
7. `Step7-synthesis-stage1-stage2.ipynb` — compare and synthesize both stages.

The step-three filename is retained from the working directory to preserve its
provenance.

## Results and evidence

- `tables/` contains compact source and derived CSVs, including the T-cell
  cancer-versus-normal cohort, classifier metrics, and donor-signal summaries.
- `figures/embeddings/` contains selected embedding and donor-effect figures.
- `figures/classifiers/` contains the two held-out confusion matrices.
- `reports/geneformer_nsclc_progress_summary.html` is a self-contained
  technical report describing verified progress, limitations, and the proposed
  T-cell perturbation study.
- `reports/artifact.json` is the canonical source used to build the HTML.

Step 7 remains a draft synthesis notebook. Its hard-coded placeholder metrics
and random fallback matrices are not treated as verified evidence in the report.

## Local assets

The notebooks expect generated data beneath the repository root and a local
Geneformer base model in a sibling `Geneformer-V2-104M` directory. Large data,
tokenized datasets, training runs, embeddings, and model weights are excluded by
`.gitignore` and must be generated or provisioned locally.

Create an environment with Python 3.12 and install the packages in
`requirements.txt`. Geneformer itself must also be installed or otherwise
available in that environment.
