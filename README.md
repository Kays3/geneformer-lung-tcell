# Geneformer NSCLC workflow

This branch contains the reproducible notebook workflow developed in the local
`Geneformer/KD` workspace for non-small-cell lung cancer (NSCLC) analysis.
Notebook outputs are cleared before version control so local paths, generated
figures, and large embedded results are not stored in Git.

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

## Local assets

The notebooks expect generated data beneath the repository root and a local
Geneformer base model in a sibling `Geneformer-V2-104M` directory. Large data,
tokenized datasets, training runs, embeddings, model weights, and figures are
excluded by `.gitignore` and must be generated or provisioned locally.

Create an environment with Python 3.12 and install the packages in
`requirements.txt`. Geneformer itself must also be installed or otherwise
available in that environment.
