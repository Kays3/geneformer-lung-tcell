# Curated figures

`embeddings/` contains selected PNG outputs from the completed Stage 1 and
Stage 2 embedding analyses:

- joint UMAP views colored by cell type and disease;
- cell-type or disease centroid-similarity heatmaps;
- donor-effect UMAPs;
- Stage 2 donor-purity and donor-centroid diagnostics.

`classifiers/` contains the held-out confusion-matrix PDFs saved by the
Geneformer classifiers.

PNG files are deterministically downscaled to a maximum 1,800-pixel dimension
at 100 DPI for Git review; plotted values, labels, and composition are otherwise
unchanged. They document representation quality and donor structure; they are not
gene-perturbation results. Step 7 draft figures are intentionally excluded
because that notebook contains placeholder metrics and random fallback data.
