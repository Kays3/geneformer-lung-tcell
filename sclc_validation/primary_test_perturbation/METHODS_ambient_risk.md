# Methods: ambient-RNA risk diagnostic

## 1. Why this is not CellBender

The planned sensitivity analysis was CellBender `remove-background`. **It cannot be run on
this dataset.** CellBender models the ambient pool by contrasting cell-containing droplets
against **empty droplets**, and estimates a per-cell contamination fraction from that
contrast. Both requirements fail here:

| Requirement | Status in this dataset |
|---|---|
| Unfiltered droplet matrix including empty barcodes | Absent. Minimum **113 UMI/cell** across all 46,140 cells; nothing in the empty-droplet range |
| All cell types from the digest | Absent. Already subset to 6 T-cell subtypes, all `filter_pass == 1` |
| Raw CellRanger output on disk | Absent on both nodes; no `raw_feature_bc_matrix` anywhere |

The CELLxGENE collection that supplied this cohort distributes only cell-type-partitioned,
filtered matrices (`raw_data_location: raw.X` denotes raw *counts* inside a filtered
matrix, not raw droplets). With no empty droplets, the ambient profile is not
identifiable, and fitting the model anyway would produce a contamination estimate with
nothing constraining it.

Running CellBender therefore requires obtaining primary CellRanger output from the HTAN
data portal — a data-access task, not a compute task, and plausibly controlled-access.
That is the only route to the originally requested analysis and remains open.

## 2. What is computed instead

The question behind the request is: *are the prioritized candidates, especially the
antigen-presentation program, artifacts of ambient contamination?* That question can be
attacked without a correction model, by asking whether each gene **behaves** like a
contaminant.

Four properties separate ambient signal from genuine expression:

| Property | Ambient behaviour | Feature used |
|---|---|---|
| Lineage foreignness | Present in cells that cannot transcribe it | anchor-set membership |
| Cell-state structure | None — the cell did not make it | `log_subtype_f`, one-way F across the 6 T-cell subtypes |
| Sample dependence | High — soup composition differs per digest | `log_donor_f`, one-way F across donors |
| Droplet-volume dependence | High — contamination scales with content | `libsize_corr`, expression vs. total UMI on ranks |
| Breadth without depth | Detected widely, never strongly | `detect_frac / mean_expr` |

`subtype_over_donor` is the ratio of the two F statistics: genuine genes are organized by
cell state, ambient genes by sample.

## 3. Calibration against known answers

The score is only meaningful if it recovers genes whose status is already known, so it is
fit rather than asserted. Two anchor sets are defined **a priori** in the script:

- **Ambient anchors (39 defined, 34 retained):** surfactant (`SFTPC`, `SFTPB`, `SCGB1A1`),
  haemoglobin (`HBB`, `HBA1/2`, `ALAS2`), myeloid (`LYZ`, `S100A8/9`, `CD68`, `MARCO`),
  epithelial (`EPCAM`, `KRT8/18/19`), stromal and endothelial (`COL1A1`, `DCN`, `PECAM1`,
  `VWF`). A T cell does not transcribe these, so detection is contamination by construction.
- **T-cell anchors (39 defined, 36 retained):** canonical T-cell-intrinsic transcripts
  (`CD3D/E/G`, `TRAC`, `TRBC1/2`, `LCK`, `ZAP70`, `IL7R`, `TCF7`, `GZMA`, `FOXP3`).

A logistic regression on the six features is fit to the anchors with balanced class
weights, and **5-fold cross-validated ROC AUC is reported before any candidate is
ranked**. The script warns if AUC < 0.8 and the score should then be discarded.

Genes are scored only if detected in ≥0.5% of cells (11,047 of 24,540 genes); below that
the features are dominated by sampling noise.

## 4. Flag threshold

A candidate is flagged only if it scores at least as ambient-like as a genuine
contaminant: the threshold is the **25th percentile of the known-ambient anchors**
(0.950). This is deliberately conservative — it will not flag a gene merely for being
above average, only for sitting inside the contaminant distribution.

## 5. Identifier handling

The HTAN h5ad is indexed by Ensembl ID and carries no gene-symbol column. Symbols are
mapped from the perturbation stats tables, which carry both identifiers, and candidates
are joined to features on **Ensembl ID, not symbol**, because symbols are not one-to-one.
An earlier run silently matched zero anchors for exactly this reason; the script now fails
loudly if the mapping cannot be built.

## 6. Reproducibility

### Environment

Executed on `thinkstation2` (NVIDIA GB10, `Linux-6.17.0-1029-nvidia-aarch64`, glibc 2.39):

```text
python 3.12.13   numpy 2.4.4    scipy 1.18.0   pandas 3.0.5
sklearn 1.9.0    anndata 0.13.2 scanpy 1.12.3
```

No GPU is used; this analysis is CPU-only and takes ~1 minute.

### Command

Run from the repository checkout, not a scratch copy, so the committed code is the code
that ran (see [`../../WORKFLOW.md`](../../WORKFLOW.md)):

```bash
./tools/sync.sh
ssh ts2 'cd ~/workspace/geneformer-lung-tcell && \
  ~/workspace/geneformer-uv-starter/.venv/bin/python \
  sclc_validation/primary_test_perturbation/scripts/ambient_risk_diagnostic.py'
```

### Inputs, all overridable by environment variable

| Variable | Default | Purpose |
|---|---|---|
| `HTAN_H5AD` | `~/workspace/KD/.../htan_sclc_luad_normal_tcells_prepared.h5ad` | 46,140 × 24,540 count matrix |
| `SCLC_PERTURBATION_ROOT` | `~/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation` | stats tables, for the symbol mapping |
| `DENOISED_CANDIDATES` | `tables/immune_cancer_candidates_with_donor_robustness.csv` | candidate list to score |
| `AMBIENT_OUT_DIR` | `tables/` | output location |

### Determinism

`LogisticRegression` with `max_iter=2000` on standardized features is deterministic, and
`cross_val_score` uses non-shuffled 5-fold splits, so no seed is required — repeated runs
on the same inputs give identical output. The anchor sets are literals in the script, not
derived from the data, so they cannot drift with the candidate list.

### Outputs

| File | Contents |
|---|---|
| `tables/ambient_risk_candidates.csv` | Per-candidate risk, percentile, flag |
| `tables/ambient_risk_by_program.csv` | Median risk and flag count per program |
| `tables/ambient_risk_all_genes.csv` | All 11,047 scored genes with features |
| `tables/ambient_risk_manifest.json` | Anchors used, AUC, threshold, cell/gene counts |
| `figures/ambient_risk/ambient_risk_distribution.png` | Candidates against both anchor distributions |

## 7. Limitations

- **This is a diagnostic, not a correction.** No counts are modified. It tells you which
  candidates are at risk, not what their de-contaminated values would be. Nothing
  downstream has been re-run on corrected data, because no corrected data exists.
- **The near-perfect anchor separation reflects deliberately extreme anchors.** Surfactant
  versus `CD3D` is an easy discrimination. AUC 1.0 shows the score is valid at the
  extremes; genes in the middle of the distribution carry real uncertainty that this
  number does not convey.
- **Broadly-induced genes can be flagged without being contaminants.** The "no cell-state
  structure" feature cannot distinguish ambient RNA from a gene genuinely expressed
  uniformly across all T-cell subtypes. Interferon-stimulated genes are the obvious case
  and should be interpreted with that in mind.
- **T-cell-only input weakens the design.** With no non-T-cell population in the matrix,
  contamination cannot be contrasted against its source compartment. The Myeloid and
  Epithelial datasets from the same collection could supply that contrast and would
  strengthen a future version.
- **Anchor sets are curated by symbol** and encode an assumption about what is
  lineage-foreign. They are literals in the script and can be audited or edited.
