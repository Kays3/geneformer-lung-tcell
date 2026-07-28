# Methods and reproducibility

## 1. Source data

GSE263196: 5 fresh-frozen SCLC 10x Visium samples (GSM8187469-GSM8187473,
labelled SCLC3/SCLC4/SCLC8/SCLC9/SCLC12), already downloaded and file-integrity
audited by `sclc_validation/audit/audit_sclc_data.py`. No normal or LUAD
control tissue is available in this cohort -- it is used only as an
orthogonal spatial check on the SCLC arm, not as a training or comparison
dataset.

## 2. Spot filtering

Each sample's raw Space Ranger output (`matrix.mtx.gz`, `features.tsv.gz`,
`barcodes.tsv.gz`, `tissue_positions_list.csv.gz`) is loaded directly (no
`spaceranger`/`squidpy` directory structure required). Spots are kept if
`in_tissue == 1` and total UMI counts >= 200. Counts are normalized to 10,000
per spot and log1p-transformed before scoring.

## 3. T-cell abundance quantification -- scope decision

GSE263196 ships with no cell-type or compartment labels. Two options were
considered: (a) a marker-gene score per spot, or (b) formal deconvolution
against a multi-cell-type single-cell reference (e.g. cell2location or NNLS
against the HTAN "Immune cells" or "Combined samples" datasets). **Marker-gene
scoring was chosen** -- explicitly a proxy for T-cell abundance, not true
cell-type proportions, but consistent with the audit's own "deconvolve or
annotate" framing and avoiding a much larger reference download plus
cross-platform (scRNA-seq reference vs. Visium) deconvolution-accuracy
uncertainty. Revisit with formal deconvolution if this proxy proves
insufficient for a stronger claim.

T-cell score: `scanpy.tl.score_genes` over a pan-T-cell marker set --
`CD3D, CD3E, CD3G, CD2, CD5, CD28, TRBC1, TRBC2, IL7R, CD8A, CD8B, CD4`.

Dysfunction score: `scanpy.tl.score_genes` over the exhaustion-marker subset
of the audit's pre-registered 21-gene panel -- `PDCD1, CTLA4, HAVCR2, LAG3,
TIGIT, TOX, LAYN`. This is deliberately a different, smaller gene set from
the T-cell abundance markers, so the correlation below is not testing a
signature against a subset of itself.

## 4. Enrichment test

Per sample: Spearman correlation between the T-cell score and the
dysfunction score across all in-tissue spots, with a 95% CI via Fisher
z-transformation (`arctanh`/`tanh`, `n-3` standard error).

Per the audit's explicit guidance ("aggregate spot statistics to five
patient-level estimates; present effect sizes and uncertainty, not only
spot-level P values"), the 5 samples are combined via an inverse-variance
(Fisher z) weighted meta-analysis, rather than pooling all spots into one
test that would let the largest sample dominate and understate
between-patient heterogeneity.

## 5. Results

| Sample | Spots | Spearman ρ | 95% CI | p |
|---|---:|---:|---|---:|
| SCLC3 | 3,849 | 0.028 | [-0.004, 0.059] | 0.086 |
| SCLC4 | 2,709 | 0.154 | [0.117, 0.190] | 8.9e-16 |
| SCLC8 | 3,030 | 0.070 | [0.035, 0.106] | 1.0e-4 |
| SCLC9 | 3,519 | 0.404 | [0.376, 0.431] | 4.1e-138 |
| SCLC12 | 2,525 | 0.116 | [0.077, 0.154] | 5.8e-9 |
| **Pooled (5 samples)** | 15,632 | **0.161** | **[0.146, 0.176]** | ~0 |

4 of 5 samples show a significant positive correlation individually; SCLC3
does not reach significance on its own (CI crosses zero) but is directionally
positive and does not pull the pooled estimate negative. Effect sizes are
weak-to-moderate and heterogeneous across patients (0.03-0.40) -- this is
reported as heterogeneity, not smoothed over by the pooled estimate alone.

## 6. Limitations

- **Marker-score proxy, not deconvolution.** "T-cell abundance" here is a
  transcriptional signature score, not an estimated cell-type proportion.
  Ambient RNA, spot-level cell mixing (Visium spots typically contain
  multiple cells), and marker specificity all affect it.
- **No normal/LUAD spatial control.** This cohort cannot test whether the
  same enrichment pattern is SCLC-specific or a general tumor-immune
  phenomenon; it only tests presence/absence of the pattern within SCLC.
- **Correlational, not causal.** A positive correlation between T-cell
  abundance and dysfunction-marker expression is consistent with, but does
  not establish, in-situ T-cell exhaustion driven by the tumor
  microenvironment -- it is equally consistent with T-cells simply being
  more detectable (proportionally) in regions where they carry an exhausted
  phenotype to begin with.
- **Between-patient heterogeneity.** Effect sizes range over an order of
  magnitude (0.03-0.40) across 5 patients; the pooled estimate should not be
  read as a uniform per-patient effect.
