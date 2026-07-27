# SCLC Data Feasibility Audit

**Audit date:** 2026-07-28
**Decision:** **GO, with conditions**

## Executive conclusion

The submitted abstract can be supported with a defensible SCLC-specific validation
program, but the present LUAD/LUSC exploratory work cannot by itself substantiate
the abstract's SCLC claims.

Use the HTAN/CELLxGENE T-cell object as the primary SCLC-versus-LUAD-versus-normal
single-cell cohort. It contains 46,140 T cells with integer raw counts: 11,791 SCLC
cells from 19 labelled donors, 29,829 LUAD cells from 22 donors, and 4,520 normal
cells from 4 donors. All donors have at least 20 cells, and 16 of 19 SCLC donors
have at least 100. All 21 prespecified dysfunction, cytotoxicity, progenitor, and
SCLC-subtype genes are present in the raw feature symbols.

Use GSE263196 as an orthogonal spatial validation cohort, not as Geneformer input.
Its five 10x Visium SCLC samples contain 15,774 in-tissue spots with intact count,
feature, barcode, position, scale-factor, and low-resolution image files. All 21
prespecified genes are present. The cohort has no normal or LUAD control and the
public archive does not supply pathologist compartment labels or high-resolution
histology; cell-type deconvolution and cautious patient-level inference are
therefore mandatory.

OMIX002441 is suitable as a secondary cross-platform sensitivity cohort. Although
it contains 1,039 T cells from 11 patients and all prespecified genes are
Geneformer-tokenable, only 769 T cells are from primary tumour, five patients have
at least 20 such cells, and only two patients have at least 20 T cells in both
primary tumour and adjacent normal tissue. It should not be the main training or
paired-comparison dataset.

GSE261348 is potentially useful as a clinical GeoMx validation cohort, but the
downloaded workbooks alone are not analysis-ready for response validation. The
expression matrix has 1,738 targets across 175 AOIs; 20 of 21 prespecified genes
are present (LAYN is absent). Segment metadata has 184 rows, 16 QC flags, and nine
segments omitted from the expression matrix. No response, benefit, survival, or
outcome column is present. Clinical endpoints must be obtained and joined before
testing association with atezolizumab-based treatment benefit.

## Audit scope and acceptance rules

The audit tests whether each source can support the abstract's comparison of
SCLC, LUAD, and normal-lung T-cell dysfunction and whether it can supply an
independent spatial or clinical validation layer.

A source passes primary single-cell feasibility when it has raw nonnegative
integer counts, explicit disease and donor labels, SCLC/LUAD/normal coverage,
at least 20 T cells per retained donor, donor-disjoint splitting capability, and
gene identifiers compatible with Geneformer. A spatial source passes technical
feasibility when all matrix dimensions agree, every matrix barcode has a spatial
position, images and scale factors are present, and the signature genes are
measured. Clinical validation additionally requires a documented outcome join.

## Evidence table

| Source | Unit and scale | Strength | Material limitation | Verdict |
|---|---:|---|---|---|
| HTAN/CELLxGENE T cells | 46,140 cells; 42 labelled donors | Direct SCLC/LUAD/normal comparison, raw counts, Ensembl IDs | Only 4 normal donors; one SCLC label is `PleuralEffusion` and needs identity verification | Primary discovery and internal validation |
| GSE263196 | 15,774 spots; 5 SCLC samples | Complete Visium count and coordinate bundles | No controls; no supplied tumour masks; spots are mixtures | Spatial validation after deconvolution |
| OMIX002441 | 1,039 T cells; 11 patients | Independent platform; all 21 markers tokenable | Low and imbalanced donor-level T-cell counts | Secondary sensitivity analysis |
| GSE261348 | 175 AOIs; 1,738 targets | Pretreatment spatially resolved immune panel | Outcome columns absent; 9 segments excluded; LAYN absent | Blocked pending clinical metadata join |

## Recommended analysis design

1. Freeze the current NSCLC analysis as exploratory provenance. Do not relabel
   LUSC as SCLC and do not reuse cell-level random splits.
2. Build the primary comparison from the HTAN T-cell object. Harmonize SCLC,
   LUAD, and normal labels; verify that `PleuralEffusion` is a biospecimen label
   rather than a donor identity; and split strictly by donor.
3. Report conventional pseudobulk results and Geneformer embeddings side by
   side. Treat donor, not cell, as the inferential replicate. Balance or weight
   discovery analyses, but never oversample validation/test donors.
4. Pre-register the 21-gene validation panel. Derive exhaustion, cytotoxicity,
   progenitor/memory, and SCLC-subtype scores in discovery data without altering
   the panel after seeing spatial results.
5. Deconvolve or annotate GSE263196 spots, quantify T-cell abundance, and test
   whether the pre-registered dysfunction score is enriched in T-cell-rich
   tumour regions. Aggregate spot statistics to five patient-level estimates;
   present effect sizes and uncertainty, not only spot-level P values.
6. Use OMIX002441 to test direction and rank concordance across platforms.
7. Do not claim treatment-response validation from GSE261348 until patient
   outcomes and AOI-to-patient mapping are documented and joined.

## Claim boundary for the conference abstract

The audit supports proceeding to generate validation evidence. It does **not**
yet validate the abstract's phrases “progressive transitions,” “terminal
exhaustion-like regions,” or “candidate regulators.” Those statements require
re-running the SCLC-inclusive discovery analysis with donor-held-out validation,
quantitative trajectory stability checks, and reproducible regulator selection.
Until then, they should be presented as preliminary findings.

## Reproducibility

Run:

```bash
/Users/kaisardauyey/workspace/research1/.venv/bin/python \
  sclc_validation/audit/audit_sclc_data.py
```

Large source files are intentionally ignored by Git. Download URLs and immutable
checksums are recorded in `results/audit_summary.json`; compact audit tables and
the executable notebook remain repository-trackable.

## Sources

- HTAN/CELLxGENE collection: <https://cellxgene.cziscience.com/collections/62e8f058-9c37-48bc-9200-e767f318a8ec>
- OMIX002441: <https://ngdc.cncb.ac.cn/omix/release/OMIX002441>
- GSE263196: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263196>
- GSE261348: <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE261348>
- Geneformer: <https://huggingface.co/ctheodoris/Geneformer>
