# Biological evaluation

## Objective

Determine whether statistically robust perturbations represent intrinsic
T-cell programs rather than ambient lung RNA, doublets, donor composition, or
generic cellular stress.

## Priority 1: ambient RNA and doublet sensitivity

The current top rankings include alveolar and epithelial genes such as
`SFTPC`, `SFTPB`, `NAPSA`, `MUC1`, `PIGR`, and keratins, plus stromal or
vascular-associated genes such as `FBLN1`, `DCN`, and `ACKR1`. Their presence
in a T-cell-selected analysis is a major alternative explanation that must be
tested directly.

- Calculate epithelial, alveolar, myeloid, stromal, and erythroid marker
  burdens per cell.
- Review existing doublet scores or calculate them if absent.
- Apply an ambient-RNA correction method and retain both corrected and
  uncorrected results.
- Repeat the perturbation ranking after excluding high-burden and high-doublet
  cells.
- Test whether gene effects correlate with contamination burden or are
  concentrated in a small set of specimens.

Suggested provisional gate: retain at least 60% of the top 20 candidates with
the same direction after the prespecified decontamination and exclusion
analysis, while requiring individual promoted candidates to remain stable.

## Priority 2: T-cell subtype specificity

- Stratify effects across CD4, CD8, regulatory, memory-like, exhausted,
  proliferating, and other supported T-cell states.
- Rebuild centroids within subtype where sample size permits.
- Distinguish a shared T-cell disease program from changes caused by different
  subtype proportions.
- Require enough donors and cells within each subtype before interpreting a
  missing effect as biological absence.

## Priority 3: program coherence

- Compare each candidate with source-versus-goal differential expression and
  established T-cell activation, exhaustion, cytotoxicity, interferon, stress,
  and metabolic programs.
- Ask whether deletion changes a coherent multi-gene T-cell program rather
  than only moving a global embedding coordinate.
- Examine leading-edge genes rather than interpreting thousands of overlapping
  enrichment terms as independent discoveries.
- Consolidate redundant GO, Reactome, and KEGG terms into modules, then test
  module scores by donor.

## Priority 4: external and orthogonal replication

- Repeat the shortlist in a donor-balanced external NSCLC T-cell cohort.
- Test a second Geneformer checkpoint or another single-cell representation.
- Compare with published or available CRISPR, Perturb-seq, or loss-of-function
  evidence in human T cells.
- For a final experimental shortlist, use targeted knockdown or CRISPR and
  measure T-cell activation, exhaustion, cytokine secretion, proliferation,
  viability, and tumor-cell killing where appropriate.

## Candidate interpretation classes

Assign every reviewed candidate to one of these classes:

1. **T-cell-intrinsic and donor-stable** — suitable for replication.
2. **T-cell-subtype-specific** — potentially useful with a restricted claim.
3. **Tissue-environment-associated** — informative about the sample context,
   but not an intrinsic T-cell target.
4. **Ambient/doublet-sensitive** — remove from target prioritization.
5. **Generic stress or viability-related** — retain only with explicit
   functional guardrails.
6. **Unstable or underpowered** — defer pending more donors or detections.

## Minimum biological deliverables

- `cell_contamination_scores.csv`
- `decontamination_rank_stability.csv`
- `tcell_subtype_effects.csv`
- `candidate_interpretation_classes.csv`
- `pathway_leading_edge_modules.csv`
- `external_replication_summary.csv`
- figures showing contamination sensitivity, subtype effects, pathway modules,
  and external replication

Write future outputs under `../results/biology/`, with links to the source
cohort and exact cell-filtering definitions.
