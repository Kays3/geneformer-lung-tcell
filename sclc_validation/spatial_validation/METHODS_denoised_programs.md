# Methods: spatial validation of the denoised perturbation programs

This is a second, independent application of the GSE263196 spatial pipeline. The
first run ([`METHODS.md`](METHODS.md)) tested a single pre-registered 7-gene
dysfunction panel. This run tests the gene programs that the all-gene in-silico
perturbation screen prioritized after technical-noise removal.

## 1. What is reused, and what is new

Reused unchanged from `spatial_validation.py`, so the two runs are comparable:

- The same 5 fresh-frozen SCLC Visium samples (GSM8187469-GSM8187473).
- The same spot QC: `in_tissue == 1`, total UMI >= 200, normalize to 10,000, `log1p`.
- The same T-cell abundance score (`CD3D, CD3E, CD3G, CD2, CD5, CD28, TRBC1, TRBC2, IL7R, CD8A, CD8B, CD4`).
- The same statistic: per-sample Spearman with a Fisher z 95% CI, combined across the
  5 patients by inverse-variance meta-analysis rather than pooling raw spots.

New in this run:

- Multiple gene programs are tested instead of one panel.
- Three guards against false positives (section 3).
- Between-patient heterogeneity is quantified with Cochran's Q and I², rather than
  described qualitatively.

## 2. Program definition

Programs are read from `primary_test_perturbation/tables/immune_cancer_candidates.csv`,
the output of the denoising step. A gene enters a program if it is concordant in at
least **2 of the 6** disease transitions, which excludes single-comparison hits. Genes
are grouped by their curated `class_label`. A program is tested only if at least
**3 genes** survive the circularity guard below.

The original 7-gene dysfunction panel is re-run alongside as a benchmark.

## 3. Guards

### 3.1 Circularity guard

Any gene shared with the T-cell abundance marker set is **removed from the tested
program**. This matters most for `T-cell identity / TCR` and `Memory / progenitor`,
which contain `CD3D`, `CD3E`, `IL7R` and similar. Without this, those programs would
correlate with T-cell abundance because they are partly the same genes, and the result
would measure nothing. Dropped genes are recorded per program in
`results_denoised_programs/program_definitions.csv`.

### 3.2 Library-size control

MHC-I genes (`HLA-A/B/C`, `B2M`) are expressed by nearly every cell type, including
tumour. A raw positive correlation with the T-cell score can therefore reflect spot
cellularity and sequencing depth rather than immune co-localization. Every program is
additionally reported as a **partial Spearman controlling for log10 total counts**,
computed on ranks. A program whose partial rho collapses toward zero is a depth artifact.

### 3.3 Expression-matched random null

With ~15,000 spots, almost any gene set reaches a small p-value, so nominal significance
carries little information. For each program, **100 random gene sets matched on mean
expression decile** are scored identically and pushed through the same meta-analysis.
This yields:

- `null_mean_rho`: what an equally-expressed random gene set achieves,
- `null_z`: standard deviations of the observed rho above that null,
- `empirical_p_vs_null`: the fraction of null draws reaching the observed rho.

`null_z` is the primary read-out of this run. A program is reported as exceeding the
null only when `empirical_p_vs_null < 0.05`.

## 4. Pipeline-equivalence check

The benchmark dysfunction panel reproduces the published first-run result exactly —
pooled rho **0.1611** with CI **[0.146, 0.176]**, and per-sample values SCLC3 0.028,
SCLC4 0.154, SCLC8 0.070, SCLC9 0.404, SCLC12 0.116, matching `METHODS.md` section 5
to three decimals. Any difference in the new programs is therefore attributable to the
gene sets, not to a changed pipeline.

## 5. Outputs

| File | Contents |
|---|---|
| `results_denoised_programs/denoised_programs_pooled.csv` | One row per program: raw rho, partial rho, null statistics, I² |
| `results_denoised_programs/denoised_programs_by_sample.csv` | Per-sample raw and partial estimates |
| `results_denoised_programs/program_definitions.csv` | Genes tested, genes dropped by the circularity guard |
| `results_denoised_programs/denoised_programs_gene_coverage.csv` | Per-sample gene presence in the Visium panel |
| `results_denoised_programs/denoised_programs_manifest.json` | Parameters, seed, program membership |
| `figures/denoised_programs_forest.png` | Forest plot: raw rho, partial rho, null mean |
| `figures/denoised_programs_null_z.png` | Program ranking by standard deviations above null |

## 6. Limitations

The limitations in `METHODS.md` section 6 all carry over unchanged, and are not
weakened by testing more programs:

- **Marker-score proxy, not deconvolution.** "T-cell abundance" remains a signature
  score, not an estimated cell-type proportion.
- **No normal or LUAD spatial control.** This cohort still cannot establish that any
  pattern is SCLC-specific rather than a general tumour-immune phenomenon. This is the
  key limitation for the MHC and interferon programs, which are plausibly pan-cancer.
- **Correlational.** Spatial co-localization is not causal evidence, and it is not
  independent confirmation of an in-silico perturbation direction — the perturbation
  screen predicts what happens when a gene is removed, which this design cannot test.
- **Between-patient heterogeneity.** Reported explicitly as I² per program.

Additional limitation specific to this run:

- **Programs are not independent.** MHC, interferon and antigen-processing gene sets
  share biology and overlap in expression; their results should not be read as separate
  confirmations of one another. No multiple-testing correction across programs is applied,
  because the matched null already sets the relevant bar per program.
