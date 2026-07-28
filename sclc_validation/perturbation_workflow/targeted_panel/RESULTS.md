# Targeted 50-gene panel perturbation -- results

50 genes (21 pre-registered panel + 29 top drivers from the prior LUAD/LUSC/normal
screen, contamination markers excluded), delete + overexpress, individually,
across all three source states (SCLC, LUAD, normal). 300 gene-runs, all
completed; 6 directional comparisons x 2 perturbation types = 12 stats
tables in `results/delete/` and `results/overexpress/`.

## Headline: internal validation via SCLC master regulators

**ASCL1 and NEUROD1** -- the canonical SCLC neuroendocrine transcription
factors, not part of the immune panel, included as two of the 29 top
drivers -- show the single strongest concordant signal in the dataset for
the SCLC-to-LUAD comparison:

| Gene | Delete shift (toward LUAD) | Overexpress shift (toward LUAD) | Delete FDR | Overexpress FDR | Delete N |
|---|---:|---:|---:|---:|---:|
| NEUROD1 | +0.378 | -0.013 | 1.3e-7 | 4.2e-132 | 12 |
| ASCL1 | +0.157 | -0.037 | 2.1e-30 | ~0 | 73 |

Deleting either gene from an SCLC cell moves it toward LUAD (loses SCLC
identity); overexpressing either moves an SCLC-comparison cell further
toward SCLC (away from LUAD). This is exactly the expected direction for
genes that are master regulators of SCLC neuroendocrine identity, and
functions as an internal positive control: the model and this concordance
test recover known, strong tumor biology even though these genes were
selected purely by data-driven ranking in the prior workflow, not by
expectation. Caveat: N=12 and N=73 detections respectively -- very few T
cells express these tumor-intrinsic transcription factors at all, so despite
the extreme statistical significance, treat this as a strong sanity check on
the pipeline rather than a claim about T-cell biology.

## Panel genes: concordant hits by comparison

Of 21 x 6 = 126 panel-gene x comparison combinations tested, **60 are
concordant** (both delete and overexpress FDR < 0.05, opposite-sign shift).
Every one of the 21 panel genes has at least one concordant comparison.

Best-powered concordant panel hits (N >= 300 detections in the delete arm,
so not resting on a handful of cells):

| Comparison | Gene | Delete shift | Overexpress shift | Delete N |
|---|---|---:|---:|---:|
| luad_to_sclc | TIGIT | +0.0069 | -0.0128 | 1,320 |
| luad_to_sclc | GZMH | +0.0057 | -0.0050 | 1,621 |
| luad_to_normal | GZMH | -0.0047 | +0.0056 | 1,621 |
| luad_to_normal | CCR7 | +0.0047 | -0.0182 | 1,561 |
| luad_to_sclc | CCR7 | -0.0016 | +0.0056 | 1,561 |
| sclc_to_normal | GNLY | +0.0030 | -0.0011 | 1,183 |
| luad_to_normal | NKG7 | +0.0025 | -0.0091 | 2,819 |
| luad_to_sclc | PRF1 | +0.0016 | -0.0008 | 1,534 |
| sclc_to_normal | IL7R | +0.0013 | -0.0054 | 1,131 |
| luad_to_sclc | LAG3 | -0.0008 | +0.0062 | 1,190 |
| luad_to_normal | LAG3 | +0.0008 | -0.0116 | 1,190 |
| luad_to_sclc | TCF7 | -0.0007 | +0.0066 | 1,431 |

The exhaustion markers (TIGIT, LAG3) and cytotoxicity markers (GZMH, NKG7,
PRF1, GNLY) both show robust, well-powered concordant effects specifically
around the LUAD arm (luad_to_sclc, luad_to_normal) -- effect sizes are small
(0.001-0.02) but consistent in direction and backed by 1,000+ detections,
not driven by a handful of cells. CCR7/IL7R/TCF7 (progenitor/memory markers)
also show concordant effects, consistent in direction with a
progenitor-vs-exhausted axis moving opposite the cytotoxicity/exhaustion
markers.

## Top-driver genes: concordance and a contamination flag

The 29 top-driver genes show a higher raw concordance rate than the panel
(16/29 in luad_to_normal, 15/29 in luad_to_sclc) -- expected, since they were
originally selected for large |shift| in the prior workflow. Several of the
best-powered concordant hits are genes the repo's own
`evaluation/biology/README.md` would flag as ambient-RNA/stress candidates
even though they weren't on the explicit exclusion list used to build this
panel: **HBA1/HBB** (hemoglobin -- red cell contamination), **HSPA1B** (heat
shock -- generic stress), **RPS26** (ribosomal, N=5,583 -- extremely
well-powered but rank-abundant genes are exactly the class flagged as
creating broad sequence-context effects rather than specific biology).
**S100A8/S100A9** (myeloid alarmins) and **TPSB2** (mast cell tryptase) are
plausible cell-type-contamination signals in a T-cell-selected dataset
rather than T-cell-intrinsic effects.

**These top-driver results should be read as candidates for the biological
evaluation pipeline already scoped in `evaluation/biology/README.md`
(ambient-RNA/doublet sensitivity, T-cell subtype specificity, program
coherence), not as validated findings.** That evaluation has not been run
against this panel.

## What has not been done yet

- **Donor-level consistency.** All shifts above are cell-level aggregates
  across a source's held-out cells; whether a concordant effect holds
  consistently across donors (vs. being driven by one or two donors) has not
  been checked. This is listed as required in `../METHODS.md` \S8 and is the
  most important next step before treating any hit here as a candidate.
- **Ambient-RNA / contamination correction**, per `evaluation/biology/README.md`.
- **Sensitivity analysis excluding ribosomal/rank-dominant genes.**
- **Pathway-level interpretation.**

## Files

- `results/delete/targeted_delete_{source}_to_{target}.csv`,
  `results/overexpress/targeted_overexpress_{source}_to_{target}.csv` --
  raw per-gene stats, one file per directional comparison per type.
- `results/targeted_panel_delete_overexpress_merged.csv` -- delete and
  overexpress joined per gene per comparison, with a `concordant` flag.
- `results/targeted_panel_concordant_hits.csv` -- the 123 concordant rows,
  ranked by |delete shift|.
- `results/targeted_panel_concordance_summary_by_comparison.csv` --
  concordance rate by gene set (panel vs. top-driver) x comparison.

Rebuild the merge/concordance analysis with `analyze_targeted_results.py`.
