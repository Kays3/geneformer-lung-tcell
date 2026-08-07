# Results: spatial validation of the denoised perturbation programs

**Cohort:** GSE263196, 5 fresh-frozen SCLC Visium samples, 15,632 in-tissue spots
**Method:** [`METHODS_denoised_programs.md`](METHODS_denoised_programs.md)
**Run:** 100 expression-matched null iterations per program, seed 43

## Pipeline-equivalence check

The benchmark 7-gene dysfunction panel reproduces the published first-run result to
three decimals — pooled rho **0.1611 [0.146, 0.176]**, per-sample 0.028 / 0.154 / 0.070 /
0.404 / 0.116. The pipeline is unchanged; differences below come from the gene sets.

## Pooled results

Ordered by `null_z`, the number of standard deviations above an expression-matched
random gene set. Raw rho is the same statistic used in the first run; partial rho
controls for log10 total counts per spot.

| Program | Genes | Raw ρ | Partial ρ | Null mean ρ | null z | Emp. p | I² |
|---|---:|---:|---:|---:|---:|---:|---:|
| Cytotoxic effector | 13 | 0.413 | 0.402 | 0.009 | **9.72** | 0.010 | 99.7% |
| Antigen presentation / MHC | 13 | 0.361 | **0.384** | −0.060 | **7.95** | 0.010 | 99.7% |
| Checkpoint / exhaustion | 6 | 0.218 | 0.204 | 0.037 | 6.68 | 0.010 | 98.7% |
| Costimulation / activation | 4 | 0.252 | 0.239 | 0.019 | 5.48 | 0.010 | 99.3% |
| Interferon / inflammatory | 17 | 0.280 | 0.295 | −0.030 | 5.40 | 0.010 | 99.5% |
| Trafficking / tissue residency | 6 | 0.274 | 0.273 | −0.002 | 4.95 | 0.010 | 99.1% |
| T-cell identity / TCR | 6 | 0.181 | 0.181 | 0.005 | 4.02 | 0.010 | 98.6% |
| *Dysfunction panel (benchmark)* | *7* | *0.161* | *0.159* | *0.010* | *3.88* | *0.010* | *98.9%* |
| Immunosuppressive metabolism / TME | 3 | 0.158 | 0.160 | 0.007 | 2.98 | 0.010 | 98.9% |
| Memory / progenitor | 3 | 0.116 | 0.111 | 0.028 | 1.67 | 0.040 | 92.5% |

Two programs were skipped for having fewer than 3 genes after filtering:
`Oncogenic / tumor suppressor` (2) and `Treg / suppressive` (2).

## What the guards caught

**The circularity guard was load-bearing.** It removed 6 of 12 genes from
`T-cell identity / TCR` — `CD2, CD3D, CD3E, CD3G, CD4, CD8A` — and `IL7R` from
`Memory / progenitor`. Without it, those programs would have been scored against a
T-cell abundance signature built from the same genes, and their correlations would have
been an artifact of construction rather than a finding.

**MHC is not a sequencing-depth artifact.** This was the main risk for the program:
`HLA-A/B/C` and `B2M` are expressed by nearly every cell, so the raw correlation could
have reflected spot cellularity. Controlling for log10 total counts *increases* the
estimate (0.361 → 0.384), which is the opposite of what a depth artifact does. The
interferon program behaves the same way (0.280 → 0.295).

**The null matters more than the p-value.** Every program, including the benchmark,
reaches empirical p ≤ 0.04, so "beats the null" does not discriminate much on its own
at this spot count. The informative quantity is the ranking by `null_z`, where the two
leading programs sit 2–2.5× above the benchmark panel.

## The major caveat: this is not 5/5 replication

I² is 92–99.7% for every program. The pooled estimates conceal large disagreement
between patients:

| Program | SCLC3 | SCLC4 | SCLC8 | SCLC9 | SCLC12 |
|---|---:|---:|---:|---:|---:|
| Cytotoxic effector | 0.159 | 0.673 | 0.316 | 0.625 | 0.186 |
| Antigen presentation / MHC | 0.205 | 0.555 | 0.281 | 0.626 | **−0.022** |
| Interferon / inflammatory | 0.176 | 0.296 | 0.223 | 0.581 | **−0.013** |
| Checkpoint / exhaustion | 0.044 | 0.392 | 0.126 | 0.334 | 0.220 |

The effect is carried by **SCLC4 and SCLC9**. SCLC12 is flatly negative for the two
programs that rank highest on `null_z` after cytotoxic effector, and SCLC3 is weak
throughout. With n=5 patients, a pooled estimate driven by two of them is a hypothesis,
not a replicated result. The `Cytotoxic effector` program is the only one positive in
all five samples.

## Interpretation limits

**Immune genes co-localize with immune cells — that is partly definitional.** The
circularity guard removes exact overlap with the T-cell marker set, but it cannot
remove the broader confound: `GZMK`, `CCL5`, `NKG7` and `CTSW` are largely restricted
to T/NK cells, so a positive correlation with a T-cell abundance score is close to
expected. The high rank of `Cytotoxic effector` should be read mostly as a positive
control that the assay works, not as a discovery.

By that same logic, **`Antigen presentation / MHC` is the more interesting result**.
MHC-I is expressed broadly rather than being immune-cell-restricted, so it has much
less of this built-in advantage, yet it ranks second and strengthens under
depth control. That is consistent with the all-gene perturbation screen, where antigen
presentation emerged as the dominant program after noise removal.

**Spatial co-localization is not confirmation of a perturbation direction.** The
perturbation screen predicts what happens to a cell's state when a gene is removed.
This design tests only where genes are expressed relative to T cells. The two are
consistent here; neither validates the other's causal claim.

**No SCLC-specificity.** This cohort has no normal or LUAD spatial control, so none of
these programs can be shown to be SCLC-specific rather than general tumour-immune
biology. This applies most strongly to the MHC and interferon programs, which are
plausibly pan-cancer.

**Programs are not independent.** MHC, interferon and antigen-processing sets overlap
biologically; their results are not separate confirmations of one another.

## Figure note

In [`denoised_programs_null_z.png`](figures/denoised_programs_null_z.png), bar colour
encodes the empirical p decision rule while the dashed reference line is at z = 1.96.
`Memory / progenitor` is coloured as passing (empirical p = 0.040) while sitting left of
the line (z = 1.67); it is the one borderline program and should be treated as such.

## Standing gap

Donor-level consistency has still not been applied to the all-gene perturbation result
that defined these programs. This spatial run does not close that gap — it tests the
programs in independent tissue, but the underlying gene lists have not yet been shown
to be reproducible across the perturbation cohort's own donors.
