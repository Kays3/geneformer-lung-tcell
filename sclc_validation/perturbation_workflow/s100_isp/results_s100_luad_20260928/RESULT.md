# S100 LUAD ISP result (2026-09-28), cleared by Stanley 17:34 JST

**Status: `negative` in both LUAD contrasts.** The four ambient-high S100 genes (S100A2, S100B, S100A13,
S100PBP) did not score higher against their own detection-matched controls than the four ambient-low genes
(S100A4, S100A6, S100A10, S100A11). Test: exact two-sided Mann-Whitney U, 4 v 4, all 70 splits enumerated on
midranks. p = 6/70 = 0.086 for LUAD -> normal and p = 34/70 = 0.49 for LUAD -> SCLC. At 4 v 4 only complete
separation (p = 2/70) can pass. Stability gate: not_positive in both contrasts (leave-one-control-out 0/160
positive; the bootstrap cannot change a non-positive unresampled test). Stanley reproduced all of this independently.

The observed ordering (the ambient-low group above) is descriptive only. It is not `opposite_direction`, which
needs p <= 0.05. The ambient-low genes' high Q (50-98) was not tested and supports no claim.

Descriptive only, not tested: Spearman rho of Q against ambient risk = -0.539 (LUAD -> normal) and -0.253
(LUAD -> SCLC).

**No-op: PASS.** 0.0 on all 7,845 per-cell values across the 10 eligible panel genes: the two forward passes
were identical, gene by gene. This does not measure bf16 error; the 0.00072 numerical floor comes from the separate
bf16-vs-fp32 canary.

**Numerical floor:** 147 of 608 completed arms have |donor-balanced shift| <= 0.00072. Two of them are primary-gene
arms: S100A11 delete -> normal and S100B overexpress -> normal. No exclusion is registered, so none is applied.
Four primary genes (S100B, S100A13, S100PBP, S100A11) have 4-14 controls within the 0.00072 bf16 floor of their own E, so their Q is imprecise; treating every such near-tie as either way leaves both contrasts negative (Stanley, 17:35 JST: most favourable extreme p = 52/70 and 68/70, opposite extreme 6/70 and 18/70).

**Retained rows:** 616 rows, none dropped. 608 are `eligible_completed`. 8 are `not_estimable_cell_count`:
S100P (40 cells) and S100A16 (48 cells), all 4 rows each, below the 50-cell gate. S100A8 and S100A9 are anchors,
so they have E but no Q. Raw E is in E_by_gene_contrast.csv as a diagnostic only.

GPU time 4.165 h. Code 295aa2c; tables 44abd57; provenance 0c31e99 (see PROVENANCE.md).
