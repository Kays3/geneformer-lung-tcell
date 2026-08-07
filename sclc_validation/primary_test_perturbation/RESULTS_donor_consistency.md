# Results: donor-level consistency of the all-gene concordant hits

This closes the gap flagged in the primary all-gene report and in the spatial validation:
every concordant hit was a **cell-level** statistic, so a candidate driven by one patient
was indistinguishable from one reproducing across patients.

Script: [`scripts/donor_consistency_allgene.py`](scripts/donor_consistency_allgene.py)

## Reconstruction is exact

The all-gene run stores raw output per cell-batch, not per gene. Each file
`..._shard{S}_dict_cell_embs_{C}batch{B}_raw.pickle` holds one cell, where `C` indexes
cells in **length-descending order within the shard**. That assumption is asserted at
startup, not trusted: the perturbed-gene count for cell `C` must equal the `C`-th longest
cell's `length - offset`.

That check earned its place — it caught a real discrepancy. The offset is **2 for
deletion** but **3 for overexpression**, because overexpression skips the gene already
ranked first (moving it to the front is a no-op). Measured as exactly constant across
every cell of every shard checked in all three sources.

Re-aggregating the per-cell values to gene level and comparing against the Geneformer
stats CSVs:

| | Result |
|---|---|
| Genes compared | 165,438 across all 12 arm × comparison tables |
| Max absolute difference | **1.11e-16** (floating-point epsilon) |
| Correlation with `Shift_to_goal_end` | **1.000000** in all 12 tables |
| Detection-count mismatches | **0** |

The donor mapping is therefore provably correct, not merely plausible.

## Overall result

| Class | Hits | % |
|---|---:|---:|
| Fully consistent (all donors agree in sign, both arms) | 1,645 | 45.9% |
| Majority consistent (≥50%, both arms) | 1,433 | 40.0% |
| Single-donor only (unassessable) | 408 | 11.4% |
| **Inconsistent (discard)** | **100** | **2.8%** |

By transition:

| Comparison | Fully | Majority | Inconsistent | Single-donor |
|---|---:|---:|---:|---:|
| luad_to_normal | 630 | 606 | 7 | 0 |
| luad_to_sclc | 597 | 627 | 15 | 0 |
| sclc_to_luad | 217 | 131 | **64** | 1 |
| sclc_to_normal | 201 | 69 | 14 | 0 |
| normal_to_luad | 0 | 0 | 0 | 155 |
| normal_to_sclc | 0 | 0 | 0 | 252 |

Every normal-source hit is unassessable: the test set has exactly **one** normal donor.
`sclc_to_luad` carries most of the inconsistency, which fits its thin 3-donor base.

## The denoising step selected for donor-robust hits

This was not designed as a test of the denoising, but it functions as one. Restricting to
assessable hits (normal-source excluded, since it cannot be scored either way):

| Set | Fully consistent |
|---|---|
| Denoised immune / cancer set | **169 / 233 = 72.5%** |
| All other concordant hits | 1,476 / 2,945 = 50.1% |

Fisher exact **odds ratio 2.63, p = 2.3e-11**.

The immune/cancer genes prioritized on biological grounds are substantially more likely to
reproduce across patients than the background of concordant hits. Since donor identity
played no part in building those programs, this is independent support for the denoising
step rather than a circular result.

## The headline candidates survive

The antigen-presentation program that led both the denoised screen and the spatial
validation holds up. No MHC gene is donor-inconsistent in any transition:

| Gene | Fully consistent | Unassessable | Inconsistent |
|---|---:|---:|---:|
| B2M | 3/3 | 0 | 0 |
| HLA-E | 3/3 | 0 | 0 |
| HLA-B | 1/1 | 0 | 0 |
| PSMB9 | 4/5 | 1 | 0 |
| HLA-C | 4/6 | 2 | 0 |
| CD74 | 3/5 | 2 | 0 |
| STAT1 | 4/5 | 1 | 0 |
| TIGIT | 3/5 | 2 | 0 |

## Discard list

Only 6 of 293 immune/cancer hits are donor-inconsistent. They should not be carried
forward despite passing FDR and concordance:

| Gene | Comparison | Program |
|---|---|---|
| CD3E | luad_to_sclc | T-cell identity / TCR |
| TGFB1 | sclc_to_luad | Treg / suppressive |
| IRF7 | sclc_to_luad | Interferon / inflammatory |
| CCL3 | sclc_to_luad | Cytotoxic effector |
| CD40LG | sclc_to_luad | Costimulation / activation |
| LDHA | sclc_to_normal | Immunosuppressive metabolism / TME |

Five of six are `sclc_to_*`, consistent with that source's weaker donor base.

## The limitation that matters most

**"Fully consistent" is a much weaker bar than it sounds**, because the held-out test set
has very few donors:

| Source | Donors | Cells | Composition |
|---|---:|---:|---|
| luad | 4 | 6,387 | RU1138 2,625 · RU1170 2,086 · RU675 1,107 · RU1137 569 |
| sclc | 3 | 2,424 | PleuralEffusion 1,816 · RU1080 396 · RU426 212 |
| normal | 1 | 566 | RU675 566 |

So "all donors agree" means **4 patients for luad and 3 for sclc**. With three donors,
agreement by chance alone is not rare, which is why the enrichment comparison against the
background is more informative than the raw 72.5%.

Two further cautions:

- **The sclc base is dominated by one sample.** `PleuralEffusion` contributes 75% of sclc
  cells, and it is a sample-type label rather than a patient identifier — it may be pooled.
  sclc-source consistency should be treated as the weakest of the three.
- **RU675 appears as both a luad and a normal donor.** That is the paired tumour-normal
  design noted in the original validation report, handled correctly by the split's
  cross-disease pinning, but it means the single normal donor is not independent of luad.

## What this does and does not settle

Settled: the concordant hits are no longer cell-level only, the reconstruction is exact,
2.8% of hits are demonstrably driven by donor disagreement and can be dropped, and the
denoised immune/cancer set is measurably more donor-robust than the background.

Not settled: donor consistency across 3–4 patients is not population-level
generalization, no normal-source hit can be assessed at all, and none of this makes the
in-silico perturbation causal. Ambient-RNA correction (CellBender) remains the other open
sensitivity analysis.
