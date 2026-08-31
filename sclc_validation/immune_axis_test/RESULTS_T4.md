# Results — T4, program-level overexpression

**Status: run.** All 273 planned GPU units completed on `thinkstation1` with
the CUDA-safe runner (`--nproc 1`): 12 primary program runs, 21 nested
exhaustion-set runs, and 240 expression-matched null runs. Every unit exited
successfully, and the analyzer found zero missing completion markers.

The raw Geneformer pickles remain on the compute host. Compact summaries are
retained in:

- [`results/t4_shift_summary.csv`](results/t4_shift_summary.csv)
- [`results/t4_donor_summary.csv`](results/t4_donor_summary.csv)
- [`results/t4_matched_null_test.csv`](results/t4_matched_null_test.csv)
- [`results/t4_nested_monotonicity.csv`](results/t4_nested_monotonicity.csv)
- [`results/t4_analysis_manifest.json`](results/t4_analysis_manifest.json)

The analysis uses the held-out test sources: Normal (566 cells, one donor),
SCLC (2,424 cells, three donors), and LUAD (6,387 cells, four donors). Cell-level
standard errors are descriptive; donor-level summaries are the biological-
replicate view. Each program has 20 matched null sets, so the smallest possible
one-sided empirical p-value is `1/21 = 0.0476`; these are screening-resolution
nulls, not high-resolution tail probabilities.

## Finding 1 — the exhaustion program moves SCLC toward LUAD, but not beyond its null

The primary program shifts (mean cosine-score change over source cells) are:

| Source | To Normal | To SCLC | **To LUAD** |
|---|---:|---:|---:|
| Normal | −0.0565 | −0.0508 | **+0.1207** |
| SCLC | −0.0210 | −0.0651 | **+0.0873** |
| LUAD | −0.0282 | −0.0280 | **+0.0504** |

The SCLC→LUAD direction is positive for all three SCLC donors (PleuralEffusion
`+0.1012`, RU1080 `+0.0508`, RU426 `+0.0366`) and remains positive in every
reported SCLC CD4/CD8 stratum. This is a robust directional observation, but
direction alone is not enough to establish a program-specific effect.

Against the expression-matched null, the SCLC exhaustion result is:

| Source→target | Observed | Null mean | Null SD | Null z | Directional p | Two-sided p |
|---|---:|---:|---:|---:|---:|---:|
| **SCLC→LUAD** | **+0.0873** | +0.0816 | 0.0265 | 0.22 | **0.4286** | **0.8095** |

The observed shift is therefore **inside the matched-null distribution**. T4
does not support the claim that co-overexpressing the exhaustion program moves
SCLC cells toward LUAD more than a same-size, expression-matched random gene set.
The primary result is a positive movement, not evidence that the movement is
specific to the curated exhaustion program.

For context, the SCLC→LUAD matched-null results for the contrast programs are:

| Program | Observed | Null mean | Directional p | Two-sided p |
|---|---:|---:|---:|---:|
| Exhaustion | +0.0873 | +0.0816 | 0.4286 | 0.8095 |
| Cytotoxicity | −0.0020 | +0.0152 | 0.3333 | 0.3333 |
| Progenitor/memory | +0.0179 | +0.0243 | 0.2857 | 0.4762 |
| SCLC subtype TF | −0.0726 | −0.0186 | **0.0476** | **0.0476** |

The subtype-TF contrast is the only SCLC→LUAD shift at the null's one-sided
resolution floor, and it is in the opposite direction. No multiplicity correction
is applied to these 36 screening comparisons.

## Finding 2 — nested exhaustion titration is near-monotone, not strictly monotone

The nested sets are ordered by held-out-test expression and contain the first
1 through 7 exhaustion genes. For SCLC→LUAD, the mean shifts are:

| Set size | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
|---|---:|---:|---:|---:|---:|---:|---:|
| Mean shift | 0.0256 | 0.0172 | 0.0231 | 0.0269 | 0.0529 | 0.0813 | **0.1008** |

Five of six steps are non-decreasing (`Spearman ρ=0.893`), but the first step
decreases. The pre-registered primary criterion of fully non-decreasing size
titration is therefore **not met**. SCLC→Normal is fully non-increasing across
all six steps, while the self-state SCLC score is not monotone; this pattern is
consistent with a directional response but does not establish a scalar dose
response.

## Finding 3 — direction is consistent across CD4/CD8 strata, with important size caveats

For SCLC-source exhaustion overexpression, the SCLC→LUAD mean shifts are:

| Stratum | Cells | Donors | Mean shift | Donor SE |
|---|---:|---:|---:|---:|
| CD4 | 1,167 | 3 | +0.1029 | 0.0211 |
| CD8 | 450 | 3 | +0.0725 | 0.0252 |
| CD4 (Treg) | 60 | 3 | +0.0238 | 0.0179 |
| Other | 747 | 3 | +0.0771 | 0.0241 |

All four strata point toward LUAD, and the corresponding SCLC→Normal shifts are
negative in all four. The `Other` category is large (747 cells) and should not
be treated as a clean subtype; the Treg category is small (60 cells). These
stratified results reduce, but do not eliminate, composition as an explanation.

## Finding 4 — other programs behave as contrasts, not interchangeable axes

The primary shifts toward LUAD from SCLC are +0.0873 (exhaustion), −0.0020
(cytotoxicity), +0.0179 (progenitor/memory), and −0.0726 (SCLC subtype TF). The
programs therefore do not share one universal response direction. In particular,
the subtype-TF set is strongly away from LUAD from SCLC and toward SCLC from
LUAD, as expected for a lineage-associated contrast rather than an immune-axis
program.

## Conclusion and limits

T4 establishes that program co-overexpression produces reproducible, source-
and subtype-stratified shifts, including a positive SCLC→LUAD exhaustion shift.
It does **not** establish that the curated exhaustion program uniquely drives
that shift: the SCLC→LUAD effect is not separated from matched expression-null
sets, and the nested titration misses the strict monotonicity criterion. The
single Normal donor remains a major limitation, and the 20-null design limits
tail resolution. No poster or talk wording was changed automatically; the
qualified T3 geometric conclusion should remain qualified pending review of
these T4 null and titration results.
