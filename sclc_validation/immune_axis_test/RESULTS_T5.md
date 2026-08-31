# Results — T5, genome-wide differential expression

**Status: run** (both T5a and T5b). Tests 1 and 2 were already run
([PLAN.md](PLAN.md), [RESULTS_T2.md](RESULTS_T2.md)); test 3 (embedding geometry) and
test 4 (program-level overexpression) remain outstanding. Per
[PLAN.md §10](PLAN.md#10-if-the-resolution-holds): **no poster or talk wording changes
here** — both T2 and T3 are required before that, and neither T3 nor this test resolves
the axis claim on its own.

## What ran

[`differential_expression.py`](differential_expression.py), on `thinkstation2` against
the same prepared H5AD T2 used, twice:

- `POPULATION=complete` (T5a): all **46,140 T cells** (train+eval+test — 29,829 LUAD/22
  donors, 11,791 SCLC/19 donors, 4,520 Normal/4 donors).
- `POPULATION=test_only` (T5b): the same **9,377 held-out test cells** T2 used (6,387
  LUAD/4 donors, 2,424 SCLC/3 donors, 566 Normal/1 donor).

Pairwise Wilcoxon DE (scanpy) across all ~24,540 genes, log1p-CP10k, BH-adjusted, for
all three state pairs. Runtime a few seconds per population, consistent with T2's
revised estimate (PLAN.md §9).

Outputs:
[`results/differential_expression_{complete,test_only}.csv`](results/),
[`results/state_pooled_expression_{complete,test_only}.csv`](results/),
[`results/pseudobulk_per_donor_{complete,test_only}.csv`](results/).
Analysis: [`analyze_differential_expression.py`](analyze_differential_expression.py).

## A pipeline bug this run surfaced, before any biology

`differential_expression.py`'s smoke tests (synthetic fixtures, per the README) checked
shapes and column names, not gene identity. The H5AD's `var_names` — and therefore this
script's `gene` column — are **Ensembl IDs** (`ENSG...`), with no symbol column in `var`
at all. `analyze_differential_expression.py` was written to match against **gene
symbols**: the seven hardcoded exhaustion genes and the existing ISP tables' `Gene_name`
column. On the first real run, every symbol-keyed lookup silently returned empty or
matched on coincidence (one Ensembl-ID string, `ENSG00000269242`, happens to also be a
fallback value in an older `Gene_name` column — this produced a spurious "1 ISP hit
overlaps" before the fix, not a real one).

Fixed by adding `ensembl_to_symbol()` to `differential_expression.py`, sourced from the
same per-comparison ISP stats tables (`Gene_name`/`Ensembl_ID` pairs) that
`measure_baseline_expression.py` already uses for the reverse direction. 15,820 of
24,540 genes (the ones the ISP screen ever tested) get a real symbol; the rest keep
their bare Ensembl ID. All three output tables now carry a `gene_symbol` column, and
`analyze_differential_expression.py` was updated to join on it. Every number below is
from the re-run, symbol-mapped output.

## Goal 1 — where do the curated exhaustion genes rank genome-wide (SCLC vs LUAD)?

| Gene | log2fc (SCLC vs LUAD) | padj | percentile of DE significance |
|---|---:|---:|---:|
| LAG3 | +0.464 | 3.2e-23 | **7.6%** |
| CTLA4 | +0.336 | 7.0e-14 | **11.2%** |
| TIGIT | +0.134 | 6.9e-06 | **18.8%** |
| TOX | +0.241 | 4.0e-05 | **20.1%** |
| HAVCR2 | −0.042 | 0.164 | 32.5% |
| PDCD1 | −0.228 | 0.251 | 33.9% |
| LAYN | −0.235 | 1.000 | 46.4% |

Four of seven (`LAG3`, `CTLA4`, `TIGIT`, `TOX`) sit in the top 20% of genome-wide DE
significance for SCLC vs LUAD, all four higher in SCLC (positive log2fc, SCLC is the
`group`, LUAD the `reference`). The other three are not significant. The pre-registered
program is not an arbitrary pick from the genome's tail, but it is not uniformly
supported either — the same heterogeneity T1 and T2 already found gene-by-gene shows up
a third way here.

## Goal 2 — curated vs. data-driven dysfunction score

| Program | Normal | SCLC | LUAD |
|---|---:|---:|---:|
| Curated exhaustion (7 genes), mean log1p(CP10k) | 0.123 | **0.286** | 0.260 |
| Data-driven top-50 LUAD-up genes, mean log1p(CP10k) | 2.123 | 1.802 | 2.376 |

By construction the data-driven set (the 50 genes most significantly *higher* in LUAD
than SCLC) shows LUAD > SCLC — that direction is definitional, not a finding. The
curated exhaustion program is the number that matters here, and on the **complete**
population it puts **SCLC above LUAD**, reversing T2's test-only ordering (Normal 0.180
< SCLC 0.186 < LUAD 0.242). This is not a rerun of T2 with more decimal places — it is a
different, larger, and differently-composed cell population, and the next section is
why the two disagree.

## Why Goal 2 disagrees with T2: donor composition, not a contradiction

[`check_donor_composition.py`](check_donor_composition.py) breaks the curated-exhaustion
score down per donor (output:
[`results/t5a_donor_composition_check.csv`](results/t5a_donor_composition_check.csv)).

**SCLC's held-out test split is 75% one donor**, and that donor is close to the bottom
of the full 19-donor distribution:

| SCLC donor | mean log1p(CP10k) | n cells | share of test-split SCLC cells |
|---|---:|---:|---:|
| RU426 | 0.463 | 212 | 8.7% |
| RU1080 | 0.230 | 396 | 16.3% |
| **PleuralEffusion** | **0.094** | **1,816** | **74.9%** |

Across all 19 SCLC donors, `PleuralEffusion` ranks 18th (lowest = 19th); the full range
is 0.062–0.598 with several donors (`RU1195` 0.598, `RU1108` 0.542, `RU426` 0.463) well
above LUAD's own maximum (`RU699`, 0.668, is the only LUAD donor that beats them). T2's
test-only SCLC estimate is therefore dominated by an atypically low donor that happens
to carry three-quarters of the held-out SCLC cell mass.

LUAD's four test donors are comparatively balanced (9–41% of test-split LUAD cells each)
and sit closer to the middle of LUAD's 22-donor range, so its test-only estimate is not
distorted the same way.

This is the same fragility PLAN.md §9 already flagged for Normal (a single test donor)
— it turns out SCLC's held-out estimate has a version of the same problem, just less
visible with 3 donors than with 1. **Neither the complete-population estimate nor the
test-only estimate is "the" answer here without a caveat**: the complete population
mixes literal different patients across train/eval/test with no guarantee the splits
are exchangeable (PLAN.md §6's donor-disjoint-not-resampled point), while the test-only
population is what the ISP screen actually used but is this donor-imbalanced. Read the
SCLC vs LUAD ordering on the curated program as **unresolved by either measurement
alone**, not as a second confirmation of T2's specific number.

## Goal 3 — ISP hit / DE overlap

Fraction of whole-genome ISP delete/overexpress-concordant hits that are also
DE-significant (padj < 0.05) in the matching comparison, T5b (test-only) population,
against the tested-background rate:

| Pair | n ISP hits | also DE-sig | % | background % (ISP universe) | background % (genome-wide) | odds ratio | Fisher p |
|---|---:|---:|---:|---:|---:|---:|---:|
| sclc_vs_luad | 1,356 | 786 | 58.0% | 16.06% | 10.49% | 10.0 | 4.9e-308 |
| sclc_vs_normal | 434 | 251 | 57.8% | 7.76% | 4.71% | 20.6 | 1.4e-169 |
| luad_vs_normal | 1,290 | 389 | 30.2% | 5.78% | 3.70% | 11.8 | 1.2e-196 |

All three pairs show strong, highly significant enrichment: ISP delete/overexpress
concordant hits are 10–21× more likely to also be genome-wide DE-significant than the
tested background. This is the independent, non-ISP signal comment 2 of the HK-gene
review report asked for (PLAN.md §6, T5b rationale) — a real statistical test, not the
earlier detection-fraction proxy — and it lands clearly positive. As PLAN.md notes going
in: DE agreement shows the input genuinely differs between states where ISP called a
hit; it does not by itself show the model's ISP shift is causally tracking that
difference, only that the two signals point the same way more than chance predicts.

## What this does not do

- Does not run T3 or T4. The axis claim remains **under test**, per PLAN.md §10.
- Does not settle whether SCLC or LUAD carries more curated-exhaustion transcript — see
  above; this is now an open disagreement between two measurements, not a confirmed
  number.
- Does not control for CD4/CD8 composition differences between cohorts or splits, which
  remains a live confounder for every program-level result in this module (T1, T2, and
  here).
- The complete-population donor pools are donor-disjoint, not resampled from the same
  patients (PLAN.md §6) — Goal 1's and Goal 3's genome-wide numbers inherit that same
  caveat even where they don't show the SCLC/LUAD flip Goal 2 does.
