# Results — T6, donor-balanced robustness of the immune axis

**Status: run.** CPU-only, from committed summary tables; no raw per-cell data, no GPU,
no rerun of any model. Script: [`donor_robustness.py`](donor_robustness.py). Tests:
[`test_donor_robustness.py`](test_donor_robustness.py). Outputs: `results/t6_*`.

**This is reconciliation and robustness only.** It does not promote the axis claim. T4's
matched null (directional empirical p = .4286) and the failed strict titration leave the
axis conclusion qualified, and nothing here changes that — the main result below makes the
qualification *stronger*, not weaker. No poster or talk wording changes.

## The question

T2 and T5a disagree on the SCLC-vs-LUAD ordering of the curated 7-gene exhaustion program:

| Measurement | Normal | SCLC | LUAD | Ordering |
|---|---:|---:|---:|---|
| T2, test-only population, pooled | 0.180 | 0.186 | 0.242 | LUAD > SCLC |
| T5a, complete population, pooled | 0.123 | **0.286** | 0.260 | SCLC > LUAD |

[RESULTS_T5.md](RESULTS_T5.md) attributed this to donor composition. T6 tests that
explanation by changing exactly one thing — **the unit that carries weight** — while
holding the population fixed.

## Replication unit

**The donor, not the cell.** Cells within a donor share a patient, a sample, a
dissociation and a sequencing run; they are not independent observations of a disease
state. Every pooled number above is *cell*-weighted, so a donor contributing 1,816 cells
counts 1,816 times as much as one contributing 212. Donor-level (pseudobulk) scores give
each donor weight 1.

Donor counts: complete population 22 LUAD / 19 SCLC / 4 Normal; test-only 4 / 3 / **1**.
**Normal has a single test donor, so the test-only Normal figure has no donor-level
replication at all** and no Normal contrast is tested here.

## Permutation floors — computed before any p-value below was read

Under donor-label permutation there are C(n₁+n₂, n₁) distinct assignments; the observed
assignment and its mirror are always at least as extreme as themselves, so no two-sided
p below 2/C(n₁+n₂, n₁) exists. For a balanced k-vs-k split that is the 2/C(2k,k) form.

| Comparison | donors | assignments | min two-sided p | α = 0.05 |
|---|---|---:|---:|---|
| test-only SCLC vs LUAD | 3 v 4 | C(7,3) = 35 | **0.0571** | **not attainable** |
| complete SCLC vs LUAD | 19 v 22 | C(41,19) = 244,662,670,200 | 8.18 × 10⁻¹² | attainable |

**At test-only donor scale a two-sided α = 0.05 test is arithmetically unpassable** — a
criterion stricter than the floor is a check that cannot pass, whatever the data say.
The complete population permits a real test.

## Result 1 — the disagreement is a weighting artifact, and it reconciles

Identical rows; only the weight differs.

| Population | State | n donors | cell-weighted | **donor-level** | largest donor's cell share |
|---|---|---:|---:|---:|---:|
| complete | SCLC | 19 | 0.286 | **0.301** | 18.7 % |
| complete | LUAD | 22 | 0.260 | **0.244** | 11.6 % |
| test-only | SCLC | 3 | 0.186 | **0.321** | **74.9 %** |
| test-only | LUAD | 4 | 0.242 | **0.267** | 41.1 % |

**Donor-level, both populations give SCLC > LUAD.** The reversal exists only in the
test-only *cell-weighted* number, and the mechanism is the 74.9 % cell share of a single
SCLC donor. T2 and T5a are not in conflict about the biology; they differ because one
population lets one patient carry three-quarters of a state's weight.

Direct confirmation from leave-one-donor-out: **dropping `PleuralEffusion` moves the
test-only cell-weighted SCLC − LUAD difference from −0.055 to +0.137** — a full reversal
from removing one donor.

## Result 2 — the reconciled direction is NOT statistically significant

**Null in force for both tests:** the donor-level exhaustion score is exchangeable between
the two state labels — a donor's score carries no information about whether that donor is
SCLC or LUAD. Two-sided alternative; no direction assumed.

| Comparison | donor-level difference | p (two-sided) | method | floor |
|---|---:|---:|---|---:|
| test-only SCLC vs LUAD | +0.0541 | 0.714 | exact enumeration, all 35 assignments | 0.0571 |
| complete SCLC vs LUAD | +0.0576 | **0.216** | Monte Carlo, B = 100,000 (resolution 10⁻⁵) | 8.18 × 10⁻¹² |

**Even where a real test is possible — 19 vs 22 donors, floor twelve orders of magnitude
below α — the SCLC > LUAD direction does not separate from the null (p = 0.216).** The
between-donor spread swamps it: donor-level SD is 0.151 (SCLC) and 0.141 (LUAD) against a
difference of 0.058, roughly 0.4 SD.

So the reconciliation succeeds and the claim still does not. **Both the T2 ordering and
the T5a ordering are within donor-level noise.** Reporting either as a finding about
disease states is unsupported, and that applies to T5a's SCLC > LUAD exactly as much as to
T2's LUAD > SCLC.

## Result 3 — LODO: donor-level is stable, cell-weighted is knife-edge

Donor-level SCLC − LUAD difference, recomputed dropping each donor in turn:

| Population | donor-level range | sign stable | cell-weighted range | sign stable |
|---|---|---:|---|---:|
| complete (41 drops) | +0.041 to +0.078 | **41/41** | −0.0003 to +0.061 | 40/41 |
| test-only (7 drops) | −0.065 to +0.154 | 6/7 | −0.098 to +0.137 | 6/7 |

- **Complete, donor-level: every one of 41 drops keeps the sign and the magnitude inside a
  narrow band.** The direction is robust to any single donor even though it is not
  significant — robustness and significance are different properties, and only the second
  licenses a claim.
- **Complete, cell-weighted: dropping `RU1195` takes the difference to −0.0003**, i.e.
  across zero. The cell-weighted statistic is fragile to one donor even at 19-vs-22.
- **Test-only, donor-level: dropping `RU426` flips the sign** (−0.065). With 3 SCLC donors,
  removing one removes a third of the evidence; this is expected and is why the test-only
  population cannot settle the ordering either way.

## Result 4 — CD4/CD8 composition controls

Available for the **test-only population only**; the complete population's per-cell CD4/CD8
assignment is not in any committed table. Getting it would need raw per-cell data, which is
a scope decision, not something to take unilaterally.

**The states differ in composition.** Donor-level mean share of each donor's cells:

| State | CD4 | CD4 (Treg) | CD8 | other |
|---|---:|---:|---:|---:|
| LUAD | 54.1 % | 3.0 % | 41.5 % | 1.4 % |
| SCLC | 45.0 % | 7.6 % | 31.8 % | 15.6 % |

**`PleuralEffusion` is a compositional outlier, not merely a low-scoring donor**: 10.7 %
CD8 and 40.2 % "other", against 41.4 %/43.4 % CD8 and 1.8 %/4.7 % "other" for the two other
SCLC donors. Since exhaustion markers score far higher in CD8 and Treg cells than in CD4,
its aberrant mix is a mechanism for its low score, not a coincidence beside it.

**The direction survives stratification.** Donor-level score within each stratum:

| Stratum | SCLC | LUAD | SCLC − LUAD |
|---|---:|---:|---:|
| CD4 | 0.240 | 0.219 | +0.021 |
| CD4 (Treg) | 0.844 | 0.837 | +0.007 |
| CD8 | 0.315 | 0.278 | +0.036 |
| other | 0.701 | 0.377 | +0.324 |

SCLC ≥ LUAD in all four strata, so the donor-level difference is **not** manufactured by
the states' differing CD4/CD8 mix. Two cautions: these are 3 vs 4 donors, below the
permutation floor, so no stratum is tested; and the "other" stratum's large gap rests
heavily on `PleuralEffusion`, which supplies most SCLC "other" cells.

## A scale discrepancy between two committed pipelines, found on the way

`baseline_expression_per_donor.csv` (T2's pipeline) reproduces `t2_program_summary.csv`
exactly — exhaustion 0.1802 Normal / 0.1863 SCLC / 0.2415 LUAD, same seven genes. The T5
pipeline's `pseudobulk_per_donor_test_only.csv` covers **the same test cells and the same
genes** but yields 0.145 / 0.148 / 0.191 — about 0.79× T2's values, and not by an exactly
constant factor (0.805 / 0.794 / 0.791).

This does not affect any ordering *within* either table, and it does not affect the
conclusions above. It does mean **the two tables are not on a common scale and must never
be pooled or differenced across**; no comparison in T6 crosses them. Root cause is not
established here — it would need the two normalisation paths compared directly, which is
outside this card. Flagging it rather than fixing it.

`pseudobulk_per_donor_complete.csv` is read by `check_donor_composition.py` but is **not
committed**; only its derived 45-donor summary is. The complete-population analysis here
therefore uses `t5a_donor_composition_check.csv`.

## What this does not do

- Does not promote the axis claim, or change poster/talk wording.
- Does not retune, re-rank, or rerun any model.
- Does not test any Normal contrast (one test donor; no donor-level replication).
- Does not establish *why* SCLC donors vary so widely (donor-level SD ≈ 0.15 on a mean of
  ≈ 0.30) — that spread, not the state difference, is the dominant signal in this data.
