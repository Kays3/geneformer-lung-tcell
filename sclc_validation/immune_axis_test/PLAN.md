# Plan — replacing the Perturb-seq analogy with an overexpression test of the disease axis

**Status:** plan plus partial execution. Tests 1, 2, 3, and 5 (T5a/T5b) are run
(numbers below and in [RESULTS_T2.md](RESULTS_T2.md), [RESULTS_T3.md](RESULTS_T3.md),
[RESULTS_T5.md](RESULTS_T5.md)). T4's deterministic manifest, GPU runner, and
analysis pipeline are implemented; its ambient/stress/ribosomal sensitivity arm
is run, while the GPU phases remain unexecuted. Nothing here has changed the
poster or the talk.

## 1. The claim under test

Poster conclusion 5, and the "OPEN QUESTION" half of talk slide 7:

> The model implies Normal < SCLC < LUAD on the checkpoint axis. This conflicts
> with the clinical tumour-level picture of SCLC as cold and ICI-resistant, so it
> remains an open question rather than a reconciled conclusion.

The evidence offered for the ordering is one gene from one source state: TIGIT
overexpression moves SCLC cells toward LUAD (+0.0256) and away from Normal
(−0.0082). Two shifts from a single starting point are read as a position on a
line. That inference is only valid if the three states actually lie on a line,
which has never been tested.

## 2. Why the current hero figure cannot test it

Fig. 5 ("One deletion at a time, every gene, every cell — a Perturb-seq screen
run in silico") is an explainer: it conveys the analogy, the scale of the sweep,
a STRING overlay, and the filtering funnel. It is illustrative, not inferential,
and it is built on the deletion arm, which is the wrong instrument for this
question on three counts.

- **Deletion is detection-limited.** It can only remove a gene from cells that
  already express it. Of the 50 panel genes, only 47 are testable in LUAD, 42 in
  SCLC and 40 in Normal, and per-gene N spans 1 to 6,281 with medians of 631,
  165 and 35 respectively. The three source states are therefore measured on
  different gene sets with power differing by orders of magnitude, and detection
  anti-correlates with apparent effect at ρ = −0.60.
- **Deletion cannot ask the immune-cold question at all.** "Immune cold" is a
  claim that the exhaustion program is *absent* from SCLC. You cannot delete what
  is not there. The counterfactual that matters — *if these SCLC T cells carried
  the exhaustion program, where would they sit?* — is a gain-of-function
  question by construction.
- **One source state cannot establish an ordering.** Positions on a line need at
  least the three pairwise contrasts to agree with each other.

## 3. Why overexpression is the right instrument

The overexpression arm is **census-complete**: `InSilicoPerturber` applies no
token filter, so every held-out cell of the source is perturbed for every gene.
N is a constant per state — 6,387 LUAD, 2,424 SCLC, 566 Normal — identical
across all 50 genes and all 6 directed comparisons. The detection confound that
dominates the deletion arm is eliminated by construction, and the design is a
complete balanced factorial: **50 genes × 3 source states × 2 goals**, already
computed and committed under
[`targeted_panel/results/overexpress/`](../perturbation_workflow/targeted_panel/results/overexpress).

That balance is what makes an axis test possible: the same edit can be applied
from all three starting states and the answers compared directly.

## 4. Test 1 — axis self-consistency (run; laptop only)

`axis_consistency.py` reads the committed stats tables and applies three
predictions that any genuine one-dimensional ordering must satisfy. Outputs land
in [`results/`](results).

```bash
python3 sclc_validation/immune_axis_test/axis_consistency.py
```

### T1a — collinearity

If the states lie on a line in the order Normal < SCLC < LUAD, then a gene that
raises the coordinate must move cells **from Normal toward both** other states
(same sign, Normal is one end), **from LUAD away from both** (same sign, LUAD is
the other end), and **from SCLC toward one and away from the other** (opposite
sign, SCLC is the middle). All three must hold for the same gene.

Fraction of genes satisfying each prediction:

| Program | from Normal, same sign | from LUAD, same sign | from SCLC, opposite sign | all three |
|---|---:|---:|---:|---:|
| Exhaustion (7 genes) | **0.000** | 0.857 | 0.714 | 0.00 |
| Cytotoxicity (6) | 0.667 | 0.000 | 0.667 | 0.00 |
| Progenitor/memory (4) | 0.250 | 0.250 | 0.750 | 0.00 |
| SCLC subtype TF (4) | 0.500 | 0.750 | 1.000 | 0.25 |
| All 50 | 0.320 | 0.340 | 0.740 | 0.02 |

**1 gene of 50 satisfies all three.** For the exhaustion program specifically —
the program the claim is actually about — the "from Normal" prediction fails for
**all seven genes**: raising PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX or LAYN in a
Normal T cell moves it toward LUAD and *away from* SCLC, every time. That is the
opposite of what SCLC-in-the-middle requires.

### T1b — reciprocity

A displacement along a fixed direction must raise similarity to one state
exactly as it lowers similarity to the other, so shift(A→B) and shift(B→A) must
carry opposite signs.

| Reciprocal pair | fraction opposite sign | Pearson r |
|---|---:|---:|
| SCLC→LUAD / LUAD→SCLC | 0.78 | −0.47 |
| LUAD→Normal / Normal→LUAD | 0.84 | −0.79 |
| **SCLC→Normal / Normal→SCLC** | **0.18** | **+0.71** |

The two LUAD-involving contrasts behave like a real axis. The SCLC–Normal
contrast does not: for 82% of genes the edit moves SCLC away from Normal *and*
Normal away from SCLC. That is geometrically impossible for a displacement lying
in the SCLC–Normal plane, and it is the signature of a displacement that is
largely **orthogonal** to that contrast. Neither number is evidence about where
SCLC and Normal sit relative to one another.

### T1c — is the second goal an independent measurement?

Regressing the second-goal shift on the first, across the 50 genes:

| Source | model | slope | R² |
|---|---|---:|---:|
| SCLC | →Normal from →LUAD | −0.20 | **0.71** |
| Normal | →SCLC from →LUAD | −0.40 | **0.62** |
| LUAD | →SCLC from →Normal | −0.56 | 0.33 |

From SCLC and from Normal, most of the "away from the third state" signal is a
fixed-slope shadow of the movement toward LUAD — a geometric consequence of
rotating away from it, not an independent readout. The −0.0082 that anchors the
poster's claim is, to first order, just −0.2 × the +0.0256.

### T1d — dimensionality

One component carries **74%** of the variance in the 6-comparison response, and
it loads almost entirely on `sclc_to_luad` (0.90). Its extreme genes are HBB,
HSPA1B, EEF1G, HBA2, S100A8/S100A9 — precisely the ambient-RNA / heat-shock /
rank-abundant class the repo already flags as
[contamination candidates](../perturbation_workflow/targeted_panel/RESULTS.md).
PC2 (19%) is the antisymmetric LUAD↔Normal contrast that behaves properly.

The single comparison the axis claim rests on is the one whose variance is
dominated by the flagged gene class, and the ambient correction is still on the
"not done yet" list.

## 5. Working resolution (to be confirmed by tests 2–4)

The three states form a **triangle, not a line**. There is one dominant
exhaustion direction in the model's embedding and it points at **LUAD**; SCLC
and Normal both sit at the low end of it and are separated from each other along
a different, roughly orthogonal direction. Reading two shifts from a single
source as a position on a line collapsed that triangle onto one axis and
produced an ordering that no other source state agrees with.

If this holds, the model's position is **SCLC ≲ Normal < LUAD on the exhaustion
coordinate** — the T cells present in SCLC are *not* checkpoint-high. That
**agrees** with the immune-cold phenotype rather than contradicting it, and it
points at a specific mechanism: exhaustion requires chronic antigen stimulation,
SCLC's low MHC-I means it is not delivered, so SCLC T cells are closer to
ignorant than to exhausted. Different mechanism of ICI resistance, same clinical
outcome.

Note what this does **not** touch: the four replicated hits (TIM-3, TIGIT,
CTLA-4, IL7R) are SCLC→Normal internal comparisons and stand or fall on their
own concordance, donor agreement and detection evidence. The axis reading is a
separate claim layered on top of them.

## 6. Tests 2–5 — what is still needed

Ordered by cost. Each is worth running independently of the others.

### T2 — measured baseline expression (no GPU, ~2 h, needs the atlas) — **run, see [RESULTS_T2.md](RESULTS_T2.md)**

The control the whole argument depends on and which has never been run: **what
is the actual per-cell checkpoint expression in SCLC vs LUAD vs Normal T cells
in this atlas?** Compute per-donor mean expression and detection rate for each
of the four programs, on the held-out test cells, stratified by donor and by
CD4/CD8.

This decides whether the model is describing the data or distorting it. If SCLC
T cells genuinely carry less exhaustion transcript than LUAD T cells, the model
is right about its input and the tension is between two different measurements —
a tumour-level immune phenotype and a per-T-cell transcript level — not between
the model and the clinic.

Runs on the compute host against `heldout_test.dataset`; needs no perturbation
output.

### T3 — embedding geometry (run; no GPU, see [RESULTS_T3.md](RESULTS_T3.md))

Load the three state centroids from
`state_embeddings/training_donor_disease_centroids.pkl`, compute the triangle
they span (pairwise cosine distances and interior angles), and project each
gene's overexpression displacement into that plane. Two shifts per gene give two
projections, which is exactly enough to place the displacement in the plane once
the centroid geometry is known.

This turns T1's sign tests into an explicitly assumed in-plane reconstruction
and produces the figure's central panel. It also lets T1c be compared against
specified displacement models. Centroid geometry alone does **not** determine a
unique complement slope: the slope depends on the displacement covariance.
T3 therefore reports full-space and centroid-plane isotropic references beside
the fitted −0.20 / −0.40 / −0.56, rather than calling either one the prediction.

The three-vector export is a few KB and is retained as a portable `.npz` beside
its provenance manifest, so this analysis and the figure can be rebuilt on the
laptop. The export contains aggregate vectors only; no cell- or donor-level
records are included.

### T4 — program-level overexpression (GPU runtime unprofiled; implementation ready)

Single-gene edits are small (|shift| ~0.01–0.03) and each gene carries its own
idiosyncrasies. The axis claim is about a *program*, so perturb the program.

`InSilicoPerturber` already takes a list for `genes_to_perturb` with `combos=0`,
which perturbs the set together — `run_targeted_panel.py` passes a
single-element list today, so this is a parameter change, not new machinery.

- **Program overexpression.** All 7 exhaustion genes as one set, from all 3
  sources. Plus cytotoxicity, progenitor and the SCLC-subtype-TF set as
  contrasts. 4 programs × 3 sources = 12 runs.
- **Dose–response by nested set size.** Geneformer's overexpression promotes a
  token to the front of the rank list; there is **no dose parameter**, so a
  graded dose is not available. The supported substitute is a size titration:
  overexpress the 1, 2, 3, … 7 highest-ranked exhaustion genes as nested sets and
  test whether the shift toward LUAD increases monotonically with set size. A
  real axis gives monotone dose–response; an artifact need not.
- **Matched null.** For each program, expression-matched random gene sets of the
  same size, so "the shift is big" can be read against a null of the same rank
  profile. The spatial validation already uses an expression-matched null; reuse
  that construction.
- **Subtype stratification.** Composition is a live confounder — if SCLC and
  LUAD source cells differ in CD8/CD4/Treg proportions, an apparent state axis
  may be a subtype axis. Report every program result split by CD4/CD8, which the
  poster's future-directions list already calls for.
- **Ambient sensitivity.** Re-run T1a–T1d with the flagged ambient/stress/
  ribosomal genes excluded, and report both. Given PC1's loadings this is not
  optional.

**Implementation status (30 August 2026).** The manifest fixes seed 43, 20
expression bins, 20 matched-null sets per program, and a cell-count-weighted
held-out-test expression profile shared across sources. This is 273 resumable
GPU units: 12 primary, 21 nested, and 240 null. Twenty null sets are deliberately
the minimum that permits a one-sided empirical p-value below 0.05 (minimum
1/21); treat it as a screening null, not high-resolution tail calibration.

The laptop-only ambient/stress/ribosomal sensitivity is complete. Excluding 19
flagged genes leaves 31/50: the all-three-sign 1-D consistency count remains
**1 gene** (1/50 before, 1/31 after), and PC1 variance changes only from 0.740 to
0.727. Under the operational reversal rule (at least half of remaining genes
satisfy all three sign predictions), the conclusion does **not** reverse. The
GPU program, nested-set, matched-null, donor, and CD4/CD8 results are not yet run.

### T5 — genome-wide differential expression for a data-driven dysfunction score (no GPU, needs the atlas) — **run, see [RESULTS_T5.md](RESULTS_T5.md)**

T2 answered "is the curated 21-gene panel differentially expressed" and found a
mixed picture: the exhaustion-program *average* sits Normal < SCLC < LUAD, but
only 3 of the 7 genes individually increase monotonically to LUAD (RESULTS_T2.md).
That leaves open a question T2 cannot answer because it never looked outside the
pre-registered gene list: **is "exhaustion," as this atlas actually measures it,
the right 7 genes at all**, and would an unbiased, genome-wide differential
expression test surface a different, better-supported dysfunction signature?

T5 has two variants with different cell populations, on purpose, not by
oversight — see the rationale below before the design details.

**Why two variants.** T2 (and an earlier draft of T5) restricted to held-out
test cells, mirroring the ISP pipeline's own decision boundary in
`build_primary_report.py` (only test-cell perturbations against training-only
reference centroids count as primary evidence). That boundary exists because
ISP shift is measured in the *fine-tuned model's embedding*, which saw the
training cells' labels — a real leakage concern for anything model-derived.
T5 never touches the model or its embedding; it reads raw counts directly from
the H5AD, so that leakage concern does not apply here. Checking donor
composition directly (read-only query against the prepared H5AD) found the
splits are **donor-disjoint but not resampled from the same donors** — each
split holds different individual patients:

| Disease | train | eval | test | total unique donors |
|---|---:|---:|---:|---:|
| Normal | 2 (RU682, RU685) | 1 (RU684) | 1 (RU675) | **4** |
| SCLC | 13 | 3 | 3 | up to 19 |
| LUAD | 14 | 4 | 4 | up to 22 |

Restricting to test-only means Normal's DE rests on **one donor** when three
more sit unused in train/eval — the exact "single-patient measurement" problem
§9 already flags as this module's biggest weakness, made worse rather than
better by carrying the ISP work's test-only rule into an analysis that doesn't
need it.

#### T5a — complete-atlas differential expression (primary variant, goals 1–2)

**Design.** Pairwise DE (SCLC vs LUAD, SCLC vs Normal, LUAD vs Normal) across
all ~24,540 genes on **all 46,140 T cells (train + eval + test combined)**,
`scanpy.tl.rank_genes_groups` (Wilcoxon), log1p-CP10k, BH-adjusted p-values.
Given T2's measured runtime (§9), "no GPU, needs the atlas" almost certainly
means seconds to low minutes even at this size, not hours.

**What it would decide.**
1. **Where do the 7 pre-registered exhaustion genes rank** among all genome-wide
   DE genes for SCLC vs LUAD? If they sit outside, say, the top 5% of DE
   significance, the exhaustion-program *selection* — not just the axis argument
   T1 already questioned — is itself weakly supported by this atlas.
2. **A data-driven dysfunction score**, built from the top-K genome-wide SCLC-vs-LUAD
   DE genes (whichever genes the data actually nominates, not a curated list),
   compared against the curated exhaustion program's score. Agreement would be
   reassuring; disagreement would say the curated program is the wrong lens for
   what this atlas calls SCLC T-cell dysfunction.

**What would overturn or weaken the current framing:** if the 7 exhaustion genes
rank in the bottom half of genome-wide DE significance for SCLC vs LUAD, or if
the data-driven top-K dysfunction score and the curated exhaustion score disagree
in direction, that is evidence the pre-registered program — not just T1's
reading of it — needs revisiting.

#### T5b — held-out-test differential expression (ISP cross-check only, goal 3)

**Design.** Same pairwise DE, same statistics, restricted to the same held-out
test cells T2 used (6,387 LUAD / 2,424 SCLC / 566 Normal, 1 Normal donor).
Kept deliberately narrow and separate from T5a: **only** for cross-referencing
against the whole-genome ISP hit lists (`primary_test_perturbation`'s
delete-vs-overexpress concordant hits and goal-vs-alt significant movers),
because those hit lists were themselves computed from test-only cells.
Comparing DE-on-the-complete-atlas against ISP-hits-from-test-only would mix
two different cell populations into one comparison — a real, avoidable
confound for this specific check, which is why T5a is not simply used for
goal 3 as well.

**What it would decide.** What fraction of ISP hits are *also* DE-significant
in the same cell population the hits were computed from, using a real
statistical test rather than the detection-fraction proxy the HK-gene review
report used. This is a stronger version of that proxy check, and a partial,
narrow answer to that report's Comment 2 (an independent, non-ISP signal to
compare hits against) — narrow because DE agreement shows the input differs
between states, not that the model's ISP shift is causally tracking it.

**Caveats to carry in before running, given T2's experience:**
- **Normal's single donor is a T5b problem, not a T5a problem.** T5b inherits
  it in full — any Normal-specific technical or biological idiosyncrasy in that
  one donor can populate a large share of "significant" hits for every
  comparison involving Normal, and cell-level Wilcoxon (the standard scanpy
  default) cannot distinguish that from a real donor-independent effect,
  structurally, no matter how small the p-value. T5a's four Normal donors don't
  eliminate this concern but substantially reduce it.
- **T5a's larger donor pools are a real donor-disjoint composition, not a
  resample of the same patients** — something nothing else in this module has
  stress-tested. CD4/CD8 proportions or other cohort characteristics could
  differ by split for reasons unrelated to disease state, since train/eval/test
  are different literal patients. Report T5a results per-donor, not just
  pooled, so a single outlier donor in the larger pool is still visible.
- **SCLC and LUAD support pseudobulk aggregation** in either variant (up to 19
  and 22 donors respectively in T5a, 3 and 4 in T5b), which scanpy's Wilcoxon
  does not do by default; a donor-level pseudobulk test (sum counts per donor,
  then a low-replicate test, or at minimum report per-donor direction of effect
  the way `donor_consistency_allgene.py` already does for ISP hits) is worth
  the extra step, though it cannot fix T5b's single-donor Normal side.
- **This still is not a closed loop, in either variant.** A gene both
  DE-significant and an ISP hit has two independent signals pointing the same
  way, which is more than either alone — but neither signal is a measured
  perturbation response. Real perturbation data remains the only thing that
  closes that gap (HK-gene review report, Comment 2).

## 7. The figure — Fig. 5-alt, "One gain-of-function, three starting states"

Replaces the Perturb-seq analogy band in the hero slot. Same full width. The
method explainer it currently carries is already covered by Fig. 1, so little is
lost; the funnel can move to the panel-B margin.

The figure makes one argument in four steps: *the design lets us start anywhere
→ the same edit points the same way from everywhere → so the states are a
triangle, not a line → and once you read the triangle, the model and the clinic
agree.*

**Panel A — why overexpression.** Two paired bars per source state: genes
testable by deletion (47/50 LUAD, 42/50 SCLC, 40/50 Normal, annotated with the
median N — 631, 165, 35) versus by overexpression (50/50, N constant at 6,387 /
2,424 / 566). One line: *deletion asks what a cell loses; only overexpression can
ask what a cold state would become.* Establishes the instrument in one glance.

**Panel B — the triangle (the load-bearing panel).** The three centroids in the
plane they span, sized by test-cell count. From each centroid, the measured
overexpression displacement vector for the exhaustion program (T3), plus faint
per-gene vectors behind it. The visual claim: **all three arrows point the same
way — at LUAD — regardless of where they start.** A dashed line shows the
Normal—SCLC—LUAD ordering the poster currently asserts, visibly not the geometry
the arrows describe. If T4 runs, the nested-set titration appears as graded arrow
lengths from the SCLC vertex.

**Panel C — the consistency test.** The 50 × 3 sign matrix from T1a, genes as
rows, the three predictions as columns, exhaustion program bold at the top. One
column of the exhaustion block is entirely negative. Right of it, the SCLC↔Normal
reciprocity scatter with the two "impossible" same-sign quadrants shaded and 82%
of genes sitting in them. Caption carries the single number that matters: *1 of
50 genes is consistent with a one-dimensional ordering.*

**Panel D — the reconciliation.** Three registers side by side on a common
exhaustion scale: **measured** per-donor transcript level per state (T2),
**model** position on the exhaustion direction (T3), and **clinical** tumour-level
phenotype as a separate labelled register — explicitly a different quantity, not
a third data point. Ends the panel with the resolution in one line: *the model
places SCLC T cells at the low end of the exhaustion axis, which is what "cold"
predicts; the earlier contradiction came from projecting a triangle onto a line.*

If poster space forces a cut, **B and D** carry the argument; A and C become
supplementary.

### Colour and build

Reuse the poster's system: teal for in-silico content, hematoxylin violet for
anything measured, eosin rose reserved for the clinical register in panel D so
it cannot be misread as a model output. Generator at
`sclc_validation/immune_axis_test/make_axis_figure.py`, following
`make_hero_figure.py`; every number read from `results/`, no literals in the
figure code.

## 8. What would overturn this

State it before running, so the tests are answerable either way.

- **T2 shows SCLC T cells at equal or higher exhaustion expression than LUAD T
  cells.** Then the model's LUAD-ward direction is not tracking exhaustion level,
  the resolution in §5 fails, and the contradiction stands — in a sharper form.
- **T3 shows the three centroids are near-collinear** (interior angle at SCLC
  approaching 180°). Then the 1-D reading was geometrically fair and T1's sign
  failures need another explanation.
- **T4's nested-set titration is non-monotone**, or the program shift falls
  inside the expression-matched null. Then neither the poster's axis nor §5's
  replacement is supported, and the honest position is that the shift metric does
  not measure program level at all.
- **T1a–T1d reverse when ambient/stress genes are excluded.** PC1's loadings make
  this a live possibility; it would mean the axis signal was contamination.

## 9. Caveats that constrain every result here

- **Normal rests on one test donor** (566 cells, 1 donor). Every Normal-source
  number is a single-patient measurement, including the collinearity failures in
  T1a. Report per-donor throughout and never state a Normal result without the
  donor count attached.
- **The ambient correction is a sensitivity analysis, not a correction.** The
  laptop arm excluded 19 ambient/stress/ribosomal-associated genes and did not
  reverse T1's conclusion (1/50 versus 1/31 genes satisfy all three signs).
  Treat the original and filtered magnitudes as descriptive, not as a causal
  correction.
- **Shift is movement toward a centroid, not a program score.** Two states differ
  in many ways at once, so "toward LUAD" is not by itself "more exhausted". T2 and
  T3 exist precisely to license that translation; without them, panel D's middle
  register is not earned.
- **T1 is exploratory.** It was run after seeing the poster's claim, on the same
  data the claim came from. It is strong enough to justify tests 2–4 and to stop
  the current wording from being asserted, not strong enough to be the new
  headline on its own.
- **The "~2h" estimate above for T2 was a pre-registered guess, not a profiled
  number, and it was off by roughly four orders of magnitude.** Measured runtime
  was ~1.3s. The estimate did not account for how far `H5AD[test_mask, gene_ids]`
  in backed mode shrinks the working set before anything is computed: the actual
  job touches 176 of 24,540 genes and 9,377 of 46,140 cells (a ~130x reduction),
  and the per-donor/CD4-CD8 stratified means and detection rates on that subset
  are a few thousand vectorized arithmetic operations, not a search or a model
  pass. The GPU stages elsewhere in this repo (the all-gene ISP screen, ~5 days
  on 2 nodes) are a different kind of cost entirely — one Geneformer forward pass
  per perturbed cell-gene pair, millions of them — and that distinction is why T2
  was scoped "no GPU" in the first place. Treat every remaining time estimate in
  this plan (T3, T4) with the same skepticism until profiled.

## 10. If the resolution holds

Poster conclusion 5 changes from an unresolved contradiction to a resolved one,
and the talk gains a better slide 7 — the confound control stays, the "OPEN
QUESTION" half becomes a geometric correction with a figure behind it. Slide 12's
Q&A hook (*"why does the model order Normal < SCLC < LUAD?"*) is retired and
replaced with the substantive open question underneath it: **if SCLC T cells are
not exhausted, what are they, and why does ICI still fail?** That is a better
question, and T2's subtype-stratified baseline is the first step toward it.

T2 and T3 are now run, but no poster or talk wording has been changed
automatically. T3 rejects the collinear-centroid premise, while T4's program,
nested-set, matched-null, donor, and CD4/CD8 GPU results remain outstanding; the
correct public state is therefore still a qualified geometric result, not a new
unconditional replacement claim.
