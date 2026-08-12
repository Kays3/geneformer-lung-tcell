# Checkpoint (ICI) and CAR-T engineering gene perturbation

Extends the targeted-panel in-silico perturbation screen
(`../perturbation_workflow/targeted_panel/`) with two applied readouts:
immune-checkpoint-inhibitor (ICI) targets and CAR-T product-engineering
edits, plus a STRING interaction-network overlay. All numbers below are
read directly from the tables in `tables/` — none are hand-typed.

## Scope and hard data limits

Two requests could not be answered as literally asked, and are recorded here
rather than silently narrowed:

- **PD-L1 (CD274) is absent from the screen.** The atlas is T-cell-only
  (CELLxGENE "T cells"); PD-L1 is expressed on tumour cells and antigen
  presenting cells, not T cells, so it was never tokenized as a source gene.
  Shown in the network figure as an explicit greyed "not in T-cell atlas"
  marker rather than omitted.
- **CAR-T tumour-antigen targets (DLL3, SEZ6, NCAM1, CD276, CEACAM5) cannot
  be perturbed here.** They are surface antigens on the tumour, not on the
  T cell being engineered; deleting them in a T cell has no biological
  meaning and none are in the 50-gene screened panel. What the screen *can*
  address is the other half of CAR-T design — edits to the T-cell product
  itself (checkpoint knockouts, persistence genes) — and that is what
  `cart_engineering_perturbation.csv` covers.
- **The screen contains no gene-to-gene edges.** It measures each gene's
  effect on cell state, not regulatory relationships. The STRING network
  (`string_network_edges.csv`) is published interaction evidence overlaid
  with the screen's per-gene effects — prior knowledge, not an inferred
  network from this dataset.

## Detection is a first-class confound, not a filter to apply silently

`delete_n` (renamed `detection_n_cells` in these tables) is the number of
source cells in which the perturbed gene was detected — i.e. the number of
cells deletion actually acts on. It is source-state specific: the same gene
has a different detection count in `sclc_to_normal` vs `luad_to_normal` vs
`normal_to_sclc`, because the three source compartments differ in size
(2,424 / 6,387 / 566 cells respectively) and in what each disease's T cells
express. Every table here joins detection on `(gene, comparison)`, not on
gene alone — an earlier draft joined on gene alone and silently applied one
comparison's counts to all six; this version corrects that.

Two consequences that change how the effect-size tables should be read:

1. **Effect size and detection are strongly anti-correlated.** Across the
   50-gene panel in `sclc_to_normal` (`cart_overexpression_vs_deletion.csv`),
   Spearman ρ = −0.600, p = 2.7×10⁻⁵ (n=42 testable genes) between
   detection count and `|delete_shift|`. The largest nominal effects come
   from the sparsest genes (the single-cell-detected hit had shift ≈0.028;
   the true top hit by magnitude, MGP, is detected in 6 cells). None of
   these are reliable estimates.
2. **8 of 50 genes have no deletion result at all** (MATR3, MMP12, PLA2G2A,
   POU2F3, PRICKLE4, S100A7, SFTPA1, YAP1) — undetected in the source
   compartment, so deletion is undefined rather than zero. Overexpression
   is still computable for these (a token can be inserted where none was
   present), so the two perturbation arms are not symmetric in coverage.
   Marked `deletion_testable=False` / "n.d." in tables and figures rather
   than dropped.

A 100-cell detection threshold (`low_detection_lt100`) is used throughout to
separate genes whose effect estimate rests on enough cells to trust from
those that don't; 18 of 50 genes fall below it in `sclc_to_normal` alone,
and the low-detection set differs by comparison (e.g. HAVCR2/CTLA4/PDCD1/
LAG3 all drop below 100 cells specifically in the small `normal_to_*`
compartments).

## Results

**Four checkpoint knockouts replicate with full donor agreement, adequate
detection, and FDR<0.05 in both arms** for SCLC→Normal — TIM-3 (HAVCR2),
TIGIT, CTLA-4, IL7R:

| gene  | delete_shift | overexpress_shift | detection (cells) |
|-------|-------------:|-------------------:|-------------------:|
| HAVCR2 (TIM-3) | +0.001868 | −0.006009 | 202 |
| TIGIT          | +0.001850 | −0.008230 | 438 |
| CTLA4          | +0.001394 | −0.001444 | 282 |
| IL7R           | +0.001349 | −0.005356 | 1,131 |

TIM-3 and TIGIT are effectively tied on deletion shift (≈1% apart) — TIM-3
is marginally first, not TIGIT "by a wide margin" as an earlier draft of
this analysis stated before correction.

**Ranked among immune/lineage genes only** (excluding ambient/technical
transcripts such as surfactant, haemoglobin, and S100 genes that dominate
the unfiltered top ranks), the four candidates sit mid-table: TIM-3 5th,
TIGIT 6th, CTLA-4 8th, IL7R 9th of the 24 immune genes detected in ≥100
cells. Real, replicated, adequately detected — not standout.

**The checkpoint programme's direction reverses between goals.** Comparing
overexpression toward Normal against overexpression toward LUAD for the
same genes shows a sign flip: overexpressing TIGIT in SCLC T cells moves
them away from Normal (−0.0082) but strongly toward LUAD (+0.0256). This is
internally consistent — the reciprocal comparison (overexpressing in LUAD
T cells moving away from SCLC) agrees in direction for 5 of 6 checkpoints
tested — and implies an ordering **Normal < SCLC < LUAD** on the
checkpoint/exhaustion axis in this model, with LUAD placed *further* along
than SCLC. This runs against a naive assumption that SCLC is the more
immune-evasive disease and is worth flagging explicitly rather than
smoothing over.

**STRING network** (`string_network_edges.csv`): started from 99 edges at
combined score ≥0.4 among the 5 ICI genes plus 11 T-cell functional/
exhaustion context genes; 45 were text-mining-only and dropped, leaving 54
edges with experimental, database, or co-expression support. PD-1, CTLA-4,
LAG-3, and TIGIT are mutually connected in this filtered set; TIM-3 has no
non-text-mining edge to any other checkpoint gene despite being the
strongest individual hit — a real asymmetry, not a plotting choice.

## Figures

- `figures/ici_cart_perturbation_network.png` — three panels: (a) all 5 ICI
  targets across 6 comparisons with donor-replication tiers, (b) all 11
  CAR-T engineering genes ranked by SCLC→Normal deletion shift, with the 4
  replicated knockouts (TIM-3, TIGIT, CTLA-4, IL7R) highlighted, (c) the
  filtered STRING network coloured by SCLC→Normal deletion shift (not
  overexpression — the colourbar label on the figure is authoritative).
- `figures/cart_overexpression.png` — (a) detection vs. |effect| for all 50
  screened genes, showing the ρ=−0.60 confound; (b) all 31 immune/lineage
  genes ranked by deletion shift, every gene labelled with its detection
  count, low-detection and undetected genes marked distinctly.
- `figures/perturbation_networks.png` — the same STRING map rendered three
  times (deletion→Normal, overexpression→Normal, overexpression→LUAD) on
  one shared colour scale and node layout, node area = detection count, so
  only colour differs panel to panel — this is what shows the Normal→LUAD
  sign flip directly.

## Tables

| file | rows | grain |
|---|---:|---|
| `ici_target_perturbation.csv` | 30 | 5 ICI genes × 6 comparisons |
| `cart_engineering_perturbation.csv` | 66 | 11 CAR-T engineering genes × 6 comparisons |
| `cart_overexpression_vs_deletion.csv` | 50 | full screened panel, sclc_to_normal only |
| `network_node_perturbation.csv` | 16 | STRING network genes, all 4 readouts |
| `string_network_edges.csv` | 54 | filtered (non-textmining) STRING edges |

## Reproducing

Regenerated from `../perturbation_workflow/targeted_panel/results/
targeted_panel_delete_overexpress_merged.csv` and
`targeted_panel_donor_consistency.csv`, plus a live STRING API call
(`string-db.org/api/tsv/network`, species 9606, required_score=400) for the
16 network genes. No cached intermediate; re-running re-fetches STRING,
which is a stable public reference and not expected to change edge
composition meaningfully between runs, though scores may shift slightly
with database updates.
