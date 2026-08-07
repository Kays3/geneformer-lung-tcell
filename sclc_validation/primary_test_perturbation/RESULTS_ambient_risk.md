# Results: ambient-RNA risk of the prioritized candidates

Methods and reproducibility: [`METHODS_ambient_risk.md`](METHODS_ambient_risk.md)

**This is a diagnostic, not a correction.** CellBender could not be run — the only
available matrix is filtered (min 113 UMI/cell, no empty droplets) and already subset to
T cells, so the ambient profile is unidentifiable. No counts were modified. What follows
scores whether candidates *behave* like contaminants, calibrated against genes whose
status is known.

## The score separates known answers perfectly

| | Value |
|---|---|
| Ambient anchors retained | 34 (surfactant, haemoglobin, myeloid, epithelial, stromal) |
| T-cell anchors retained | 36 (CD3D/E/G, TRAC, LCK, ZAP70, IL7R, GZMA, FOXP3, …) |
| **5-fold cross-validated AUC** | **1.000 ± 0.000** |
| Median risk, ambient anchors | 0.968 |
| Median risk, T-cell anchors | 0.011 |
| Flag threshold (25th pct of ambient anchors) | 0.950 |

Genes scored: 11,047 of 24,540, after requiring detection in ≥0.5% of cells.

The two anchor groups are cleanly separated, which is what licenses using the score on
candidates. It is also an easy discrimination by construction — see the limitation below.

## 109 of 110 candidates behave like genuine T-cell genes

| Flag | Genes |
|---|---:|
| ok | **109** |
| AMBIENT_RISK | **1** |

The single flag is **MX2** (Interferon / inflammatory, risk 0.960, 81st percentile,
concordant in 2 transitions).

## Antigen presentation is the *least* ambient-like program

Median ambient risk per program, ascending:

| Program | Genes | Median risk | Flagged |
|---|---:|---:|---:|
| **Antigen presentation / MHC** | 15 | **0.0011** | 0 |
| Cytotoxic effector | 16 | 0.0031 | 0 |
| T-cell identity / TCR | 14 | 0.0055 | 0 |
| Immunosuppressive metabolism / TME | 3 | 0.0113 | 0 |
| Treg / suppressive | 5 | 0.0116 | 0 |
| Memory / progenitor | 5 | 0.0237 | 0 |
| Trafficking / tissue residency | 7 | 0.0409 | 0 |
| Costimulation / activation | 8 | 0.1168 | 0 |
| Checkpoint / exhaustion | 12 | 0.1648 | 0 |
| Oncogenic / tumor suppressor | 5 | 0.1943 | 0 |
| Interferon / inflammatory | 20 | 0.2238 | 1 |

This is the result that matters most. Antigen presentation — the program that led the
denoised screen, ranked second in the spatial validation, and survived donor consistency
with no inconsistent hit — has the **lowest** ambient risk of any program, a median of
0.0011 against a contaminant median of 0.968. Its risk is below that of the canonical
T-cell anchor set.

That is also biologically coherent rather than merely reassuring: MHC-I is genuinely
transcribed by all nucleated cells, T cells included, so ambient contamination was never
the leading explanation for it. The diagnostic agrees with the biology.

## The one flag is probably a method artifact, not contamination

`MX2` is interferon-stimulated and therefore induced broadly and fairly uniformly across
cell states. The diagnostic treats "no cell-state structure" as ambient-like, so a
genuinely ubiquitous transcript is indistinguishable from a contaminant on that feature.
The interferon program's elevated median (0.224, the highest of the eleven) has the same
cause.

`MX2` should be treated as unresolved rather than discarded: this design cannot separate
"contaminant" from "uniformly induced" without a non-T-cell compartment to contrast
against.

## What this settles

The ambient concern raised in the original validation report — that hits such as
haemoglobin and epithelial genes were contamination — was already handled by the denoising
step, which excluded those classes outright. This diagnostic tests the harder remaining
question: whether the *surviving* candidates are contamination in disguise. They are not.
109 of 110 sit in the T-cell-intrinsic part of the distribution, and the headline program
sits at the very bottom of the risk scale.

## What this does not settle

- **No correction was applied.** If a candidate *were* contaminated, this would not fix it,
  and nothing downstream has been re-run on decontaminated counts.
- **AUC 1.000 overstates precision.** Surfactant versus `CD3D` is an easy call. The score
  is trustworthy at the extremes, where the candidates happen to fall, but a gene scoring
  0.5 would carry real uncertainty this number does not express.
- **A true CellBender run is still outstanding** and still requires primary CellRanger
  output from the HTAN portal, which is a data-access task rather than a compute one.
- **The design is weakened by T-cell-only input.** The Myeloid and Epithelial datasets in
  the same CELLxGENE collection would supply the missing contrast compartment and would
  let a future version distinguish "uniformly induced" from "ambient" — resolving `MX2`.
