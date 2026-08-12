# Interpretation: how the four replicated hits map onto the literature

Concise notes tying each of the four detection-adequate, donor-replicated
knockouts (TIM-3, TIGIT, CTLA-4, IL7R) to published CAR-T / checkpoint
evidence, plus one tension the model result raises that the literature does
not resolve. Not a systematic review — targeted checks against the claims
this screen makes.

## TIM-3 (HAVCR2) — supported

TIM-3 knockout/blockade is an active, evidence-backed CAR-T engineering
strategy in solid tumors. A patient-derived prostate-cancer CAR-T study
found HAVCR2/TIM-3-high CD8+ cells specifically marked the poor-outcome
cell state, and dual knockout of its upstream regulators reduced the
TIM-3+ population while improving tumor control. Preclinical work
targeting TIM-3 directly in an ovarian-cancer CAR-T model reports that
TIM-3 negatively regulates T-cell responses through promoting exhaustion. Consistent with the
screen's finding that deleting HAVCR2 shifts SCLC T cells toward the
normal state and that overexpressing it moves them away from it.

## TIGIT — supported, most extensively validated of the four

TIGIT knockout is the best-replicated finding in the applied CAR-T
literature. In BCMA-CAR-T for multiple myeloma, <cite index="4-1,4-6">both TIGIT knockout and anti-TIGIT antibody enhanced CAR-T proliferation, degranulation and cytotoxicity, and TIGIT-knockout CAR-T cells reduced tumor burden and exhaustion in vivo</cite>. In CD19-CAR-T for lymphoma, <cite index="1-3">TIGIT inhibition alone improved CAR-T cell efficacy</cite>, and TIGIT was identified as <cite index="1-1">a marker of exhaustion distinguishing poor from favorable responders</cite>. This is the strongest external support of the four candidates.

## CTLA-4 — supported as a target class, less CAR-T-specific evidence

CTLA-4 blockade is clinically established in lung cancer generally — <cite index="36-1">ipilimumab and other CTLA-4 antibodies are among the ICI classes used against PD-1/PD-L1 in NSCLC</cite> — but the search surfaced comparatively little CAR-T-specific knockout literature for CTLA-4 relative to TIM-3/TIGIT; most CAR-T checkpoint-editing work in the retrieved results centers on PD-1, TIM-3, and TIGIT. By deletion shift alone CTLA-4 (0.001394) is not the smallest of the four — IL7R (0.001349) is slightly smaller — so the weaker case here rests on the literature gap, not on effect size.

## IL7R — supported, but as persistence engineering, not exhaustion reversal

IL7R is mechanistically different from the other three: it is not a
co-inhibitory checkpoint, so the finding here is a persistence gene, not
another checkpoint hit. <cite index="17-1,17-6">CAR-T cells armored with IL-7 show increased in vivo persistence</cite>, and <cite index="20-1,20-3">IL7R-mediated signaling is a well-established route to CAR-T persistence, motivating engineering strategies that supply a sustained IL7Rα signal</cite>. This supports IL7R's presence in the candidate list, but the poster should not describe it as a "checkpoint" alongside the other three.

## A tension the model result raises, not resolved here

The screen's SCLC↔LUAD reciprocal check implies checkpoint/exhaustion
programme expression is **higher in LUAD than in SCLC** T cells. The
clinical literature runs the other direction at the tumor level: <cite index="33-3,33-7">the clinical efficacy of checkpoint blockade in SCLC is far less pronounced than in NSCLC, and comparative analyses indicate SCLC is even more immunodeficient than NSCLC</cite>, with <cite index="35-3">a "cold" tumor phenotype including low MHC class I, T-cell exhaustion, and a profoundly immunosuppressive microenvironment</cite> cited as the mechanism. These are not directly comparable measurements — the clinical literature describes tumor-level immune infiltration and antigen presentation, largely orthogonal to the per-cell checkpoint-gene level in the T cells that are present — but the directionality is opposite enough that this should be stated as an open question on the poster, not smoothed into either narrative. It does not undermine the four replicated hits, which are internal-comparison results (T cell vs. T cell) independent of this tumor-level literature.

## Bottom line for the poster

Cite TIGIT and TIM-3 as the two best-supported, both mechanistically
(exhaustion checkpoints with direct CAR-T knockout precedent) and
statistically (detection-adequate, donor-replicated, ranked 5th–6th among
immune genes). CTLA-4 and IL7R are real, replicated signals but weaker on
independent support — CTLA-4 for CAR-T-specific evidence, IL7R because it
is a different mechanistic class (persistence, not exhaustion) and should
be labelled as such rather than grouped with the checkpoints.
