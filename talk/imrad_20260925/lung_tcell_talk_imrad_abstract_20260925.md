# Abstract (structured, ~300 words) — IMRaD restructure, 2026-09-25

**Working title:** A Foundation-Model T-cell Dysfunction Screen: Translational Candidates, Then a
Full Audit

**Background:** Foundation models fine-tuned on single-cell transcriptomes can be probed for
candidate disease-state drivers via in-silico perturbation, but a classifier trained on real tissue
also learns whatever else correlates with its labels — ambient RNA, detection artifacts, cohort
composition — and a screen alone cannot say which one it found.

**Methods:** A donor-disjoint SCLC/LUAD/normal T-cell classifier (91.9% accuracy, macro F1 0.903 --
a third of that macro-F1 from a single held-out normal-class donor, 566 cells; 46,140 cells, 42
donors overall) supported an in-silico deletion/overexpression screen, scored by significance and by
bidirectional concordance between the two operations. A screen-level ambient-risk audit and a
donor-weighted reanalysis were applied after the screen.

**Results:** Four originally identified checkpoint/persistence candidates (TIGIT, TIM-3, CTLA-4,
IL7R) replicate and validate independently in spatial tissue data (antigen-presentation programme,
rho=0.361, 7.95 sigma above the null). A checkpoint-axis ordering the original work had already
flagged as an open, pending question is completed here: donor-level SCLC-vs-LUAD shows no
significant difference (p=.216, two-sided Monte Carlo, 100,000 replicates, 19 v 22 donors) once a
single donor's 74.9% cell share is accounted for. A subsequent, expanded screen's two most
significant hits by FDR (S100A8, S100A9; FDR~1e-224) fail bidirectional concordance in all six
comparisons tested; screen-wide, half of the top 120 candidates by effect are ambient-flagged, a
curated ambient anchor, or both. One candidate, NBEAL1, survives every check applied. Model precision
(bf16 vs fp32) replicates cleanly (rho>0.998) and explains none of these results; the pre-registered
104M sign-agreement gate nonetheless still reads FAIL, final, because a bitwise-deterministic
fp32-vs-fp32 rerun made its noise floor zero rather than a usable sign threshold -- a floor
technicality, not evidence that bf16 disagrees with fp32, and the domain reviewer's ruling left the
frozen gate unamended.

**Discussion:** A perturbation screen can support real translational candidates and still require
gene-by-gene auditing for confounds; significance and effect size alone are not sufficient evidence
of a cell-intrinsic driver, and an ordering claim built on a cell-weighted cross-cohort comparison
needs a donor-level check before it is reported either way.

(323 words as drafted, counted directly rather than restated -- now over the 250-300-word default the
2026-09-24 abstract used, after the bf16 sentence was brought to slide 15's fuller construction
(2026-09-25 fix: dropped "unrelated", which is not in the source and understated that the gate reads
FAIL, final, unamended). Not trimmed to force the count back under 300 -- accuracy took the word count
over the chosen default, and that is a consequence for the human's venue/word-limit decision to see,
not something to compress back out quietly. Author list: Kaisar
Dauyey, Shinji Nakaoka — matching `poster_final` and `JSDP_P25_talk` exactly; no AI-co-authorship
line, pending the human's ruling, per the standing instruction. Venue, duration and any abstract word
limit remain unconfirmed with the human.)
