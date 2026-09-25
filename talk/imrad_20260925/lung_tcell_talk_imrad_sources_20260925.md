# Sources note — IMRaD-restructured oral presentation, 2026-09-25

**Timestamp convention (human's standing ruling, 2026-09-25):** every timestamp in this note and the
deck is given JST (UTC+9) first, with the UTC/source form in brackets — e.g. `18:57 JST (09:57Z)`.
Git commit times in this repo are recorded `+0900`, i.e. already JST; the bracketed form is the
derived UTC. Calendar-only dates (no time-of-day) are unaffected.

Generators: `lung_tcell_talk_imrad_figures_20260925.py` (stages the six JSDP image assets; the four
2026-09-24 figures are referenced directly, not copied), `make_lung_tcell_talk_imrad_20260925.py`
(the deck). Every number is pulled from source at build time via `git_show(ref, path)` — nothing is
typed from memory. Base branch: `origin/main` at `bfeb82fc7437c4b991fe89524bffd777136e9245`,
**2026-09-24 18:57 JST (09:57Z)** — fetched fresh at build time, not assumed from a prior message.

## Correction checklist — six numbers, five surfaces each

Per the standing rule ("a correction is not made until it has reached every artifact that repeats
the number"), these six items must agree across all five surfaces: **abstract, slides, this sources
note, figure captions, generator constants.** Checked for all six at build time; re-check this table
first if any of these six numbers is ever revised.

| # | Number / claim | Abstract | Slides | Sources note | Figure caption | Generator constant |
|---|---|---|---|---|---|---|
| 1 | 91.9% / macro F1 0.903 + 1-test-donor caveat | ✓ (with caveat) | ✓ slide 5 (with caveat) | ✓ below | n/a (confusion.png is JSDP-original, uncaptioned beyond title) | n/a (cited, not computed) |
| 2 | S100A8/S100A9 FDR + concordance failure | ✓ (with caveat) | ✓ slide 12 (with caveat) | ✓ below | `fig_concordance.png` title states 9/12, 0/12 | `N_SAME_SIGN`, `N_CONCORDANT`, `FDR_EXAMPLE` |
| 3 | T6 donor-level ordering + p=0.216 + PleuralEffusion/HTA8_2001 | ✓ (p stated, no "significant") | ✓ slide 11 | ✓ below | `fig_t6_reversal.png` (2026-09-24, unchanged) | `P_COMPLETE`, `PE_SHARE` |
| 4 | bf16 + frozen-FAIL 104M gate | ✓ (gate stated) | ✓ slide 15 (gate in same sentence) | ✓ below | `fig_bf16_canary.png` (2026-09-24, unchanged) | cited from `RESULTS_BF16.md`, not recomputed here |
| 5 | 42-vs-45 donor reconciliation | n/a (abstract states 42 only, no need for 45) | ✓ slide 5 | ✓ below | n/a | cited from `METHODS.md` |
| 6 | Any normal-class claim, per-contrast | ✓ (no blanket claim made) | ✓ slides 12, 18 | ✓ below | n/a | n/a |

## Base material

`talk/JSDP_P25_talk.{pptx,pdf}` and `poster/poster_final.{html,pdf}` — both gitignored (`talk/`,
`poster/` per `.gitignore`), local filesystem artifacts, not git-citable objects; measured/quoted,
not modified. `poster_final` confirmed canonical over `poster_draft_12` (2026-09-24 finding,
unchanged). Both present the SCLC/LUAD/normal classifier (the one behind this deck) exclusively —
zero real mentions of S100A8/S100A9 in either (one case-insensitive false-positive inside
`poster_final.html`'s base64-embedded image data, confirmed by inspecting the surrounding bytes;
same noise class as the earlier LuCA/LUSC false positive).

**`rho=0.361` / `7.95` sigma (slide 9, R1) are TALK-ONLY.** Verified: 0 hits for either string in
`poster_final.html`; 4 hits each in `pdftotext talk/JSDP_P25_talk.pdf`. Cited to the talk, not the
poster, on that slide and here.

**Quoted talk text, verbatim, page/slide as marked:**
- Poster (visible text, HTML entities rendered): "The model implies Normal < SCLC < LUAD on the
  checkpoint axis. This conflicts with the clinical tumour-level picture of SCLC as cold and
  ICI-resistant, so it remains an open question rather than a reconciled conclusion."
- Talk, slide 7, under a shape reading "TESTED SINCE THE POSTER": "The states are not on a line, so
  there is no ordering to contradict... Working reading: a triangle — SCLC ≲ Normal < LUAD, which
  agrees with SCLC being cold... Test 1 of 4 — sign directions, not centroid geometry. Pending T2 and
  T3." **Symbol substitution on slide 10:** the original "≲" (U+2272, less-than-or-approximately) has
  no glyph in the deck's base Helvetica font and rendered as a broken box; substituted with "≤"
  (U+2264, less-than-or-equal) for legibility. Semantically close, not letter-for-letter identical to
  the source string — flagged here rather than silently substituted.
- Talk, slide 11: "Still open: if SCLC T cells are not exhausted, what are they — and why does ICI
  still fail?" (closing slide, D2/D3).
- Talk, slide 10: STRING mechanism captions (TIGIT/TIM-3/CTLA-4/IL7R evidence lines, "54 STRING edges
  among 16 context genes", "TOX, LAYN: no edge above threshold", "PD-L1 absent from the T-cell atlas")
  — quoted for slide 16 (R8).

**Image assets** (talk/assets/*.png, gitignored, staged by filesystem copy with a size check, not
`git show` — see `lung_tcell_talk_imrad_figures_20260925.py`): `screen_b.png` (I1/slide 3),
`spatial.png` (R1/slide 9), `detection.png` (R2/slide 10), `confusion.png` (M1/slide 5),
`network_c.png` (R8/slide 16, the designated cut), `qr.png` (title/slide 1, closing/slide 18).

## Per-slide citations

| Slide | Section | Claim | Source |
|---|---|---|---|
| 1 | Title | — | — |
| 2 | Background | (framing only) | — |
| 3 | Intro I1 | 4 original hits, tissue validation | `talk/JSDP_P25_talk.pptx/pdf` |
| 4 | Intro I2 | Poster/talk hedges, quoted | `poster/poster_final.html`; `talk/JSDP_P25_talk.pdf` p.7 |
| 5 | Methods M1 | 91.9%/F1 0.903 + 1-donor caveat; 42-vs-45 | `README.md`; `perturbation_workflow/METHODS.md` lines 57–67, 111–112, 33–39 |
| 6 | Methods M2 | ISP mechanics, six comparisons | `README.md` |
| 7 | Methods M3 | Concordance criterion, verbatim | `primary_test_perturbation/scripts/build_delete_overexpress_shift_report.py` @ `b99365b` |
| 8 | Methods M4 | Ambient/donor-weighted methodology | `ambient_risk_diagnostic.py` line 198 @ `b99365b`; `immune_axis_test/donor_robustness.py` @ `9ec518c` |
| 9 | Results R1 | 4 hits + spatial validation, no caveat | `talk/JSDP_P25_talk.pptx/pdf` |
| 10 | Results R2 | JSDP's own QC + pending-test flag | `talk/JSDP_P25_talk.pdf`, slide 7 |
| 11 | Results R3 | p=0.216 Monte Carlo (B=100,000); PleuralEffusion/HTA8_2001; 2/35 exact-enumeration floor | `t6_weighting_reconciliation.csv`, `t6_permutation_tests.csv` @ `9ec518c`; `METHODS.md` lines 29–32 |
| 12 | Results R4 | S100A8/S100A9 FDR + concordance failure | `isp_plausibility_top_candidates.csv` @ `b99365b` |
| 13 | Results R5 | Ambient audit, 60/120 union | same CSV + `ambient_risk_manifest.json` @ `b99365b` |
| 14 | Results R6 | NBEAL1 | same CSV @ `b99365b` |
| 15 | Results R7 | bf16 + frozen-FAIL gate | `committed_run_stats/*` @ `origin/main`; `RESULTS_BF16.md` @ `afa0564` |
| 16 | Results R8 (designated cut) | STRING mechanism | `talk/JSDP_P25_talk.pdf`, slide 10 |
| 17 | Discussion D1 | Synthesis | (rolls up the above) |
| 18 | Discussion D2/D3 | Limitations, closing | `perturbation_workflow/METHODS.md` lines 24–25; `talk/JSDP_P25_talk.pdf`, slide 11 |

## Designated cut

**Slide 16 (R8, mechanism/STRING overlay) is the designated cut for a shorter slot.**
`INCLUDE_R8_MECHANISM = True` near the top of `make_lung_tcell_talk_imrad_20260925.py` — set to
`False` and rebuild to drop it; every other slide's numbering, content and citations are untouched
(page numbers/footers use the `PAGES` constant, would need updating to 17 if the slide is cut for a
real submission, not just previewed).

## What was NOT independently re-verified for this deck (unchanged from 2026-09-24 where applicable)

- Tissue-of-origin skew and the single-site (`HTA8`) confirmation (Discussion, slide 18): Kevin's Q4
  live investigation, relayed via Michael, not yet a committed table in the repo — same caveat as the
  2026-09-24 deck.
- Classifier-quality/fine-tune timing numbers in `RESULTS_BF16_316M.md` prose: excluded from this
  deck, no raw table backs them.
- The 316M same-model precision reference (rho=0.9998) and the ~2,500x/3.4-orders derivation:
  recomputed directly from committed CSVs in the 2026-09-24 pass; reused unchanged here.

## Word count / pacing

`pdftotext` word count over the full rendered 18-page PDF (includes eyebrows, footers, citations —
not spoken content only): **2,117 words.** At 130–150 wpm, treating this as an upper bound on spoken
content, that's roughly **14–16 minutes** — closer to a 12+3 target than the outline-stage estimate
(2,510 words / 16.7–19.3 min) suggested, because several slides' actual prose came out tighter than
their word budget. Dropping the designated cut (R8, ~150 words) brings it to 17 slides / ~1,970
words / 13–15 min. **Venue, duration and abstract word limit remain unconfirmed with the human** —
this deck is built to the JSDP-matched default, not a stated constraint; report the measured number
rather than assume it settles the venue question.
