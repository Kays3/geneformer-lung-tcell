# Sources note — IMRaD-restructured oral presentation, 2026-09-25

**Timestamp convention (human's standing ruling, 2026-09-25):** every timestamp in this note and the
deck is given JST (UTC+9) first, with the UTC/source form in brackets — e.g. `18:57 JST (09:57Z)`.
Git commit times in this repo are recorded `+0900`, i.e. already JST; the bracketed form is the
derived UTC. Calendar-only dates (no time-of-day) are unaffected.

Generators: `lung_tcell_talk_imrad_figures_20260925.py` (stages the six JSDP image assets; the four
2026-09-24 figures are referenced directly, not copied), `make_lung_tcell_talk_imrad_20260925.py`
(the deck), `make_lung_tcell_talk_imrad_abstract_20260925.py` (the abstract, added 2026-09-25 -- see
below). Every number is pulled from source at build time via `git_show(ref, path)` — nothing is
typed from memory. Base branch: `origin/main` at `bfeb82fc7437c4b991fe89524bffd777136e9245`,
**2026-09-24 18:57 JST (09:57Z)** — fetched fresh at build time, not assumed from a prior message.

**The abstract's word count is now generator-computed, not hand-typed.** It was wrong three times in a
row -- 197 (undercounted before a caveat was added and never rechecked), 300 (stale after the bf16
sentence grew), and 323 (a fresh recount that was itself wrong: truncated by an off-by-one line range
in a throwaway `sed 6,34p` check, one line short of the paragraph's actual end, never re-verified
against the file before being written down). Michael caught the 323 figure not matching his own
independent `wc -w` (332/329) and named the structural problem: every other number in this deck is
pulled through `git_show()` at build time, and the abstract's self-reported count was the one number
still typed by hand -- the one value outside the rule the rest of the deck follows, and the one that
had been wrong three times. `make_lung_tcell_talk_imrad_abstract_20260925.py` now defines the four
section paragraphs as the single source of the abstract's text and computes the word count from that
same text at generation time (**329 words without the four section-label words, 333 with them** --
both over the 250-300-word default), so the number in the file can never again be stale relative to
the text it describes.

## Correction checklist — six numbers, five surfaces each

Per the standing rule ("a correction is not made until it has reached every artifact that repeats
the number"), these six items must agree across all five surfaces: **abstract, slides, this sources
note, figure captions, generator constants.** Checked for all six at build time; re-check this table
first if any of these six numbers is ever revised.

**A presence check is not a strength check.** The abstract's bf16 sentence (row 4) originally passed
this table's own "✓ (gate stated)" column while stating a materially weaker, wrong-implication version
of the caveat ("an unrelated ... technicality" — "unrelated" is not in the source and understates that
the gate reads FAIL, final, unamended). A tickbox for presence cannot see that. **Added a "Strongest
wording" column**: it names which surface states the caveat most completely and precisely, so the
others can be checked against that wording, not against a checkmark, the next time any of these six
is touched.

| # | Number / claim | Abstract | Slides | Sources note | Figure caption | Generator constant | Strongest wording |
|---|---|---|---|---|---|---|---|
| 1 | 91.9% / macro F1 0.903 + 1-test-donor caveat | ✓ (with caveat) | ✓ slide 5 (with caveat) | ✓ below | n/a (confusion.png is JSDP-original, uncaptioned beyond title) | n/a (cited, not computed) | **Slide 5** — full line citation (METHODS.md 57–67, 112) and the raw cell/donor counts; abstract necessarily compresses this |
| 2 | S100A8/S100A9 FDR + concordance failure | ✓ (with caveat) | ✓ slide 12 (with caveat) | ✓ below | `fig_concordance.png` title states 9/12, 0/12 | `N_SAME_SIGN`, `N_CONCORDANT`, `FDR_EXAMPLE` | **Slide 12** — carries the FDR-is-certainty-not-magnitude point explicitly, which the abstract's shorter form only implies |
| 3 | T6 donor-level ordering + p=0.216 + PleuralEffusion/HTA8_2001 | ✓ (p stated, no "significant") | ✓ slide 11 | ✓ below | `fig_t6_reversal.png` (2026-09-24, unchanged) | `P_COMPLETE`, `PE_SHARE` | **Slide 11** — keeps the exact-enumeration 2/35 floor as a visibly separate test; the abstract states the Monte Carlo side only |
| 4 | bf16 + frozen-FAIL 104M gate | ✓ (gate stated, matched to slide 15's construction — fixed 2026-09-25, see below) | ✓ slide 15 (gate in same sentence) | ✓ below | `fig_bf16_canary.png` (2026-09-24, unchanged) | cited from `RESULTS_BF16.md` lines 99, 118–125, not recomputed here | **Slide 15** — "both are true at once" names the non-contradiction explicitly; keep the abstract's wording matched to this, not to a shorter paraphrase |
| 5 | 42-vs-45 donor reconciliation | n/a (abstract states 42 only, no need for 45) | ✓ slide 5 | ✓ below | n/a | cited from `METHODS.md` | **Slide 5** — names all three dual-ID donors; the abstract has no equivalent detail and doesn't need one at this length |
| 7 | Four-hits cell-set caveat (targeted panel, pre-2026-09-22 runner) | ✓ (in the four-hits sentence) | ✓ slide 9 (card, replacing the blanket "no caveat" wording) | ✓ below | n/a (`screen_b.png` is JSDP-original and deletion-only) | `PANEL_OVER_N`, `PANEL_DEL_N`, `PANEL_FIX_DATE`, all asserted | **Slide 9** — states both arms' cell counts and that the panel was not rerun; the abstract compresses to one clause |
| 6 | Any normal-class claim, per-contrast | ✓ (no blanket claim made) | ✓ slides 12, 18 | ✓ below | n/a | n/a | **Slides 12, 18** — the only surfaces that state the per-contrast split (LUAD-pairing benefit vs. SCLC's tissue skew) at all; the abstract correctly makes no claim rather than a compressed, wrong one |

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
| 9 | Results R1 | 4 hits + spatial validation; cell-set caveat on the concordance half (2026-09-25) | `talk/JSDP_P25_talk.pptx/pdf`; `targeted_panel_delete_overexpress_merged.csv` @ `6882627` (overexpress_n 2,424 constant; delete_n 202/438/282/1,131); runner fix `66d235b` (2026-09-22); no later commit to the results table |
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

## Correction 2026-09-25 23:55 JST (14:55Z) — the four hits' concordance half was scored on two cell sets

**The finding (Phyllis, 2026-09-25; ruled by Michael the same night).** Before 2026-09-22 the targeted
50-gene panel runner overexpressed each gene into **every** held-out source cell, while deletion touched
only cells where the gene's token was present. In `targeted_panel_delete_overexpress_merged.csv` @ `6882627`
`overexpress_n` is constant within each contrast (2,424 SCLC / 566 normal / 6,387 LUAD) while `delete_n`
varies per gene (HAVCR2 202, TIGIT 438, CTLA-4 282, IL7R 1,131). That constant-vs-varying pattern is the
detector, and it is visible in the table without reading runner code. The runner was fixed on 2026-09-22
(`66d235b`, paired cell sets by construction); the panel results were **never rerun** — the merged table
has a single commit. All of this is asserted at build time in the generator.

**What this deck does and does not source from the panel.** Checked at Michael's request before any edit:
the generator reads **no** targeted-panel file for any printed number. Every computed value comes from
`isp_plausibility_top_candidates.csv` @ `b99365b` (whole-genome, `genes_to_perturb="all"`, which
overexpresses only genes already present in a cell — same cell set both arms, so **9/12 same-sign, 0/6
concordant, the 60/120 ambient union and NBEAL1 are unaffected**), the T6 tables @ `9ec518c`, and
`RESULTS_BF16.md`. Panel-derived content enters only as two JSDP images, and **both are within-arm**:
`detection.png` (slide 10) plots `delete_n` against `|delete_shift|`, so **rho = -0.60 is unaffected**;
`screen_b.png` (slides 3, 9) plots the deletion shift only.

**What changed.** The phrase "four replicated hits" inherits a cross-arm call, and concordance is
load-bearing in this deck — slide 7 makes "clear both arms and reverse sign" the criterion, and that
criterion is what eliminates S100A8/A9 on slide 12. So slide 9 now states what replication does mean
(deletion FDR < 0.05, same deletion sign in all 3 SCLC donors), names both arms' cell counts, and says the
panel has not been rerun; its card no longer claims "no caveat on this slide". The abstract's four-hits
sentence carries the same clause. **No number moved and no conclusion moved** — this is a qualifier, not a
retraction, and the deletion-arm evidence is what `screen_b.png` actually shows.

**Delivered-before-known.** The three genetics-audience decks (`lung_tcell_genetics_talk_{,v2_,v3_}20260925.pdf`)
print the cross-arm comparison directly and had already been delivered to the human. They were corrected and
re-delivered first, with the reason stated; see their own sources notes.

## Word count / pacing

`pdftotext` word count over the full rendered 18-page PDF (includes eyebrows, footers, citations —
not spoken content only): **2,117 words.** At 130–150 wpm, treating this as an upper bound on spoken
content, that's roughly **14–16 minutes** — closer to a 12+3 target than the outline-stage estimate
(2,510 words / 16.7–19.3 min) suggested, because several slides' actual prose came out tighter than
their word budget. Dropping the designated cut (R8, ~150 words) brings it to 17 slides / ~1,970
words / 13–15 min. **Venue, duration and abstract word limit remain unconfirmed with the human** —
this deck is built to the JSDP-matched default, not a stated constraint; report the measured number
rather than assume it settles the venue question.

## Correction 2026-09-25 23:55 JST (14:55Z) — the four hits' concordance half was scored on unpaired cell sets

Michael's check (card `2026-09-25T14-39-23Z`) and ruling (`14-43-30Z`): no slide in this deck *prints* a
delete-vs-overexpress comparison from the targeted 50-gene panel (the ρ = −0.60 figure on slide 10 is
`delete_n` vs `|delete_shift|`, deletion arm only; `screen_b.png` is deletion-only). But the words "four
replicated hits" inherit a cross-arm call: the four are the `concordant == True` rows of the panel table
@ `6882627`, whose runner overexpressed into every held-out SCLC cell (`overexpress_n` = 2,424, constant
across all 50 genes) while deleting only in token-positive cells (`delete_n` HAVCR2 202, TIGIT 438,
CTLA-4 282, IL7R 1,131). The runner was fixed on 2026-09-22 (`66d235b`); the panel has not been rerun.
Ruled to count, because concordance is load-bearing here: slide 7 makes it the criterion and slide 12
applies it to eliminate S100A8/A9, so the four headline hits must be shown to the same standard.

**Changed:** slide 9's "No caveat on this slide" card now carries the cell-set sentence (numbers pulled
and asserted at build: `PANEL_OVER_N`, `PANEL_DEL_N`, `PANEL_FIX_DATE`, single commit on the results
table), the replication definition is stated as deletion-arm in the same frame, and the abstract's
Results sentence carries the same clause. Checklist row 7 above. **No number moves, no conclusion
moves**; the de