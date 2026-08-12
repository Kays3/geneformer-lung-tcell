# Oral presentation blueprint — 12 slides, 16:9

**Talk:** Foundation-model perturbation analysis identifies T-cell dysfunction programs in small cell lung cancer
**Speakers:** Kaisar Dauyey · Shinji Nakaoka — Laboratory of Mathematical Biology, Faculty of Advanced Life Science, Hokkaido University
**Source:** 24JSDP poster P25 (`poster/poster_final.pdf`)
**Audience:** clinicians and wet-lab biologists, not foundation-model specialists

---

## Numbers corrected from the briefing

Three figures in the deck brief did not match the result files. The values below are what the deck uses.

| Item | Briefed | **Correct** | Where the briefed number came from |
|---|---|---|---|
| Test accuracy | 95.1% | **91.9%** | — |
| Macro F1 | 0.887 | **0.903** | — |
| Pooled spatial ρ | 0.361 [0.328, 0.393] | **0.161 [0.146, 0.176]** | 0.361 is the antigen-presentation/MHC program, a different gene set |
| Visium specimens | 11 | **5** | 11 is the ICI/CAR-T candidate count on the STRING map |

The ρ mix-up is the one that matters most: 0.161 and 0.361 are both yours, but they measure different gene sets. Slide 9 presents both, correctly attributed — the contrast strengthens the story rather than weakening it.

## Where the built deck differs from this blueprint

`JSDP_P25_talk.pptx` is the deliverable; this document is the plan it was built
from. Three deliberate departures:

- **Slide 2 is the poster's "What is an *in silico* perturbation?" panel**, not
  the clinical-problem contrast planned below. The panel is redrawn natively
  (ranked gene list, TIGIT highlighted, delete/overexpress, the state map), so
  it stays sharp and editable. The clinical-problem framing now lives only in
  slide 1's speaker notes — add a slide back if you want it on screen.
- **Slide 3 is the hit rule**, not a second pass over the four stages, which
  would have duplicated slide 2. It covers bidirectional concordance, the
  replication requirement, and what a "shift" actually is.
- **Slide 9 shows one figure, not two.** The planned standalone forest plot is
  panel C of the Visium composite already on the slide — the same per-specimen
  effects and the same pooled ρ. The reclaimed width went to the ρ callouts and
  the programme chart.

## Timing (13–15 min talk + 5 min Q&A)

| Slides | Section | Minutes |
|---|---|---|
| 1–3 | Hook, problem, method intuition | 3.0 |
| 4–5 | Cohort and model credibility | 2.5 |
| 6–8 | Screen, confound, CAR-T | 4.5 |
| 9–10 | Tissue validation, mechanism | 3.0 |
| 11–12 | Limits, close | 2.0 |

**Design system.** Reuse the poster's palette so the deck and board read as one project: navy `#12253f` headers, teal `#0b6f6a` for *in silico* content, hematoxylin violet + eosin rose for anything measured in tissue. One idea per slide; figures at 60%+ of slide area; no slide over 4 bullets.

---

# Slide 1 — Title & Take-home

### Slide Title
**A foundation model points to antigen presentation in SCLC T-cell dysfunction — and independent tissue agrees**

### Core Narrative / Speaker Notes
- Open with the clinical stake: SCLC is the most aggressive thoracic malignancy and barely responds to checkpoint blockade, despite a high mutational burden that *should* make it visible to the immune system.
- State the one-sentence result up front so the audience knows where you are going: a foundation model, perturbed *in silico*, nominates dysfunction candidates in SCLC T cells; the programme it points to is associated with T-cell abundance in independent tumour tissue.
- Name the honest framing early — this is a **hypothesis generator validated against tissue**, not a knockout experiment. Saying it now buys credibility for everything that follows.

### Visual Layout & Graphic Instructions
- Full-bleed navy title block, top third. Title left-aligned, authors and Hokkaido affiliation beneath in smaller weight.
- Centre: the poster's take-home strip as a three-node horizontal flow —
  `model-derived dysfunction candidates in SCLC` → `antigen-presentation programme` → `associated with T-cell abundance in independent tumour tissue`.
- Bottom right: 24JSDP / P25 chip and the repo QR (same code as the poster).

### On-Slide Text
- Geneformer-V2-104M, fine-tuned on HTAN SCLC T cells
- Bidirectional *in silico* perturbation → candidate dysfunction genes
- Tested in intact tissue: 5 Visium SCLC sections, 15,632 spots
- **Hypothesis generator, tissue-validated — not a knockout screen**

---

# Slide 2 — The clinical problem and the computational gap

### Slide Title
**SCLC resists checkpoint blockade, and bulk differential expression cannot see why**

### Core Narrative / Speaker Notes
- The biology: T-cell dysfunction in the tumour microenvironment is a candidate mediator of SCLC's immune resistance, but its transcriptional basis is poorly characterised at single-cell resolution.
- The methodological gap — this is the sentence to land slowly: conventional differential expression compares *group means*. It tells you a gene differs between tumour and normal, but not what happens to a cell's state **if you change that gene**. It misses non-linear, context-dependent regulatory relationships.
- The move: a foundation model pretrained on millions of cells encodes those relationships inside each cell's full transcriptomic rank order, so you can ask a counterfactual question of a single cell instead of a group average.

### Visual Layout & Graphic Instructions
- Two-panel split with a vertical divider.
- **Left, "What we can do today":** two overlapping distribution curves (tumour vs normal) with a Δmean arrow between them; grey, deliberately flat and unexciting.
- **Right, "What we need to ask":** one cell icon with a question mark and the counterfactual in italic — *if this gene changed, would this cell look more like a different disease state?* Rendered in the teal *in silico* colour.
- Do not put a figure here. This slide is a concept contrast; keep it clean.

### On-Slide Text
- SCLC: high mutational burden, poor ICI response
- T-cell dysfunction is a candidate mediator — poorly characterised
- DE captures **mean shifts**, misses context-dependent regulation
- Need a **per-cell counterfactual**, not a group difference

---

# Slide 3 — What an *in silico* perturbation actually is

### Slide Title
**Four steps: rank a cell's genes, edit one, and measure where the cell moves**

### Core Narrative / Speaker Notes
- Walk the four steps in plain language, pausing on step 2 — the key mental model is that Geneformer reads a cell as an **ordered list of its genes**, most-expressed first. It is a rank order, not expression values.
- Deleting a gene drops it out of the list; overexpressing moves it to the top. Re-reading the modified cell gives a predicted change in cell state, averaged over thousands of real cells.
- Say the caveat out loud while the schematic is on screen: **nothing is edited in a laboratory here.** The edit is to the model's input; the result is the model's predicted change. A gene counts as a hit only when deleting and overexpressing move cells in **opposite** directions — that bidirectional requirement is your internal consistency check.

### Visual Layout & Graphic Instructions
- Reuse the poster's inline SVG explainer (`poster/poster_template.html`, the "What Is an *In Silico* Perturbation?" block) — it is vector, so it scales to 16:9 without loss.
- Four numbered stages left to right: ① One T cell · ② Rank its genes · ③ Edit one gene (DELETE drops out / OVEREXPRESS to top) · ④ Measure the move on a map of cell states.
- **Build this across 4 clicks.** This is the conceptual keystone for a non-modelling audience; do not reveal it all at once.
- Caveat line in a muted box pinned to the bottom edge, present from the first click.

### On-Slide Text
- A cell = an **ordered list** of its genes
- Delete → drop out · Overexpress → move to top
- Re-read the cell: which disease state does it now resemble?
- Hit = delete and overexpress move it in **opposite** directions

---

# Slide 4 — Cohort architecture and the donor guard

### Slide Title
**42 individuals, zero donor leakage — so accuracy is measured on people the model never saw**

### Core Narrative / Speaker Notes
- Provenance first: HTAN MSK collection via the CELLxGENE "T cells" dataset — 46,140 CD4+CD8 T cells × 24,540 features from 42 individuals.
- Explain the 45-vs-42 arithmetic before anyone asks, because it looks like an error: the split table has one row per **donor × disease**, so a donor contributing both tumour and normal tissue appears twice. The three-row difference is exactly the three paired tumour–normal donors (RU675, RU682, RU684).
- The credibility point: splits are donor-disjoint and audited — zero leakage across all 42 individuals, with a cross-disease guard for the paired donors. Every performance number downstream is therefore measured on people the model never saw. This is checked by assertion at build time, not by trusting a stored "PASS" string.

### Visual Layout & Graphic Instructions
- **Left 55%:** the cohort table (Disease × Train/Eval/Test/Donors) — SCLC 7,037 / 2,330 / 2,424 / 19 · LUAD 17,831 / 5,611 / 6,387 / 22 · Normal 2,334 / 1,620 / 566 / 4. Bold the SCLC row.
- **Right 45%:** a donor-split schematic — 42 person icons partitioned into three non-overlapping blocks, with the 3 paired donors highlighted in eosin rose and a lock icon labelled "cross-disease guard".
- Footer strip: `46,140 cells · 24,540 features · median 390–1,172 cells/donor · leakage audit PASS`.

### On-Slide Text
- HTAN MSK · CELLxGENE `6fde3ad9` — 42 individuals
- 45 donor×disease groups = 42 people + 3 paired tumour–normal
- Donor-disjoint splits, **zero leakage** (audited)
- Held-out test: 9,377 cells from 8 donors

---

# Slide 5 — Model performance

### Slide Title
**91.9% accuracy on held-out donors — with an SCLC→LUAD confusion worth noticing**

### Core Narrative / Speaker Notes
- Headline the metrics: 91.9% accuracy, macro F1 0.903 on 9,377 held-out cells from 8 donors. Donor-held-out performance separates SCLC, Normal and LUAD well enough to support the downstream perturbation screen.
- **Do not skip the diagonal.** Normal recall is 98.6% and LUAD 95.9%, but SCLC is 80.0% — and the 20% that is missed goes almost entirely to LUAD (484 cells). Flag it yourself before a reviewer does.
- Preview the tension: that SCLC→LUAD confusion is not just a performance detail. It reappears on slide 7 as a directional result the model produces, and it is the one thing in this work you have not reconciled. Planting it here makes slide 7 land as rigour rather than as a hole.

### Visual Layout & Graphic Instructions
- **Top strip:** four metric callout boxes — `91.9%` accuracy · `0.903` macro F1 · `9,377` held-out cells · `8` test donors.
- **Left 60%:** `sclc_confusion_matrix.png` (row-normalised, counts + percentages).
- **Right 40%:** fine-tuning setup as four compact key–value chips — Geneformer-V2-104M · 1 epoch, LR 5e-5, 6 frozen layers · NVIDIA GB10, 119 GiB unified · ~48 h GPU, ~5.88M cell-gene perturbations.
- Draw a thin eosin-rose ring around the SCLC→LUAD cell (484) and label it "revisited on slide 7".

### On-Slide Text
- **91.9%** accuracy · macro F1 **0.903** · donor-held-out
- Normal 98.6% · LUAD 95.9% · **SCLC 80.0%**
- SCLC misses go to LUAD — remember this
- ~48 h GPU · 1 epoch · 6 frozen layers

---

# Slide 6 — The applied checkpoint screen

### Slide Title
**Four edits replicate across every donor: TIM-3, TIGIT, CTLA-4, IL7R**

### Core Narrative / Speaker Notes
- Scope the screen: 50 checkpoint and T-cell engineering genes, tested on the SCLC → Normal transition. Of 123 concordant hits genome-wide, 43 are fully donor-consistent; within this targeted panel, four edits survive every filter.
- State the conjunction the four hits must satisfy, because it is stricter than "biggest effect": concordant in **both** perturbation arms, FDR < 0.05 in both, same sign in **all three** SCLC donors, and detected in at least 100 source cells. These rows are selected in code by that conjunction, so the table cannot drift from the criteria printed beside it.
- Read the biology plainly: deleting TIM-3, TIGIT, CTLA-4 or IL7R moves SCLC T cells toward a normal T-cell state; overexpressing moves them away. Then set up slide 7 — "but before you believe the ranking, one control changes how you read these numbers."

### Visual Layout & Graphic Instructions
- **Top strip:** four stat boxes — `50` genes screened · `4` replicated edits · `3/3` donors agree · `8` genes with no deletion result.
- **Centre:** the four-hit table, one row per gene, columns Gene / Delete→Normal / Overexpress / Detected / Donors — HAVCR2 (TIM-3) +0.0019 / −0.0060 / 202 · TIGIT +0.0018 / −0.0082 / 438 · CTLA-4 +0.0014 / −0.0014 / 282 · IL7R +0.0013 / −0.0054 / 1,131.
- **Right margin:** a vertical "filter funnel" — 50 genes → concordant both arms → FDR < 0.05 both arms → all 3 donors → ≥100 detected → **4**.
- Colour delete/overexpress columns in opposing hues so the sign flip is visible without reading numbers.

### On-Slide Text
- 50 checkpoint / engineering genes, SCLC → Normal
- Survives: both arms, FDR < 0.05, **all 3 donors**, ≥100 cells detected
- **TIM-3 · TIGIT · CTLA-4 · IL7R**
- Replicated and detection-adequate — *not* the largest effects

---

# Slide 7 — The detection confound and the unresolved axis

### Slide Title
**Sparse genes fake big effects (ρ = −0.60) — and one model result still contradicts the clinic**

### Core Narrative / Speaker Notes
- The control that reframes the screen: across the 50-gene panel, detection count and nominal effect size are strongly anti-correlated, Spearman ρ = **−0.60** (p = 2.7×10⁻⁵, n = 42). The fewer cells a gene is detected in, the larger its apparent shift. Effect magnitude alone is therefore **not** a ranking criterion — which is exactly why slide 6 ranked on replication and detection instead.
- Reinforce with the gaps: 8 of 50 genes have no deletion result at all — undetected in SCLC source cells, so deletion is undefined, not zero — and 18 fall below 100 detected cells.
- Then the honest part, and give it real time: on the exhaustion axis the model orders states **Normal < SCLC < LUAD**. TIGIT overexpression moves SCLC cells away from Normal (−0.0082) but *toward* LUAD (+0.0256). That runs opposite to the clinical picture of SCLC as a cold, ICI-resistant tumour. Present it as an **open question, not a reconciled result** — a tumour-level clinical phenotype and a T-cell-intrinsic model axis are not the same measurement, and you are not going to pretend they are.

### Visual Layout & Graphic Instructions
- **Left 55%:** `cart_overexpression.png` — detection (x) vs absolute deletion effect (y), with the sparse/high-effect corner circled and annotated `ρ = −0.60`.
- **Right 45%, visually separated by a rule:** a one-dimensional axis graphic, `Normal —— SCLC —— LUAD`, with the TIGIT overexpression arrow pointing from SCLC toward LUAD. Beneath it, two contrasting labels: "model axis (T-cell intrinsic)" vs "clinic (tumour level: SCLC cold, ICI-resistant)".
- Header the right half **OPEN QUESTION** in eosin rose. Do not soften this into a result — it is a discussion magnet and your strongest credibility signal.

### On-Slide Text
- Detection vs effect size: **ρ = −0.60** across the panel
- 8/50 genes: no deletion result · 18/50: <100 cells
- Effect size alone is **not** evidence
- **Open:** model says Normal < SCLC < LUAD; clinic disagrees

---

# Slide 8 — CAR-T relevance, and what this design can and cannot test

### Slide Title
**Two of the four hits are live CAR-T engineering targets — the tumour antigens are out of reach here**

### Core Narrative / Speaker Notes
- Separate the four hits by how interpretable they are, because they are not equivalent: **TIGIT** has the strongest external CAR-T support; **TIM-3** has direct knockout precedent in solid-tumour models; **CTLA-4** is a real signal with less CAR-T-specific evidence; **IL7R** is a persistence/fitness target, not an exhaustion checkpoint at all. Grouping them as "four checkpoints" would be wrong.
- The structural limit, stated as design not apology: this is a **T-cell atlas**, so it can only test T-cell-intrinsic edits. The CAR-T tumour antigens — DLL3, SEZ6, NCAM1, CD276, CEACAM5 — are not perturbable here, because they are not expressed by T cells. Likewise PD-L1 (CD274) is absent from the atlas, so this screen cannot speak to PD-L1 biology at all.
- Land the translational read: the actionable axis this work supports is **checkpoint knockout in the T-cell product** (TIGIT, HAVCR2) plus persistence engineering (IL7R) — which is precisely the CRISPR experiment on slide 12.

### Visual Layout & Graphic Instructions
- **Left 45%:** a 2×2 interpretability grid, x = evidence strength, y = mechanism type (exhaustion checkpoint vs persistence). Place TIGIT and TIM-3 upper right, CTLA-4 mid, IL7R in the persistence quadrant. Node colour = SCLC→Normal deletion shift.
- **Right 55%:** a "what this design can / cannot test" two-column card — green ticks for T-cell-intrinsic edits (checkpoints, persistence receptors); grey crosses for tumour-side antigens (DLL3, SEZ6, NCAM1, CD276, CEACAM5) and PD-L1, each with the one-word reason "not in T cells".
- Keep the crossed-out column visually calm, not alarming — this is scope, not failure.

### On-Slide Text
- **TIGIT, TIM-3** — strongest CAR-T / knockout precedent
- **IL7R** — persistence target, not a checkpoint
- Tumour antigens (DLL3, SEZ6…) **not testable** in a T-cell atlas
- PD-L1 absent from the atlas — no PD-L1 claim made

---

# Slide 9 — Spatial validation in independent tissue

### Slide Title
**In intact tumour tissue, dysfunction tracks T-cell abundance — and antigen presentation tracks it hardest**

### Core Narrative / Speaker Notes
- Set up the independence, which is the whole point: GSE263196 is a different cohort, a different assay and a different lab — five fresh-frozen SCLC sections, 15,632 of 15,774 in-tissue spots retained (99.1%), all 21 marker genes present in every specimen. No donor, cell or read is shared with the discovery atlas.
- Explain the measurement in one breath: at every spatially indexed spot you score two **disjoint** marker panels — T-cell abundance (CD3D/E/G, CD2, CD5, CD28, TRBC1/2, CD8A/B, CD4) and T-cell dysfunction (PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX, LAYN) — then correlate them within each specimen. The panels share no genes, so the correlation cannot be a self-correlation.
- Give both numbers, correctly attributed. The benchmark 7-gene dysfunction panel gives pooled ρ = **0.161 [0.146, 0.176]**, significant in 4 of 5 specimens. The 13-gene **antigen-presentation/MHC** programme is considerably stronger at ρ = **0.361**, rising to **0.384** when you control for sequencing depth, and sitting **7.95 σ** above an expression-matched random null. Depth control *raising* the estimate is the opposite of what a technical artifact does.
- Be straight about effect size: these are modest correlations with high between-specimen heterogeneity (I² ≈ 99%), consistent with real differences in tumour microenvironment composition rather than measurement noise.

### Visual Layout & Graphic Instructions
- **Top 60%:** `spatial_tissue_validation_panel.png` framed in the poster's virtual slide mount, with the frosted label strip (`GSE263196 · 10x Visium · fresh-frozen · 5 specimens · 15,632 spots · 55 µm spot / 100 µm pitch`). Row A = T-cell identity score, row B = dysfunction score.
- **Bottom 40%, left:** the per-specimen forest plot with the pooled diamond.
- **Bottom 40%, right:** a small programme-ranking bar chart — Cytotoxic effector 0.413 · **Antigen presentation 0.361** · Interferon 0.280 · Checkpoint/exhaustion 0.218 · *dysfunction benchmark 0.161*. Highlight antigen presentation; it is the programme in your take-home.
- Label the two ρ values distinctly on screen so they are never confused: "7-gene dysfunction benchmark" vs "13-gene antigen-presentation programme".

### On-Slide Text
- Independent cohort: 5 SCLC sections, **15,632 spots**, 99.1% retained
- Disjoint marker panels — correlation is not self-correlation
- Dysfunction benchmark **ρ = 0.161** [0.146, 0.176], 4/5 specimens
- Antigen presentation **ρ = 0.361 → 0.384** depth-controlled, **7.95σ**

---

# Slide 10 — Mechanism and literature context

### Slide Title
**The hits sit in one connected interaction neighbourhood — from prior knowledge, not inferred here**

### Core Narrative / Speaker Notes
- Show that the four replicated edits are not scattered: all 11 ICI/CAR-T candidates map onto a connected STRING neighbourhood, 54 filtered edges among 16 context genes, with TIM-3, TIGIT, CTLA-4 and IL7R sitting close to the cytotoxic and interferon machinery.
- State the epistemics clearly and without hedging: this map is a **literature overlay**, prior interaction evidence — it is **not** a gene-to-gene network inferred from these cells. It provides context for the hits; it is not independent evidence for them.
- Point out the informative absences, which show the overlay is doing real work rather than decoration: TOX and LAYN carry no STRING edge above threshold, and PD-L1 is not in the T-cell atlas at all. TIM-3 is the strongest individual hit but has no supported non-text-mining edge to another checkpoint here — an explicit asymmetry, not a plotting omission.

### Visual Layout & Graphic Instructions
- **Left 60%:** panel **c** of `ici_cart_perturbation_network.png` — the STRING map with nodes coloured by SCLC→Normal deletion shift and the four replicated hits bold.
- **Right 40%:** a four-row evidence card, one per hit, each with a strength meter — TIGIT (CAR-T support: strong) · TIM-3 (solid-tumour KO precedent: strong) · CTLA-4 (real signal, less CAR-T-specific) · IL7R (persistence, not exhaustion).
- Pin a persistent caption bar across the bottom: **"Prior interaction evidence — not a network inferred from these cells."**
- Crop tightly. At full width the node labels in `perturbation_networks.png` are ~5 pt; use the cropped panel **c**, not the three-panel composite.

### On-Slide Text
- 11 candidates · 54 STRING edges · 16 context genes
- Hits cluster with cytotoxic and interferon machinery
- **Literature overlay — not inferred from these cells**
- TOX, LAYN: no edge above threshold; PD-L1 absent

---

# Slide 11 — Scope limits

### Slide Title
**What this study does not claim**

### Core Narrative / Speaker Notes
- Deliver these as deliberate scope, briskly and without apology — a clean limits slide in a talk pre-empts the hostile version of every Q&A question, and this audience will respect it.
- The one that matters most: perturbation shifts are **model-derived interventions on a rank-value encoding, not experimental knockouts**. Everything upstream is a prioritised hypothesis. Say this sentence verbatim.
- The measurement caveats, quickly: detection is source-state specific, so absence of a result is not evidence of absence; the spatial work uses marker-score proxies rather than formal deconvolution; there are no Normal or LUAD sections in the validation cohort, so **SCLC-specificity of the spatial association is untested**; and the design is correlational, in archival tissue.

### Visual Layout & Graphic Instructions
- Single centred column of four short lines, generous leading, muted grey-on-white. No figure, no icons, no decoration.
- Set the title in the same weight as the results slides — this reads as confidence, not retreat.
- Optional: a faint watermark of the pipeline strip at 5% opacity to keep visual continuity.
- Hold this slide briefly. It is a spoken slide; the audience should be listening, not reading ahead.

### On-Slide Text
- Model-derived predictions, **not experimental knockouts**
- Marker-score proxy, not deconvolution
- No Normal/LUAD sections — SCLC-specificity untested
- Correlational design, archival tissue

---

# Slide 12 — Summary and what happens next

### Slide Title
**A tissue-validated hypothesis, and the CRISPR experiment that would test it**

### Core Narrative / Speaker Notes
- Close the loop in three beats: a donor-clean fine-tuned model (91.9%, F1 0.903) → four replicated, detection-adequate candidate edits (TIM-3, TIGIT, CTLA-4, IL7R) → an antigen-presentation programme that tracks T-cell abundance in independent tumour tissue (ρ 0.361, 7.95σ).
- Name the falsifying experiment concretely, because it converts the talk from a computational exercise into a research programme: **CRISPR knockout of TIGIT and HAVCR2 in an SCLC co-culture**, and Perturb-seq in primary T cells read directly against these *in silico* predictions. That is what turns predicted shifts into measured ones.
- Leave the open question open on the final slide — the Normal < SCLC < LUAD axis. Inviting the room to argue about it is a better close than a tidy conclusion, and it is the honest state of the work.

### Visual Layout & Graphic Instructions
- **Top third:** the same three-node take-home flow as slide 1, now with each node ticked — visually closing the loop for anyone who joined late.
- **Middle:** two side-by-side next-step cards. **Experimental** — CRISPR KO of TIGIT/HAVCR2 in SCLC co-culture; Perturb-seq in primary T cells vs *in silico*. **Computational** — cell2location deconvolution; Normal/LUAD specificity controls; scGPT cross-model replication.
- **Bottom right:** large QR to `github.com/Kays3/geneformer-lung-tcell` (code, result tables, notebooks, poster), plus contact `snakaoka@sci.hokudai.ac.jp`.
- **Bottom left:** a single boxed line — *Open: why does the model order Normal < SCLC < LUAD?* — as the deliberate Q&A hook.

### On-Slide Text
- Donor-clean model → 4 replicated edits → tissue-validated programme
- Next: **CRISPR KO of TIGIT / HAVCR2** in SCLC co-culture
- Then: Perturb-seq vs *in silico* predictions
- Code + data: QR

---

## Delivery notes

- **The three sentences to get right, verbatim:** "Nothing is edited in a laboratory here" (slide 3) · "Effect size alone is not evidence" (slide 7) · "Model-derived predictions, not experimental knockouts" (slide 11). Each pre-empts a specific objection.
- **If you must cut to 10 slides:** fold slide 11 into slide 12 as a half-panel, and merge slide 10 into slide 8 (both are interpretability/context). Never cut slide 7 — the confound control is what makes slide 6 credible.
- **Likely Q&A:** (1) Why is SCLC recall only 80%? → 484 cells go to LUAD; same axis as the open question on slide 7. (2) Is ρ = 0.16 too small to matter? → modest but 3.88σ above a matched null, and the MHC programme reaches 7.95σ; heterogeneity is between-specimen, I² ≈ 99%. (3) Could the spatial correlation be a depth artifact? → depth control *raises* partial ρ from 0.361 to 0.384. (4) Why not PD-L1? → absent from the T-cell atlas; no claim made.
- **Figure sources:** all referenced PNGs live in `sclc_validation/perturbation_workflow/figures/`, `sclc_validation/spatial_validation/figures/` and `sclc_validation/checkpoint_cart_perturbation/figures/`. The slide-3 explainer is inline SVG in `poster/poster_template.html`.
