# JSDP poster — A0 portrait

Abstract P25 — Dauyey K & Nakaoka S, Laboratory of Mathematical Biology,
Faculty of Advanced Life Science, Hokkaido University.

Sheet size **A0 portrait, 841 × 1189 mm**, styled around the digital
pathology narrative (single-cell perturbation → spatial validation in tissue).

## Files

| File | Purpose |
|---|---|
| `generate_poster.ipynb` | Run this. Reads the result files, fills the template, writes the poster. |
| `poster_template.html` | A0 layout and wording. Edit to move panels or change prose. |
| `poster_template_90x210.html` | Previous 90 × 210 cm layout, kept for the full JSDP board. |
| `poster_<draft>.html` | Generated output. Not committed by default. |
| `../tools/make_hero_figure.py` | Builds Fig. 4, the full-width key-message band. |
| `../tools/make_snapshots.py` | Writes the tracked PNG preview of the final PDF. |
| `../sclc_validation/checkpoint_cart_perturbation/scripts/make_network_figure.py` | Builds Fig. 6. |

## Sheet size note

A0 (841 × 1189 mm) is **53% of the JSDP board area** (900 × 2100 mm). A0 is the
standard print size and is what this template targets. To fill the whole board
instead, swap in `poster_template_90x210.html`.

## Digital pathology theme

- **Palette** derived from H&E staining — hematoxylin violet (`--hema`) for
  nuclear/spatial content, eosin rose (`--eosin`) for tissue-level results.
  Computational sections keep the neutral navy/teal scheme, so the eye can
  separate *in silico* from *in tissue* at a glance.
- **Virtual slide mount** frames the Visium figure like a glass slide on a
  stage, with a frosted label strip carrying accession, platform and specimen
  count.
- **Slide tray** renders one card per specimen; the bar length is proportional
  to that specimen's rho and is coloured only when it individually reached
  p < 0.001, so a non-significant specimen cannot read as a strong result.
- **Pipeline strip** shows block → cryosection → capture → spot transcriptome →
  score map.

Panel accents `hema` and `eosin` are available to `EXTRA_PANELS` for new
pathology content.

## Final — key-message band, light-ink theme

The current final poster (`DRAFT_LABEL = "Final"`, built after Draft 12) differs
from every earlier draft in three ways that change how you edit it:

**The body is two band grids, not three columns.** `.body` is a block containing
`.band` (three columns) → `.hero` (full width) → `.band` (three columns). Two
consequences. First, **only the tallest column in each band sets that band's
height**, so trimming a short column frees nothing — the admin panels dropped in
the final draft saved 0 mm on their own, and it took rebalancing both bands to
recover the space. Second, the two bands each align to a shared baseline, so a
column much shorter than its neighbours shows as dead space.

**The one-page check is not the page count.** WeasyPrint clips band overflow
instead of paginating it: the export reported one A0 page while the bottom band
was cut off mid-panel. Cell 6 now re-renders onto a 4000 mm page, sums the
laid-out band heights and asserts against 1189 mm. Trust that number, and look
at the render.

**The theme is light-ink.** Header, footer, section bars, table headers and the
slide mount are tinted grounds with dark text rather than solid navy or
hematoxylin fills — an A0 sheet of solid dark areas is a lot of toner and dries
unevenly. Colour is carried by rules, left borders and text.

Earlier drafts' notes below are kept as history; where they disagree with this
section, this section is current.

## History — superseded drafts

Everything below records how the poster got here. It is kept because the
reasoning still explains why parts of the layout are the way they are, but
**where it disagrees with the section above, it is out of date** — the column
widths, `--fit` values and figure inventory quoted in these notes are all
superseded.

### What Draft 4 added

Two panels, both aimed at readers who are not foundation-model people:

**"What Is an *In Silico* Perturbation?"** (left column, above the cohort) — a
plain-language explainer, because every number in sections B and C is a
perturbation shift and most JSDP attendees are clinicians and wet-lab
biologists. It lives directly in `poster_template.html`, not in
`EXTRA_PANELS`, because it is static prose and needs to sit near the top of
the column rather than appended to the end. The four-step schematic is
**inline SVG**, so it stays sharp at A0 and needs no figure file. Two things
to know before editing it:

- Its text is sized in viewBox units, so it scales with the column width
  (fixed at 248 mm) and **not** with `--fit`. At the current geometry the
  smallest labels land near 10 pt; check them against the 9 pt floor if the
  column width ever changes.
- The SVG's rendered height follows its viewBox aspect ratio, so slack
  viewBox units become dead millimetres in the tallest column on the sheet.
  The viewBox height is trimmed to just past the lowest element for that
  reason. Also note elements are painted in document order — the stage-4 box
  is drawn after the stage-3 captions and will cover anything that runs
  under it.

The applied checkpoint/CAR-T read-out now spans the middle and right columns,
using the result package in `sclc_validation/checkpoint_cart_perturbation/`.
The middle column carries the detection-aware four-hit table; the right column
carries literature context, the STRING overlay, and the unresolved checkpoint
axis.

The four rows (TIM-3, TIGIT, CTLA-4, IL7R) are selected in Cell 1 by the
conjunction the panel claims — concordant in both arms, `tier == "all donors
agree"`, and not `low_detection_lt100` — rather than typed in, so the table
cannot drift from the criteria printed beside it.

The updated result package makes detection a first-class part of the claim:
counts are joined on `(gene, comparison)`, eight genes have no deletion result,
and detection count is strongly anti-correlated with nominal effect size
(Spearman ρ = −0.60). Draft 5 therefore presents the four hits as replicated,
detection-adequate candidates rather than as the largest effects in the panel.

The panel also carries the SCLC↔LUAD sign flip as an **open question, not a
result**, per that analysis's own interpretation notes: the model orders these
states Normal < SCLC < LUAD on the exhaustion axis, which runs opposite to the
clinical picture. Do not let a later edit smooth this into either narrative —
the two are not the same measurement.

### What Draft 5 changed

Draft 5 is the conference-facing revision of Draft 4. It shortens the title,
adds a single take-home strip to the header, and renames the section bars so a
reader can follow the story as question/cohort → model/screen → tissue
validation → conclusions. It replaces the former genome-scale/ambient prose
with detection-aware checkpoint results, scope limits, literature context and
new figures. The template uses `--fit: 0.88` so the new panels still
strip still exports as one A0 page; the smallest labels remain above the
documented readability floor.

### Draft 6 Minimal Picture-led Version

`poster_template.html` is now the picture-led source used by
`generate_poster.ipynb`. It uses relative paths to the current classifier,
detection-confound, spatial-validation and checkpoint-network figures, and
keeps only the four-hit table. The rendered `poster_draft_6.html` and
`poster_draft_6.pdf` are generated deliverables.

### Draft 7 Reduced Draft 5

Draft 7 restores the Draft 5 visual system through `generate_poster.ipynb`,
removes the classifier, paired-donor, specimen-QC and validation tables, and
keeps the four-hit table plus the spatial, detection-confound, forest and
checkpoint-network figures.

### Superseded: three-column rebalance (Draft 10)

Draft 10 was Draft 9's content in a rebalanced grid. It is worth knowing why the
grid changed, because the obvious edit reintroduces the problem.

Through Draft 8 the middle column carried **both** section B and section C, and
the extra panels that padded the outer columns were dropped one by one for
space. Draft 9 finished that trend and put its two new denoising figures in the
middle column as well. The result measured **left 59.5% / middle 95.0% /
right 59.5%** of sheet height — the middle column nearly overran A0 while the
outer two sat 40% empty, and `--fit` had to sit at the 0.88 readability floor
to keep one page.

The fix moves **section C to the top of the right column**, so the columns read
A | B | C over D. That frees the middle column for the two denoising figures at
`WIDTH_MM["middle_full"] = 260` and leaves room in the left column for the
restored detection-confound panel. Measured after the move: **left 88.3% /
middle 88.4% / right 94.0%**, one A0 page, with `--fit` back up to 0.96.

Two things to keep in mind before editing this layout:

- `WIDTH_MM["spatial"]` is derived from `_RIGHT_COL`, which is now the column
  section C actually sits in. Before the move it was sized for a column it was
  not in; `object-fit:contain` letterboxed the image rather than distorting it,
  so the figure rendered at the same 233 mm either way and the bug was
  invisible. If section C ever moves again, re-derive this key from its new
  column or the figure silently shrinks inside its slide mount.
- Balance is measurable, so measure it rather than eyeballing the render. The
  per-column fill percentages quoted here come from finding the lowest
  non-background pixel row within each column's x-range of the snapshot, with
  the header and footer bands excluded. A one-page assertion pass does **not**
  imply a balanced sheet — Draft 9 passed it.

## Figures

Eight figures, numbered in reading order: top band left→right, then the
key-message band, then the bottom band left→right. **Keep them consecutive** —
an earlier draft dropped the forest plot and left the sheet numbering 1, 2, 4–9
with no Fig. 3, which reads as a mistake to anyone at the poster.

| Figure | Panel / column | Source data | Generator |
|---|---|---|---|
| 1 Confusion matrix | Classifier, middle top | `test_confusion_matrix.csv` | `tools/make_result_figures.py` — see "Confusion-matrix fix" |
| 2 Qualified genes heatmap | Classifier, middle top | `primary_test_perturbation/tables/` | `primary_test_perturbation/scripts/build_denoised_notebook.py` |
| 3 Spatial score maps | Section C, right top | `spatial_validation/results/*` | pre-existing Visium tissue panel |
| 4 Key-message band | Full-width band | several committed tables | **`tools/make_hero_figure.py`** |
| 5 Detection/effect confound | Checkpoint screen, left bottom | `checkpoint_cart_perturbation/tables/cart_overexpression_vs_deletion.csv` | no committed generator |
| 6 Perturbation network | Network context, middle bottom | `checkpoint_cart_perturbation/tables/{network_node_perturbation,string_network_edges}.csv` | **`checkpoint_cart_perturbation/scripts/make_network_figure.py`** |
| 7 ICI / CAR-T candidates | ICI/CAR-T, middle bottom | `checkpoint_cart_perturbation/tables/` | no committed generator |
| 8 Denoised bidirectional | Denoised programs, right bottom | `primary_test_perturbation/tables/immune_cancer_candidates.csv` | `build_denoised_notebook.py` (`_plot_bidirectional`) |

Figure numbers appear both as captions and as cross-references in panel prose,
in `poster_template.html` **and** in `EXTRA_PANELS` in Cell 1. Renumber both
together or the references drift from the captions.

Two figures still have no committed generator (5 and 7). Fig. 6 had none either
until it needed resizing and could not be rebuilt — that is what the
`make_network_figure.py` docstring records. If either of these two ever needs to
change, write the generator rather than editing the PNG.

**Fig. 4 carries a constraint.** Its panel B paints per-gene effects on cell
state onto prior-knowledge STRING edges; the screen contains no gene-to-gene
measurements. The edges are deliberately pale and unweighted and the caveat is
inside the axes. `make_hero_figure.py` has a do-not-remove note explaining this;
restyling those edges would make the figure claim something never measured.

## Cohort panels

Two panels describe the data provenance, both computed from source tables at
run time rather than typed in:

**Discovery cohort (left column)** — HTAN MSK, CELLxGENE "T cells" dataset
`6fde3ad9`. Reads `htan_donor_split_assignment.csv`,
`htan_tcell_summary_by_disease.csv` and `cellxgene_collection_inventory.csv`.

The donor arithmetic is worth understanding before editing it: the split
table has one row per **donor × disease**, so a donor contributing both
tumour and normal tissue appears twice. That is why 45 donor×disease groups
map to 42 individuals — the 3-row difference is exactly the paired
tumour-normal donors (RU675, RU682, RU684), which the panel names explicitly.
Cell 4 asserts this arithmetic closes, and separately asserts no individual
spans two splits, rather than trusting the stored `donor_leakage_check: PASS`
string. If either assertion fires, the poster's "measured on individuals the
model never saw" claim is no longer supported — fix the split, don't relax
the assertion.

**Validation cohort (right column)** — GSE263196, from
`gse263196_spatial_file_audit.csv` joined to the correlation results. Reports
per-specimen in-tissue vs. analysed spot counts (15,632 of 15,774 retained,
99.1%), genes detected, and marker-panel coverage (21/21 in every specimen).
Cell 4 asserts the audit's spot total matches the results table's, so a
mismatch between the QC file and the analysis surfaces as an error rather
than two different numbers on the same poster.

Adding these panels cost roughly 107 mm of column height and pulled the
one-page `--fit` ceiling from 1.60 down to 1.12. See the fit table below.

### Confusion-matrix fix

`FIG_CONFUSION` previously pointed at
`current_workflow/visuals/final_tcell_confusion_matrix.png` — a leftover
from the pre-SCLC-pivot workflow showing **LUAD/LUSC/Normal** counts. It sat
directly next to a table built from `test_confusion_matrix.csv`, which
correctly shows **SCLC/Normal/LUAD** — the image and the table it was
captioned to describe disagreed. `sclc_confusion_matrix.png` was generated
directly from `test_confusion_matrix.csv` so the figure and its caption
agree; regenerate it with `tools/make_result_figures.py` if the underlying
metrics change, not by pointing back at the old file.

## Making a new draft

1. Open `generate_poster.ipynb`.
2. Edit **Cell 1** only:
   - `DRAFT_LABEL` — names the output file, prints on the poster.
   - `FIG_*` — repoint at regenerated figures.
   - `EXTRA_PANELS` — add a dict per new analysis (see below).
   - `CONCLUSIONS` / `LIMITATIONS` / `FUTURE` — edit the prose.
3. Run all cells.
4. Output: `poster/poster_<draft_label>.html`.

Numbers on the poster (accuracy, macro F1, concordant hit counts, spatial rho,
cohort table, confusion matrix, per-sample spatial table) are read from the
result files at run time. They are never typed by hand, so a rerun after new
analysis picks up the new values automatically.

## Adding a panel for a new analysis

Append to `EXTRA_PANELS` in Cell 1:

```python
EXTRA_PANELS = [
    {
        "column": "middle",              # left | middle | right
        "title": "Genome-scale screen",
        "accent": "sky",                 # "" | sky | mauve | amber | green | mid
        "text": "All-gene screen across held-out SCLC cells.",
        "figure": REPO / "path/to/figure.png",     # or None
        "caption": "Fig. 4 - Top genes by goal-state shift.",
        "table": REPO / "path/to/table.csv",       # or None
        "table_rows": 10,
    },
]
```

A figure or table path that does not exist yet renders as a dashed "Pending"
box, so the poster still builds while an analysis is in flight.

There is also an `"html"` key, rendered verbatim just after `"text"`. Use it
when a panel needs a curated table — chosen columns, formatted numbers, units
in the header — rather than the whole-CSV dump `"table"` produces. The
checkpoint/CAR-T panel uses it; build the row markup in Cell 1 from the
dataframe so the numbers still come from the result files.

Every column a figure panel can sit in needs its own `WIDTH_MM` entry. If the
key is missing, `render_panel` passes `width_mm=None`, the image gets no inline
height, and WeasyPrint silently re-splits the poster across several
mostly-blank pages instead of erroring. Cell 3 defines `left_full` and
`right_full`, and the cell after it adds `middle_full`; add a matching entry
before putting a figure in any new column.

## Exporting a PDF

**Preferred: Cell 6 in the notebook.** Renders the HTML with WeasyPrint,
which reads the template's `@page { size: A0 portrait }` rule directly, so
the output is exactly 841 × 1189 mm with no print-dialog steps. Needs
WeasyPrint's system libraries (Pango, Cairo, GDK-Pixbuf) — a bare
`pip install weasyprint` is usually not enough:

```
conda install -c conda-forge weasyprint pango cairo gdk-pixbuf
```

The cell asserts the export landed on exactly one A0 page and fails loudly
if not — see "Why every figure has an inline height" below.

**Fallback: browser print.** Open the generated HTML in Chrome:
Print → Paper **A0** → Margins **None** → Scale **100%** →
**Background graphics ON** → Save as PDF. Background graphics is not
optional: with it off, every coloured panel, header and slide mount prints
white.

The organizers supply the 20 × 20 cm presentation-number panel; the poster
reserves matching space at the top left showing the abstract number.

### Why every figure has an inline height

WeasyPrint 69 mis-paginates this template's `<img>` tags when they have
`width:100%` and no matching height — it silently splits one page of content
into several mostly-blank pages instead of raising an error. This has
nothing to do with content actually overflowing A0; a browser renders the
identical markup on one page.

The fix lives in `img_tag(..., width_mm=...)` in Cell 3: it reads each
figure's real pixel dimensions and writes an explicit `height` in
millimetres, computed from the known column width in
`poster_template.html`'s CSS grid. If you resize a column or add a new
figure panel, add its width to `WIDTH_MM` in Cell 3 so the image gets a
matching inline height — an image without one will silently reintroduce the
multi-page bug, which is why Cell 6 asserts the page count rather than
trusting the render.

## If the content does not fit the sheet

The poster is one fixed A0 page. Rather than clipping silently, an overrun is
made visible: a **red dashed rule** marks the A0 bottom edge on screen.

- Content sits entirely above the rule → it fits, print as is.
- Anything appears below the rule → lower `--fit` in `poster_template.html`
  (`:root` block, near the top) and regenerate.
- Noticeable blank space above the rule → raise `--fit` slightly.

`--fit` scales **every** font size on the poster — all 46 declarations,
including the title, the P25 box and the large stat figures:

| `--fit` | title | body | smallest label | fits 1 A0 page? |
|---|---|---|---|---|
| `0.88` | 35 pt | 15.4 pt | 9.2 pt ← practical floor | yes |
| `0.90` | 36 pt | 15.8 pt | 9.5 pt ← **current default (Final)** | yes — 1150 mm of 1189 mm, +39 mm margin |
| `0.92` | 37 pt | 16.1 pt | 9.7 pt | yes, but only ~0.5 mm margin — too close |
| `0.96` | 38 pt | 16.8 pt | 10.1 pt | no — overruns with the key-message band |
| `0.99` | 40 pt | 17.3 pt | 10.4 pt | yes — measured with Draft 4 content |
| `1.00` | 40 pt | 17.5 pt | 10.5 pt | yes — measured with Draft 4 content |
| `1.12` | 45 pt | 19.6 pt | 11.8 pt | Draft 3 ceiling; not re-measured for Draft 4 |
| `1.15` | 46 pt | 20.1 pt | 12.1 pt | no — spilled to a 2nd page with Draft 3 content |

Draft 4 went *up* from 0.96 to 1.00, which looks backwards for a draft that
added content. It added ~291 mm of column height (the explainer and the
checkpoint panel) and paid for it by dropping Figs 4 and 5, which are worth
~331 mm together — so the sheet came out with slack, and the type went back up
to use it. The Draft 4 ceiling above 1.00 was never measured.

The Final layout sits at `0.96` with its tallest column (the right, carrying
C over D) filling 94.0% of the sheet. That leaves roughly 6% of headroom, so
the ceiling here is tight — a step to `1.00` is about a 4% type increase and
would very likely land past the A0 edge. Re-measure before raising it, and
prefer moving a panel between columns over shrinking type if content grows.

Do not go below about `0.88`: the smallest labels fall under 9 pt and stop
being readable at normal poster viewing distance.

The ceiling is **content-dependent**, not a fixed limit, and it moves a lot:
it was `1.60` before the two cohort panels were added and is `1.12` now. It
is found by rendering at each value and checking the page count Cell 6
reports. Any edit to panel text, table rows, or figures changes how much
room is left, so re-measure rather than assuming a value here still holds.
If it overruns at your target `--fit`, cut content rather than shrinking
below the 0.88 floor.

The red rule is hidden in print output.

## Large files and what gets committed

Generated posters are big — a rendered HTML carries ~10 MB of embedded figures,
and the exported PDF is similar. Committing every revision would bloat the
repository permanently, so they are git-ignored:

| Path | Tracked? |
|---|---|
| `poster/poster_template*.html`, `generate_poster.ipynb`, `README.md` | yes — source |
| `poster/poster_<draft>.html` | no — regenerate from the notebook |
| `output/pdf/*.pdf`, `tmp/` | no — exported deliverables |
| `snapshots/*.png` | **yes** — small PNG previews |

So the repository still shows what was produced, render a preview after
exporting a PDF:

```
python tools/make_snapshots.py           # write previews into snapshots/
python tools/make_snapshots.py --check   # report stale previews, exit 1 if any
```

Each preview is ~0.7 MB with the long edge at 2000 px, enough to read the panel
layout and headings.

One deliberate exception: `spatial_tissue_validation_panel.png` is 7.5 MB but
stays tracked, because `generate_poster.ipynb` embeds it as an input — dropping
it would break poster generation from a fresh clone. Shrinking it is not an
option either; getting it under 5 MB requires palette quantisation, which
introduced errors up to 61/255 on the expression colour map.

## Requirements

```
pip install pandas qrcode Pillow pypdfium2
```

On macOS, WeasyPrint also needs the Homebrew libraries on the loader path:

```
brew install pango gdk-pixbuf libffi
DYLD_FALLBACK_LIBRARY_PATH=/opt/homebrew/lib python -c 'import weasyprint'
```

For the WeasyPrint PDF export (Cell 6), a conda environment with the system
libraries is more reliable than pip alone:

```
conda create -n pdfrender -c conda-forge python=3.11 weasyprint pango cairo gdk-pixbuf
conda activate pdfrender
pip install pandas qrcode Pillow pypdfium2
```
