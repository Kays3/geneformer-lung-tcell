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

## Exporting a PDF

Open the generated HTML in Chrome:

Print → Paper **A0** → Margins **None** → Scale **100%** →
**Background graphics ON** → Save as PDF.

Background graphics is not optional: with it off, every coloured panel,
header and slide mount prints white.

The organizers supply the 20 × 20 cm presentation-number panel; the poster
reserves matching space at the top left showing the abstract number.

## If the content does not fit the sheet

The poster is one fixed A0 page. Rather than clipping silently, an overrun is
made visible: a **red dashed rule** marks the A0 bottom edge on screen.

- Content sits entirely above the rule → it fits, print as is.
- Anything appears below the rule → lower `--fit` in `poster_template.html`
  (`:root` block, near the top) and regenerate.
- Noticeable blank space above the rule → raise `--fit` slightly.

`--fit` scales **every** font size on the poster — all 46 declarations,
including the title, the P25 box and the large stat figures:

| `--fit` | title | body | smallest label |
|---|---|---|---|
| `1.00` | 40 pt | 17.5 pt | 10.5 pt |
| `0.94` | 38 pt | 16.4 pt | 9.9 pt |
| `0.88` | 35 pt | 15.4 pt | 9.2 pt ← practical floor |

Do not go below about `0.88`: the smallest labels fall under 9 pt and stop
being readable at normal poster viewing distance. If it still overruns at
`0.88`, cut content rather than shrinking further.

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
