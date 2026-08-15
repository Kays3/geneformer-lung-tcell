#!/usr/bin/env python3
"""Build the in-silico perturbation explainer (poster Fig. 1).

Four steps, in plain language, for readers who are not foundation-model people:
rank a cell's genes, edit one of them, re-read the cell, measure how far it
moved. Every number on the rest of the sheet is a perturbation shift, so this
is the figure that makes the others readable.

    python tools/make_explainer_figure.py

This figure carries no data. It is a schematic, and the gene names in it are
illustrative - deliberately chosen from the panel the poster actually screens
(CD3D, IL7R, TIGIT, GZMB, TOX) so the example is not invented, but the bar
lengths are drawn for legibility and are not expression values. Nothing here
should be read as a measurement; the caveat under the panel says so on the
sheet as well.

It replaces an inline SVG that lived in poster_template.html. The SVG rendered
well but was the one figure on the poster with no generator, and its text was
sized in viewBox units, so it silently rescaled whenever the column width
changed. As a PNG it goes through the same img_tag path and the same explicit
inline height as every other figure.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch

REPO = Path(__file__).resolve().parent.parent
OUT = REPO / "sclc_validation" / "perturbation_workflow" / "figures"
OUT.mkdir(parents=True, exist_ok=True)

NAVY = "#10243c"
SKY = "#17608f"
TEAL = "#0b6f6a"
HEMA = "#3d2f6b"
EOSIN = "#b8465e"
INK = "#16202c"
INK2 = "#4a5768"
RULE = "#d5dee7"
WASH = "#e9f2fa"

# Illustrative only - see the module docstring.
GENES = [("CD3D", 1.00), ("IL7R", 0.78), ("TIGIT", 0.56), ("GZMB", 0.40), ("TOX", 0.26)]
EDITED = "TIGIT"

# Sized for the poster's left column (242 mm inner). Text is in points and does
# not scale with the canvas, so a smaller canvas renders larger on the sheet:
# at 9 in wide this scales up ~1.06x, putting the 11 pt labels near 11.6 pt.
FIGSIZE = (9.0, 2.65)

X0, X1 = 0.0, 100.0
Y0, Y1 = 0.0, 30.0
STAGES = [2.0, 27.5, 53.0, 78.5]      # left edge of each stage
WIDTH = 19.5


def _badge(ax, x, y, number, title, colour):
    ax.add_patch(Circle((x + 1.0, y), 1.15, facecolor=colour, edgecolor="none",
                        zorder=4))
    ax.text(x + 1.0, y, str(number), ha="center", va="center", color="white",
            fontsize=8.5, fontweight="bold", zorder=5)
    ax.text(x + 3.0, y, title, ha="left", va="center", color=INK,
            fontsize=10.5, fontweight="bold")


def _arrow(ax, x, y):
    ax.add_patch(FancyArrowPatch((x, y), (x + 3.6, y), arrowstyle="-|>",
                                 mutation_scale=11, linewidth=1.5, color="#93a3b5"))


def _gene_rows(ax, x, y_top, *, edited=None):
    """Five ranked gene rows. `edited` greys and strikes one out."""
    edited_index = next((i for i, (g, _) in enumerate(GENES) if g == edited), None)
    for index, (gene, weight) in enumerate(GENES):
        y = y_top - index * 3.15
        struck = gene == edited
        ax.add_patch(FancyBboxPatch((x, y - 1.25), WIDTH, 2.5,
                                    boxstyle="round,pad=0.12",
                                    facecolor="#eef1f5" if struck else "#f6f8fb",
                                    edgecolor="#c9d3de" if struck else RULE,
                                    linestyle="--" if struck else "-",
                                    linewidth=0.9))
        ax.text(x + 0.9, y, gene, ha="left", va="center",
                fontsize=9.5, fontweight="bold" if gene == EDITED else "normal",
                color="#93a3b5" if struck else (EOSIN if gene == EDITED else INK))
        if struck:
            ax.plot([x + 0.7, x + 6.6], [y, y], color=EOSIN, linewidth=1.4, zorder=5)
            ax.text(x + WIDTH - 1.4, y, "✕", ha="right", va="center",
                    fontsize=10, color=EOSIN, fontweight="bold")
        elif edited is None:
            bar = 8.5 * weight
            ax.add_patch(FancyBboxPatch((x + 9.5, y - 0.62), bar, 1.25,
                                        boxstyle="round,pad=0.05",
                                        facecolor=EOSIN if gene == EDITED else SKY,
                                        edgecolor="none"))
        elif edited_index is not None and index > edited_index:
            ax.add_patch(FancyArrowPatch((x + WIDTH - 3.0, y - 0.9),
                                         (x + WIDTH - 3.0, y + 0.9),
                                         arrowstyle="-|>", mutation_scale=8,
                                         linewidth=1.2, color=TEAL))


def main() -> Path:
    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.set_xlim(X0, X1)
    ax.set_ylim(Y0, Y1)
    ax.set_axis_off()

    # ---- 1. one cell -----------------------------------------------------
    x = STAGES[0]
    _badge(ax, x, 28.0, 1, "One T cell", SKY)
    cx, cy = x + WIDTH / 2 - 1.0, 13.0
    ax.add_patch(Ellipse((cx, cy), 15.5, 15.0, facecolor=WASH, edgecolor=SKY,
                         linewidth=1.6))
    ax.add_patch(Ellipse((cx, cy), 7.0, 6.4, facecolor="#cfe1f1", edgecolor=SKY,
                         linewidth=1.1))
    for dx, dy in [(-5.2, 4.6), (5.0, 4.2), (-5.4, -4.4), (5.2, -4.6), (0, 6.2)]:
        ax.add_patch(Circle((cx + dx, cy + dy), 0.55, facecolor=SKY, edgecolor="none"))
    ax.text(cx, cy, "nucleus", ha="center", va="center", fontsize=7.5, color=NAVY)

    # ---- 2. rank the genes ----------------------------------------------
    _arrow(ax, STAGES[0] + WIDTH + 1.0, 13.0)
    x = STAGES[1]
    _badge(ax, x, 28.0, 2, "Rank its genes", SKY)
    ax.text(x, 24.4, "most expressed", fontsize=7.5, color=INK2)
    _gene_rows(ax, x, 22.0)
    ax.text(x, 4.4, "⋮", fontsize=9, color="#93a3b5")
    ax.text(x + 5.0, 4.4, "least expressed", fontsize=7.5, color=INK2)

    # ---- 3. edit one gene ------------------------------------------------
    _arrow(ax, STAGES[1] + WIDTH + 1.0, 13.0)
    x = STAGES[2]
    _badge(ax, x, 28.0, 3, "Edit one gene", EOSIN)
    ax.text(x, 24.4, "DELETE — drop it out", fontsize=8, fontweight="bold", color=EOSIN)
    _gene_rows(ax, x, 22.0, edited=EDITED)
    ax.text(x, 4.6, "the rest move up a rank", fontsize=7.5, color=INK2)
    ax.text(x, 1.9, "OVEREXPRESS", fontsize=8, fontweight="bold", color=TEAL)
    ax.text(x + 11.5, 1.9, "— move it to the top", fontsize=7.5, color=INK2)

    # ---- 4. measure the move --------------------------------------------
    _arrow(ax, STAGES[2] + WIDTH + 1.0, 13.0)
    x = STAGES[3]
    _badge(ax, x, 28.0, 4, "Measure the move", TEAL)
    ax.add_patch(FancyBboxPatch((x, 3.2), WIDTH, 21.0, boxstyle="round,pad=0.2",
                                facecolor="white", edgecolor=RULE, linewidth=1.0))
    ax.text(x + WIDTH / 2, 22.6, "model's map of cell states", ha="center",
            fontsize=7.5, color=INK2)
    # SCLC sits high enough to clear the caption underneath it.
    ax.add_patch(Ellipse((x + 6.4, 10.8), 10.6, 7.6, facecolor="#f0edf9",
                         edgecolor=HEMA, linewidth=1.2, linestyle="--"))
    ax.text(x + 6.4, 10.8, "SCLC", ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=HEMA)
    ax.add_patch(Ellipse((x + 13.6, 17.4), 10.6, 8.0, facecolor="#e9f7f1",
                         edgecolor=TEAL, linewidth=1.2, linestyle="--"))
    ax.text(x + 13.6, 17.4, "Normal", ha="center", va="center", fontsize=8.5,
            fontweight="bold", color=TEAL)
    ax.add_patch(Circle((x + 7.8, 13.0), 0.75, facecolor=HEMA, edgecolor="none",
                        zorder=6))
    ax.add_patch(FancyArrowPatch((x + 8.4, 13.4), (x + 12.2, 15.6),
                                 arrowstyle="-|>", mutation_scale=11,
                                 linewidth=1.8, color=EOSIN, zorder=6))
    ax.text(x + WIDTH / 2, 4.6, "shift toward Normal", ha="center",
            fontsize=7.5, style="italic", color=EOSIN)

    fig.subplots_adjust(left=0.004, right=0.996, top=0.99, bottom=0.01)
    out = OUT / "insilico_perturbation_explainer.png"
    fig.savefig(out, dpi=260)
    plt.close(fig)
    return out


if __name__ == "__main__":
    print(f"wrote {main()}")
