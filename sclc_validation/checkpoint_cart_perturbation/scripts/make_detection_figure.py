#!/usr/bin/env python3
"""Rebuild the detection-versus-effect figure (poster Fig. 6).

Two panels, both making the same point from different angles: in this screen the
largest nominal deletion effects come from the genes detected in the fewest
cells, so effect magnitude on its own is not a ranking criterion.

    a  detection against |deletion shift| for the whole 50-gene panel, log-log,
       with the 100-cell threshold drawn
    b  every immune / lineage gene ranked by deletion shift, each labelled with
       its detection count

Reads one committed table:

    tables/cart_overexpression_vs_deletion.csv   50 genes, sclc_to_normal

Run from anywhere:

    python sclc_validation/checkpoint_cart_perturbation/scripts/make_detection_figure.py

PROVENANCE. Written after the fact to reproduce a figure whose notebook was
never committed. Every quantity is recomputed from the table rather than
transcribed, and one headline number in the previously committed PNG does not
survive that. Its gene count is right - the table does hold 31 immune / lineage
genes - but it reported "TIM-3 5th, TIGIT 6th of 24". Applying the threshold the
caption states (detection >= 100 cells) selects 20 genes and ranks TIM-3 4th and
TIGIT 5th. The poster's ICI/CAR-T panel already quotes 4th and 5th of 20, so
this brings the figure into line with the sheet rather than the reverse. The
same stale rank appears in this package's README.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

LOW_DETECTION = 100          # cells; matches low_detection_lt100 in the tables

BLUE = "#17608f"             # replicated in all donors
BLUE_OPEN = "#2b7cb3"        # other CAR-T candidate
GREY = "#b7c2ce"             # other immune gene
ORANGE = "#c2691a"           # detected in too few cells
INK = "#16202c"
INK2 = "#4a5768"

# Sized for the poster's left column (242 mm inner). Text is in points and does
# not scale with the canvas, so canvas width sets the rendered text size: at
# 13.5 in the 31 gene labels in panel b came out near 7 pt on the sheet, under
# the poster's 9 pt readability floor. 9.5 in renders them ~1:1 instead. Panel b
# also needs real height - 31 labels at >=9 pt cannot fit in less than about
# 100 mm however wide the canvas is.
FIGSIZE = (9.5, 6.3)


def _load() -> tuple[pd.DataFrame, list[str]]:
    panel = pd.read_csv(TABLES / "cart_overexpression_vs_deletion.csv")
    cart = pd.read_csv(TABLES / "cart_engineering_perturbation.csv")
    sclc = cart[cart.comparison == "sclc_to_normal"]
    replicated = (sclc[sclc.concordant
                       & sclc.tier.eq("all donors agree")
                       & ~sclc.low_detection_lt100]
                  .sort_values("delete_shift", ascending=False))
    return panel, list(replicated.Gene_name)


def _panel_scatter(ax, panel: pd.DataFrame, replicated: list[str]) -> None:
    testable = panel[panel.deletion_testable]
    rho, pval = spearmanr(testable.detection_n_cells, testable.delete_shift.abs())

    for _, row in testable.iterrows():
        x, y = row.detection_n_cells, abs(row.delete_shift)
        if row.Gene_name in replicated:
            ax.scatter(x, y, s=190, c=BLUE, zorder=5, linewidths=0)
        elif row.low_detection_lt100:
            ax.scatter(x, y, s=130, facecolors="none", edgecolors=ORANGE,
                       linewidths=2.0, zorder=4)
        elif row.is_cart_candidate:
            ax.scatter(x, y, s=130, facecolors="none", edgecolors=BLUE_OPEN,
                       linewidths=2.0, zorder=4)
        else:
            ax.scatter(x, y, s=80, c=GREY, zorder=2, linewidths=0)

    ax.axvline(LOW_DETECTION, color=ORANGE, linestyle="--", linewidth=1.6, zorder=1)
    ax.text(LOW_DETECTION * 1.18, testable.delete_shift.abs().max() * 0.72,
            f"detection\nthreshold\n{LOW_DETECTION} cells",
            fontsize=11, color=ORANGE, va="top")

    # Label the four replicated hits and the sparse genes that dominate the top,
    # which together are the whole argument of the panel.
    sparse = (testable[testable.low_detection_lt100]
              .assign(_abs=lambda d: d.delete_shift.abs())
              .nlargest(5, "_abs"))
    # Stagger the offsets: the four hits sit at similar x and y and a single
    # fixed offset stacked their labels on one another.
    offsets = [(10, 9), (10, -16), (-12, 11), (-12, -18)]
    labelled = pd.concat([testable[testable.Gene_name.isin(replicated)], sparse])
    for index, (_, row) in enumerate(labelled.iterrows()):
        hit = row.Gene_name in replicated
        dx, dy = offsets[index % len(offsets)] if hit else (9, 7)
        ax.annotate(row.Gene_name, (row.detection_n_cells, abs(row.delete_shift)),
                    xytext=(dx, dy), textcoords="offset points",
                    ha="left" if dx > 0 else "right",
                    fontsize=12 if hit else 11,
                    fontweight="bold" if hit else "normal",
                    color=INK if hit else ORANGE)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(f"Detection: cells expressing the gene "
                  f"(of {int(panel.source_cells.iloc[0]):,} SCLC source cells)", fontsize=13)
    ax.set_ylabel("|Deletion shift|", fontsize=13)
    ax.set_title("Largest effects come from\nthe least-detected genes\n"
                 f"(Spearman ρ = {rho:.2f}, p = {pval:.0e}, n = {len(testable)})",
                 fontsize=10.5, loc="left")
    ax.tick_params(labelsize=11)
    ax.grid(alpha=0.15)
    ax.text(-0.085, 1.19, "a", transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top")


def _panel_rank(ax, panel: pd.DataFrame, replicated: list[str]) -> None:
    immune = panel[~panel.technical_or_ambient].copy()
    immune = immune.sort_values("delete_shift")
    detected = immune[immune.detection_n_cells >= LOW_DETECTION]
    ranked = detected.sort_values("delete_shift", ascending=False).reset_index(drop=True)
    rank = {g: i + 1 for i, g in enumerate(ranked.Gene_name)}

    y = np.arange(len(immune))
    for pos, (_, row) in zip(y, immune.iterrows()):
        testable = bool(row.deletion_testable)
        if not testable:
            ax.text(0.0, pos, "  n.d.", va="center", fontsize=10,
                    style="italic", color="#9aa8b8")
            continue
        ax.plot([0, row.delete_shift], [pos, pos], color="#dbe2ea",
                linewidth=1.4, zorder=1)
        detection = 0.0 if pd.isna(row.detection_n_cells) else row.detection_n_cells
        size = 60 + 260 * (detection / immune.detection_n_cells.max()) ** 0.6
        if row.Gene_name in replicated:
            ax.scatter(row.delete_shift, pos, s=size, c=BLUE, zorder=3, linewidths=0)
        elif row.low_detection_lt100:
            ax.scatter(row.delete_shift, pos, s=size, facecolors="none",
                       edgecolors=ORANGE, linewidths=1.8, zorder=3)
        elif row.is_cart_candidate:
            ax.scatter(row.delete_shift, pos, s=size, facecolors="none",
                       edgecolors=BLUE_OPEN, linewidths=1.8, zorder=3)
        else:
            ax.scatter(row.delete_shift, pos, s=size, c=GREY, zorder=2, linewidths=0)

    labels, weights, colours = [], [], []
    for _, row in immune.iterrows():
        count = row.detection_n_cells
        shown = "0" if pd.isna(count) else f"{int(count):,}"
        labels.append(f"{row.Gene_name} ({shown})")
        hit = row.Gene_name in replicated
        weights.append("bold" if hit else "normal")
        colours.append(INK if hit else INK2)
    ax.set_yticks(y, labels, fontsize=9)
    for tick, weight, colour in zip(ax.get_yticklabels(), weights, colours):
        tick.set_fontweight(weight)
        tick.set_color(colour)

    ax.axvline(0, color="#8b96a5", linewidth=1.0)
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_xlabel("Deletion shift, SCLC → Normal (symlog)\n(right = toward normal state)",
                  fontsize=13)
    top = ", ".join(f"{g} {rank[g]}{'st' if rank[g]==1 else 'nd' if rank[g]==2 else 'rd' if rank[g]==3 else 'th'}"
                    for g in replicated[:2])
    # Keep each line under ~48 characters: the axes is about 4.3 in wide here
    # and matplotlib does not wrap or clip a title, it just runs off the canvas.
    ax.set_title(f"All {len(immune)} immune / lineage genes\n"
                 f"(count in brackets)\n"
                 f"{len(detected)} detected ≥{LOW_DETECTION} cells — {top}",
                 fontsize=10.5, loc="left")
    ax.tick_params(axis="x", labelsize=11)
    ax.grid(axis="x", alpha=0.15)
    ax.text(-0.34, 1.19, "b", transform=ax.transAxes, fontsize=15,
            fontweight="bold", va="top")


def main() -> Path:
    panel, replicated = _load()
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE,
                             gridspec_kw={"width_ratios": [1.0, 1.0], "wspace": 0.34})
    _panel_scatter(axes[0], panel, replicated)
    _panel_rank(axes[1], panel, replicated)

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=11,
               markerfacecolor=BLUE, markeredgecolor="none",
               label=f"replicated in all 3 donors"),
        Line2D([], [], marker="o", linestyle="none", markersize=11,
               markerfacecolor="none", markeredgecolor=BLUE_OPEN, markeredgewidth=2,
               label="other CAR-T candidate"),
        Line2D([], [], marker="o", linestyle="none", markersize=10,
               markerfacecolor=GREY, markeredgecolor="none", label="other immune gene"),
        Line2D([], [], marker="o", linestyle="none", markersize=11,
               markerfacecolor="none", markeredgecolor=ORANGE, markeredgewidth=2,
               label=f"detected in <{LOW_DETECTION} cells — effect unreliable"),
        Line2D([], [], marker="|", linestyle="none", markersize=11, color="#9aa8b8",
               label="n.d. — undetected, deletion undefined"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=2, frameon=False,
               fontsize=9.5, bbox_to_anchor=(0.5, -0.005))

    fig.subplots_adjust(left=0.105, right=0.985, top=0.80, bottom=0.22)
    out = FIGURES / "cart_overexpression.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    print(f"wrote {main()}")
