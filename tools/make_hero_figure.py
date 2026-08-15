#!/usr/bin/env python3
"""Build the poster's full-width key-message figure.

One banner spanning all three poster columns, in three acts, answering the
question a passer-by actually asks: what did you do, and why should a
computational screen be believed?

    A  The Perturb-seq analogy - what a pooled CRISPR screen does in cells,
       and what this does to the model's input instead.
    B  Where the effect lands - each deletion's effect on cell state, painted
       onto the known interaction map for the checkpoint genes.
    C  What survives - the funnel from 5.9M cell-level perturbations down to
       the four edits that replicate on every criterion.

Every number is read from a committed table; none is typed here.

    python tools/make_hero_figure.py

HONESTY CONSTRAINT, do not remove. Panel B paints *per-gene effects on cell
state* onto *prior-knowledge STRING edges*. This screen contains no gene-to-gene
measurements: it never estimates whether deleting gene X changes gene Y. Drawing
effects on a network invites exactly that misreading, which is why the edges are
drawn pale and unweighted, the panel is titled for the overlay rather than for
inference, and the caveat is printed inside the axes. A future edit that styles
the edges as though they were inferred, or that relabels this as a regulatory
network learned from the data, would make the figure claim something the
analysis did not measure.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
from matplotlib.colors import TwoSlopeNorm
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

REPO = Path(__file__).resolve().parent.parent
PRIM = REPO / "sclc_validation" / "primary_test_perturbation"
CKPT = REPO / "sclc_validation" / "checkpoint_cart_perturbation"
OUT = PRIM / "figures"
OUT.mkdir(parents=True, exist_ok=True)

# Spans the full A0 body width (817 mm inner). Everything is scaled up ~1.34x
# on the sheet, so a 9 pt label here prints near 12 pt. Height is the knob of
# last resort for the whole sheet: the band is full width, so every millimetre
# off this height is a millimetre off the poster, and --fit is already near
# its readability floor.
FIGSIZE = (24, 3.48)

INK = "#16202c"
INK2 = "#4a5768"
NAVY = "#10243c"
TEAL = "#0b6f6a"
EOSIN = "#b8465e"
HEMA = "#3d2f6b"
RULE = "#d5dee7"


# --------------------------------------------------------------------------
# numbers
# --------------------------------------------------------------------------
def _numbers() -> dict:
    arm = pd.read_csv(PRIM / "tables" / "primary_arm_summary.csv")
    hits = pd.read_csv(PRIM / "tables" / "allgene_concordant_hits_with_donor_robustness.csv")
    immune = pd.read_csv(PRIM / "tables" / "immune_cancer_candidates_with_donor_robustness.csv")
    cart = pd.read_csv(CKPT / "tables" / "cart_engineering_perturbation.csv")

    sclc_to_normal = cart[cart.comparison == "sclc_to_normal"]
    replicated = sclc_to_normal[
        sclc_to_normal.concordant
        & sclc_to_normal.tier.eq("all donors agree")
        & ~sclc_to_normal.low_detection_lt100
    ]
    return {
        # Each row of the arm summary is one arm x comparison; genes_tested is
        # the gene axis, so the sum is every gene-level test performed.
        "tests": int(arm.genes_tested.sum()),
        "genes_min": int(arm.genes_tested.min()),
        "genes_max": int(arm.genes_tested.max()),
        "hits": len(hits),
        "immune_genes": int(immune.Gene_name.nunique()),
        "programs": int(immune.class_label.nunique()),
        "replicated": len(replicated),
        "replicated_names": list(replicated.sort_values("delete_shift", ascending=False).Gene_name),
        "source_cells": int(cart.source_cells.iloc[0]),
    }


# --------------------------------------------------------------------------
# panel A - the analogy
# --------------------------------------------------------------------------
def _chip(ax, x, y, w, h, text, *, face, edge, fontsize=8.5, weight="normal",
          color=INK):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                facecolor=face, edgecolor=edge, linewidth=1.0))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
            fontsize=fontsize, fontweight=weight, color=color, linespacing=1.35)


def _arrow(ax, x0, x1, y, color="#93a3b5"):
    ax.add_patch(FancyArrowPatch((x0, y), (x1, y), arrowstyle="-|>",
                                 mutation_scale=13, linewidth=1.4, color=color))


def _panel_analogy(ax, num: dict) -> None:
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    ax.text(0, 0.965, "A", fontsize=14, fontweight="bold", color=INK, va="top")
    ax.text(0.038, 0.965, "Same question, no cells", fontsize=11.5,
            fontweight="bold", color=INK, va="top")

    steps_wet = ["guide RNA\nlibrary", "knock out\none gene", "scRNA-seq\nread-out",
                 "expression\nchange"]
    steps_sil = ["rank genes\nby expression", "drop one\ngene token", "Geneformer\nre-reads cell",
                 "cell-state\nshift"]

    width, gap = 0.178, 0.080
    for row, (label, steps, face, edge, tint) in enumerate([
        ("Perturb-seq  ·  in cells", steps_wet, "#f2f5f8", "#c3cedb", INK2),
        ("This work  ·  in silico", steps_sil, "#e9f2fa", "#8fbcd9", NAVY),
    ]):
        y = 0.58 - row * 0.34
        ax.text(0, y + 0.235, label, fontsize=9.5, fontweight="bold", color=tint)
        for index, step in enumerate(steps):
            x = index * (width + gap)
            _chip(ax, x, y, width, 0.185, step, face=face, edge=edge,
                  fontsize=8.5, color=tint)
            if index < len(steps) - 1:
                _arrow(ax, x + width + 0.012, x + width + gap - 0.012, y + 0.0925)

    # The contrast that makes the point: scale, not novelty of the question.
    ax.text(0.5, 0.02,
            f"A pooled CRISPR screen perturbs one gene per cell. Deleting a token costs "
            f"nothing, so every gene is deleted in every cell:\n"
            f"{num['genes_min']:,}–{num['genes_max']:,} genes × {num['source_cells']:,} "
            f"held-out cells × 6 disease transitions × 2 directions.",
            ha="center", va="bottom", fontsize=8.6, color=INK2, linespacing=1.5)


# --------------------------------------------------------------------------
# panel B - effect painted onto the known interaction map
# --------------------------------------------------------------------------
def _panel_network(ax) -> None:
    nodes = pd.read_csv(CKPT / "tables" / "network_node_perturbation.csv")
    edges = pd.read_csv(CKPT / "tables" / "string_network_edges.csv")

    graph = nx.Graph()
    graph.add_nodes_from(nodes.gene)
    for edge in edges.itertuples():
        graph.add_edge(edge.gene_a, edge.gene_b)

    degree = dict(zip(nodes.gene, nodes.string_degree))
    candidate = dict(zip(nodes.gene, nodes.is_candidate))
    connected = sorted([g for g in nodes.gene if degree[g] > 0],
                       key=lambda g: (not candidate[g], -degree[g], g))
    isolated = [g for g in nodes.gene if degree[g] == 0]

    positions = {}
    for index, gene in enumerate(connected):
        angle = np.pi / 2 - 2 * np.pi * index / len(connected)
        positions[gene] = (np.cos(angle), 0.84 * np.sin(angle))
    for index, gene in enumerate(isolated):
        positions[gene] = ((index - (len(isolated) - 1) / 2) * 0.6, -1.38)

    shifts = nodes.delete_sclc_to_normal.to_numpy(dtype=float)
    bound = float(np.abs(shifts).max())
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
    cmap = plt.get_cmap("RdBu_r")
    detection = nodes.detection_n_cells.to_numpy(dtype=float)
    sizes = 90 + 900 * (detection / detection.max())

    # Edges deliberately faint and unweighted: they are prior knowledge, not a
    # measured quantity, and must not read as an inferred regulatory graph.
    nx.draw_networkx_edges(graph, positions, ax=ax, edge_color="#dbe2ea",
                           width=0.7, alpha=0.9)
    nx.draw_networkx_nodes(
        graph, positions, ax=ax, nodelist=list(nodes.gene), node_size=sizes,
        node_color=[cmap(norm(v)) for v in shifts],
        edgecolors=["#1f2937" if c else "#9ca3af" for c in nodes.is_candidate],
        linewidths=[1.1 if c else 0.6 for c in nodes.is_candidate],
    )
    for gene, size, cand, alone in zip(nodes.gene, sizes, nodes.is_candidate,
                                       [g in isolated for g in nodes.gene]):
        x, y = positions[gene]
        if alone:
            ux, uy = 0.0, -1.0
        else:
            radius = float(np.hypot(x, y)) or 1.0
            ux, uy = x / radius, y / radius
        offset = np.sqrt(size / np.pi) + 5
        ax.annotate(gene, (x, y), xytext=(ux * offset, uy * offset),
                    textcoords="offset points",
                    ha="left" if ux > 0.35 else ("right" if ux < -0.35 else "center"),
                    va="bottom" if uy > 0.35 else ("top" if uy < -0.35 else "center"),
                    fontsize=7.6 if cand else 6.8,
                    fontweight="bold" if cand else "normal",
                    color=INK if cand else "#8b96a5")

    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-2.05, 1.5)
    ax.set_axis_off()
    ax.text(-1.75, 1.45, "B", fontsize=14, fontweight="bold", color=INK, va="top")
    ax.text(-1.60, 1.45, "Where the effect lands", fontsize=11.5,
            fontweight="bold", color=INK, va="top")
    ax.text(0, -1.78,
            "Colour = measured shift toward Normal when that gene is deleted.\n"
            "Edges are published STRING interactions shown for context — this\n"
            "screen measures each gene against cell state, not gene against gene.",
            ha="center", va="top", fontsize=7.6, color=INK2, linespacing=1.5)

    bar = ax.inset_axes([0.60, 0.90, 0.37, 0.030])
    colorbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=bar,
                            orientation="horizontal")
    colorbar.set_ticks([-bound, 0, bound])
    colorbar.set_ticklabels(["away", "0", "toward Normal"])
    colorbar.ax.tick_params(labelsize=6.8, pad=1.5)
    colorbar.outline.set_linewidth(0.4)


# --------------------------------------------------------------------------
# panel C - the funnel
# --------------------------------------------------------------------------
def _panel_funnel(ax, num: dict) -> None:
    stages = [
        (num["tests"], "gene × transition × direction tests",
         "every expressed gene, both arms", "#cddcea"),
        (num["hits"], "concordant hits",
         "deletion and overexpression disagree in sign, FDR < 0.05 in both", "#9dbfda"),
        (num["immune_genes"], f"genes across {num['programs']} immune / cancer programs",
         "after removing ribosomal, mitochondrial and ambient classes", "#5f93bd"),
        (num["replicated"], "replicate on every criterion",
         " · ".join(num["replicated_names"]), TEAL),
    ]

    ax.set_axis_off()
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.text(0, 0.965, "C", fontsize=14, fontweight="bold", color=INK, va="top")
    ax.text(0.033, 0.965, "What survives", fontsize=11.5, fontweight="bold",
            color=INK, va="top")

    # Label above, bar below. Putting the count and headline inside the bar
    # only works while the bar is wide, and this funnel drops five orders of
    # magnitude - the last bar cannot hold its own label. Text above, bar as a
    # pure magnitude strip, keeps every stage legible at the same size.
    #
    # Widths are log-scaled for the same reason: linearly, 4 against 165,438
    # would be a bar less than a pixel wide.
    top = np.log10(stages[0][0])
    for index, (value, headline, detail, colour) in enumerate(stages):
        base = 0.78 - index * 0.215
        width = 0.16 + 0.82 * (np.log10(value) / top)

        ax.text(0.052, base, f"{value:,}", va="baseline", ha="left",
                fontsize=12.5, fontweight="bold", color=NAVY)
        ax.text(0.052 + 0.008 + 0.115 * (len(f"{value:,}") / 7 + 0.35), base,
                headline, va="baseline", ha="left", fontsize=8.8,
                fontweight="bold", color=NAVY)
        ax.add_patch(FancyBboxPatch((0.052, base - 0.058), width, 0.032,
                                    boxstyle="round,pad=0.003",
                                    facecolor=colour, edgecolor="none"))
        ax.text(0.052, base - 0.098, detail, va="baseline", ha="left",
                fontsize=7.6, color=INK2)

    ax.annotate("", xy=(0.022, 0.60), xytext=(0.022, 0.80),
                arrowprops=dict(arrowstyle="-|>", color="#b6c2d0", linewidth=1.3))
    ax.annotate("", xy=(0.022, 0.17), xytext=(0.022, 0.38),
                arrowprops=dict(arrowstyle="-|>", color="#b6c2d0", linewidth=1.3))


def main() -> Path:
    num = _numbers()
    fig = plt.figure(figsize=FIGSIZE)
    grid = fig.add_gridspec(1, 3, width_ratios=[1.32, 0.86, 1.06],
                            left=0.012, right=0.988, top=0.94, bottom=0.05,
                            wspace=0.10)
    _panel_analogy(fig.add_subplot(grid[0, 0]), num)
    _panel_network(fig.add_subplot(grid[0, 1]))
    _panel_funnel(fig.add_subplot(grid[0, 2]), num)

    out = OUT / "insilico_perturbseq_hero.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    path = main()
    print(f"wrote {path}")
