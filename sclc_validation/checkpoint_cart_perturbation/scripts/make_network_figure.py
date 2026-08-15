#!/usr/bin/env python3
"""Rebuild the STRING interaction-map figure (poster Fig. 5).

Renders `figures/perturbation_networks.png`: one STRING map drawn three times
under three perturbation read-outs, on a shared node layout, a shared colour
scale and a shared size scale, so that only colour differs panel to panel.
That is the whole point of the figure - it is what shows the SCLC->Normal
versus SCLC->LUAD sign flip directly.

Reads only committed tables, so it reproduces without the remote compute
artifacts and without re-querying STRING:

    tables/network_node_perturbation.csv   16 genes, all four read-outs
    tables/string_network_edges.csv        54 non-text-mining edges

Run from anywhere:

    python sclc_validation/checkpoint_cart_perturbation/scripts/make_network_figure.py

PROVENANCE. This script was written after the fact, to reproduce a figure that
was originally produced by an uncommitted notebook. It regenerates every
encoding from the tables - colour, node area, rings, labels, the edge filter
and the colour-scale bounds are all derived, not transcribed. What it cannot
recover is the original node layout: those coordinates were never saved. The
ring layout here is deterministic (see `_layout`), so the figure is stable run
to run, but it is NOT pixel-identical to the previously committed PNG. Every
quantity the figure encodes is identical; only where a node sits differs, and
node positions carry no meaning in this figure - the same layout is reused
across all three panels precisely so that position is held constant.
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
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
FIGURES.mkdir(parents=True, exist_ok=True)

# The three read-outs drawn, in panel order. delete_sclc_to_luad is present in
# the table but deliberately not shown: three panels already carry the sign-flip
# comparison, and a fourth makes every node too small to read at poster size.
READOUTS = [
    ("delete_sclc_to_normal", "Deletion → Normal"),
    ("overexpress_sclc_to_normal", "Overexpression → Normal"),
    ("overexpress_sclc_to_luad", "Overexpression → LUAD"),
]

LOW_DETECTION = 100      # cells; matches low_detection_lt100 in the tables

# Sized for the poster's middle column (299 mm inner width in Draft 12).
# Node areas and label sizes are in points, so they do NOT scale with the
# figure: shrinking the canvas is what makes them render larger once the PNG
# is scaled to the column width. 13 in put the gene labels near 10 pt on the
# A0 sheet; 10.5 in puts them near 12.5 pt, which is what the middle column's
# spare height was spent on.
FIGSIZE = (10.5, 6.4)


def _load() -> tuple[pd.DataFrame, pd.DataFrame]:
    nodes = pd.read_csv(TABLES / "network_node_perturbation.csv")
    edges = pd.read_csv(TABLES / "string_network_edges.csv")
    return nodes, edges


def _build_graph(nodes: pd.DataFrame, edges: pd.DataFrame) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from(nodes.gene)
    for edge in edges.itertuples():
        # Edges are already filtered upstream to drop text-mining-only support;
        # combined_score is kept only to weight the drawn line.
        graph.add_edge(edge.gene_a, edge.gene_b, weight=float(edge.combined_score))
    return graph


def _layout(nodes: pd.DataFrame) -> dict[str, tuple[float, float]]:
    """Ring layout: connected genes on an ellipse, isolated genes beneath it.

    A force-directed layout is the obvious choice and the wrong one here. With
    16 nodes and 54 dense edges it converges to a single tight clump in which
    the large detection-scaled nodes overlap and their labels are unreadable -
    the whole figure depends on being able to read node colour per gene. A ring
    is deterministic without a seed, guarantees no node overlaps another, and
    keeps every label in clear space. Candidates are placed first so they
    occupy a contiguous arc rather than being interleaved with context genes.
    """
    degree = dict(zip(nodes.gene, nodes.string_degree))
    connected = [g for g in nodes.gene if degree[g] > 0]
    isolated = [g for g in nodes.gene if degree[g] == 0]

    candidate = dict(zip(nodes.gene, nodes.is_candidate))
    connected.sort(key=lambda g: (not candidate[g], -degree[g], g))

    positions: dict[str, tuple[float, float]] = {}
    count = len(connected)
    for index, gene in enumerate(connected):
        # Start at the top and run clockwise so the highest-degree candidate
        # sits top-left, as in the original figure.
        angle = np.pi / 2 - 2 * np.pi * index / count
        positions[gene] = (1.00 * np.cos(angle), 0.86 * np.sin(angle))

    # Isolated genes get their own row below the ring; they are part of the
    # result (TOX and LAYN have no qualifying STRING edge) and are labelled as
    # such rather than dropped.
    for index, gene in enumerate(isolated):
        offset = (index - (len(isolated) - 1) / 2) * 0.62
        positions[gene] = (offset, -1.34)
    return positions


def _node_sizes(nodes: pd.DataFrame) -> np.ndarray:
    """Node area proportional to detection, so area reads as cell count."""
    detection = nodes.detection_n_cells.to_numpy(dtype=float)
    return 120 + 2600 * (detection / detection.max())


def main() -> Path:
    nodes, edges = _load()
    graph = _build_graph(nodes, edges)
    positions = _layout(nodes)

    # One symmetric colour scale across all three panels - panels are only
    # comparable if they share it. Bound is the largest absolute shift shown.
    values = nodes[[column for column, _ in READOUTS]].to_numpy(dtype=float)
    bound = float(np.abs(values).max())
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
    cmap = plt.get_cmap("RdBu_r")

    sizes = _node_sizes(nodes)
    order = list(nodes.gene)
    low_detection = nodes.low_detection_lt100.to_numpy(dtype=bool)
    is_candidate = nodes.is_candidate.to_numpy(dtype=bool)
    isolated = nodes.string_degree.to_numpy(dtype=int) == 0

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE)

    for ax, (column, title), letter in zip(axes, READOUTS, "abc"):
        shifts = nodes[column].to_numpy(dtype=float)

        nx.draw_networkx_edges(
            graph, positions, ax=ax, edge_color="#cbd5e1", alpha=0.75,
            width=[0.5 + 2.4 * graph[u][v]["weight"] ** 3 for u, v in graph.edges()],
        )
        # Ring encodes two things that are otherwise invisible: a candidate
        # gene versus network context, and detection too low to trust.
        edge_colors = ["#c2691a" if low else ("#1f2937" if cand else "#9ca3af")
                       for low, cand in zip(low_detection, is_candidate)]
        edge_widths = [2.4 if low else (1.5 if cand else 0.8) for low, cand
                       in zip(low_detection, is_candidate)]
        nx.draw_networkx_nodes(
            graph, positions, ax=ax, nodelist=order,
            node_size=sizes, node_color=[cmap(norm(v)) for v in shifts],
            edgecolors=edge_colors, linewidths=edge_widths,
        )

        for gene, size, cand, count, alone in zip(order, sizes, is_candidate,
                                                  nodes.detection_n_cells, isolated):
            x, y = positions[gene]
            # Push each label radially outward from the ring centre rather than
            # straight down. Straight-down collided along the bottom of the
            # ring, where several nodes sit at similar y - "TCF7" landed on
            # "HAVCR2". Radially, every label leaves in its own direction.
            if alone:
                ux, uy = 0.0, -1.0          # isolated row: no meaningful radius
            else:
                radius = float(np.hypot(x, y)) or 1.0
                ux, uy = x / radius, y / radius
            offset = np.sqrt(size / np.pi) + 7
            ax.annotate(
                f"{gene}\n{int(count):,}", (x, y),
                xytext=(ux * offset, uy * offset), textcoords="offset points",
                ha="left" if ux > 0.35 else ("right" if ux < -0.35 else "center"),
                va="bottom" if uy > 0.35 else ("top" if uy < -0.35 else "center"),
                fontsize=11 if cand else 9.5,
                fontweight="bold" if cand else "normal",
                color="#111827" if cand else "#6b7280",
            )

        ax.set_title(title, fontsize=13, pad=14)
        ax.text(0.0, 1.06, letter, transform=ax.transAxes,
                fontsize=15, fontweight="bold", va="bottom")
        ax.set_axis_off()
        ax.margins(0.20)

    note = ", ".join(sorted(nodes.gene[isolated]))
    if note:
        fig.text(0.5, 0.165, f"{note}: no STRING edge above threshold",
                 ha="center", fontsize=10.5, style="italic", color="#6b7280")

    fig.suptitle(
        "Same STRING interaction map, three perturbation read-outs\n"
        "node area = detection (cells expressing the gene, of "
        f"{int(nodes.source_cells.iloc[0]):,})  ·  label shows the count",
        fontsize=11.5, y=0.995, linespacing=1.35,
    )

    # ---- shared colour bar -------------------------------------------------
    bar_ax = fig.add_axes([0.06, 0.075, 0.34, 0.022])
    colorbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=bar_ax, orientation="horizontal")
    colorbar.set_ticks([-bound, 0, bound])
    colorbar.set_ticklabels([f"{-bound:.3f}", "0", f"+{bound:.3f}"])
    colorbar.ax.tick_params(labelsize=10)
    colorbar.set_label("Goal-state shift  (blue = away · red = toward goal)",
                       fontsize=10.5)

    # ---- shared size key ---------------------------------------------------
    detection_max = float(nodes.detection_n_cells.max())
    handles, labels = [], []
    for count in (10, 100, 1000):
        handles.append(Line2D([], [], marker="o", linestyle="none",
                              markersize=np.sqrt(120 + 2600 * (count / detection_max)) * 0.62,
                              markerfacecolor="#dfe6ee", markeredgecolor="#9ca3af",
                              markeredgewidth=0.8))
        labels.append(f"{count:,}")
    handles.append(Line2D([], [], marker="o", linestyle="none", markersize=8,
                          markerfacecolor="#dfe6ee", markeredgecolor="#c2691a",
                          markeredgewidth=2.2))
    labels.append(f"<{LOW_DETECTION} cells")
    legend = fig.legend(handles, labels, loc="lower right",
                        bbox_to_anchor=(0.99, 0.025), ncol=4, frameon=False,
                        fontsize=10.5, handletextpad=0.5, columnspacing=1.4,
                        title="detection (cells expressing gene)")
    legend.get_title().set_fontsize(10.5)

    fig.subplots_adjust(left=0.045, right=0.955, top=0.86, bottom=0.22, wspace=0.34)
    out = FIGURES / "perturbation_networks.png"
    fig.savefig(out, dpi=220)
    plt.close(fig)
    return out


if __name__ == "__main__":
    path = main()
    print(f"wrote {path}")
