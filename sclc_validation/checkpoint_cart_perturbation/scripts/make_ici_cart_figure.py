#!/usr/bin/env python3
"""Rebuild the ICI / CAR-T candidate figure (poster Fig. 8).

Three views of the same 11-gene candidate set:

    a  the 5 ICI drug targets across all 6 disease transitions, marked by how
       far their donor support goes - which is where "8 of 17 concordant calls
       rest on a single donor" comes from
    b  restricted to SCLC -> Normal, every CAR-T engineering gene ranked by
       deletion shift, with the four that replicate across all 3 donors marked
    c  the same 11 candidates placed on the STRING interaction map, coloured by
       their SCLC -> Normal deletion shift

Reads three committed tables:

    tables/ici_target_perturbation.csv          5 ICI genes x 6 comparisons
    tables/cart_engineering_perturbation.csv    11 CAR-T genes x 6 comparisons
    tables/string_network_edges.csv             54 non-text-mining edges

Run from anywhere:

    python sclc_validation/checkpoint_cart_perturbation/scripts/make_ici_cart_figure.py

PROVENANCE. Written after the fact; the notebook that first produced this figure
was never committed. Every value and count is recomputed from the tables. As in
make_network_figure.py, the original node coordinates for panel c were never
saved; it uses the same deterministic ring layout, so positions differ from the
previously committed PNG while every encoded quantity is identical.

PD-L1 (CD274) is drawn as an explicit greyed marker rather than omitted. It is a
requested ICI target that cannot appear in this screen at all: the atlas is
T-cell only and PD-L1 is expressed on tumour and antigen-presenting cells, so it
was never tokenized as a source gene. Silently dropping it would read as "not a
hit" instead of "not measurable here".
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

GOAL = "sclc_to_normal"

BLUE = "#17608f"
BLUE_LT = "#8fbcd9"
ORANGE = "#c2691a"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

# Display names for the drug targets; the tables carry HGNC symbols.
DRUG_NAME = {"CTLA4": "CTLA-4", "PDCD1": "PD-1", "LAG3": "LAG-3",
             "HAVCR2": "TIM-3", "TIGIT": "TIGIT"}

# Sized for the poster's middle column (299 mm inner width).
FIGSIZE = (16, 5.6)


def _replicated(cart: pd.DataFrame) -> list[str]:
    goal = cart[cart.comparison == GOAL]
    hits = goal[goal.concordant & goal.tier.eq("all donors agree")
                & ~goal.low_detection_lt100]
    return list(hits.sort_values("delete_shift", ascending=False).Gene_name)


def _panel_ici(ax, ici: pd.DataFrame) -> pd.DataFrame:
    """Every ICI target across every transition, marked by donor support."""
    genes = sorted(ici.Gene_name.unique(),
                   key=lambda g: -ici[ici.Gene_name == g].delete_shift.abs().max())
    style = {
        "all donors agree": dict(s=170, c=BLUE, linewidths=0, zorder=5),
        "majority agree":   dict(s=150, c=BLUE_LT, linewidths=0, zorder=4),
        "single donor":     dict(s=150, facecolors="none", edgecolors=ORANGE,
                                 linewidths=2.2, zorder=4),
        "donors disagree":  dict(s=110, c="#9aa8b8", marker="X", linewidths=0, zorder=3),
        "not concordant":   dict(s=26, c=GREY, linewidths=0, zorder=2),
    }
    for position, gene in enumerate(genes):
        for _, row in ici[ici.Gene_name == gene].iterrows():
            ax.scatter(row.delete_shift, position,
                       **style.get(row.tier, style["not concordant"]))

    ax.axvline(0, color="#8b96a5", linewidth=1.0)
    ax.set_yticks(range(len(genes)),
                  [f"{DRUG_NAME.get(g, g)}\n({g})" for g in genes], fontsize=11)
    ax.set_xlabel("Change in similarity to the goal state after deleting the gene\n"
                  "(6 source→goal transitions per gene;  right = toward goal)",
                  fontsize=11.5)
    single = int((ici.tier == "single donor").sum())
    concordant = int(ici.concordant.sum())
    ax.set_title("Checkpoint deletion moves T cells\nbetween disease states; "
                 f"{single} of {concordant}\nconcordant calls rest on one donor",
                 fontsize=11, loc="left")
    ax.tick_params(axis="x", labelsize=10)
    ax.xaxis.set_major_locator(plt.MaxNLocator(5))
    ax.grid(axis="x", alpha=0.15)
    ax.text(-0.28, 1.14, "a", transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="top")
    return ici


def _panel_rank(ax, cart: pd.DataFrame, replicated: list[str]) -> None:
    goal = cart[cart.comparison == GOAL].sort_values("delete_shift")
    positions = np.arange(len(goal))
    for position, (_, row) in zip(positions, goal.iterrows()):
        ax.plot([0, row.delete_shift], [position, position],
                color="#dbe2ea", linewidth=1.6, zorder=1)
        hit = row.Gene_name in replicated
        ax.scatter(row.delete_shift, position, s=190 if hit else 90,
                   c=BLUE if hit else GREY, linewidths=0, zorder=3)

    best = goal.iloc[-1]
    ax.annotate(f"{best.delete_shift:+.4f}",
                (best.delete_shift, len(goal) - 1), xytext=(10, 0),
                textcoords="offset points", va="center",
                fontsize=11.5, fontweight="bold", color=BLUE)

    ax.axvline(0, color="#8b96a5", linewidth=1.0)
    ax.set_yticks(positions, goal.Gene_name, fontsize=11)
    for tick, gene in zip(ax.get_yticklabels(), goal.Gene_name):
        if gene in replicated:
            tick.set_fontweight("bold")
            tick.set_color(INK)
    ax.set_xlabel("Change in similarity to the normal-T-cell state\n"
                  "after deleting the gene in SCLC T cells  (right = toward normal)",
                  fontsize=11.5)
    ax.set_title("SCLC T cells → normal T-cell state:\n"
                 f"{len(replicated)} knockouts replicate across\nall 3 SCLC donors",
                 fontsize=11, loc="left")
    ax.tick_params(axis="x", labelsize=10)
    ax.grid(axis="x", alpha=0.15)
    ax.text(-0.24, 1.14, "b", transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="top")


def _panel_network(ax, nodes: pd.DataFrame, edges: pd.DataFrame,
                   candidates: set[str]) -> None:
    graph = nx.Graph()
    graph.add_nodes_from(nodes.gene)
    for edge in edges.itertuples():
        graph.add_edge(edge.gene_a, edge.gene_b)

    degree = dict(zip(nodes.gene, nodes.string_degree))
    connected = sorted([g for g in nodes.gene if degree[g] > 0],
                       key=lambda g: (g not in candidates, -degree[g], g))
    isolated = [g for g in nodes.gene if degree[g] == 0]

    positions = {}
    for index, gene in enumerate(connected):
        angle = np.pi / 2 - 2 * np.pi * index / len(connected)
        positions[gene] = (np.cos(angle), 0.86 * np.sin(angle))
    for index, gene in enumerate(isolated):
        positions[gene] = (-2.05, 0.45 - index * 0.55)

    shifts = nodes.delete_sclc_to_normal.to_numpy(dtype=float)
    bound = float(np.abs(shifts).max())
    norm = TwoSlopeNorm(vmin=-bound, vcenter=0.0, vmax=bound)
    cmap = plt.get_cmap("RdBu_r")

    nx.draw_networkx_edges(graph, positions, ax=ax, edge_color="#dbe2ea",
                           width=0.8, alpha=0.9)
    nx.draw_networkx_nodes(
        graph, positions, ax=ax, nodelist=list(nodes.gene), node_size=760,
        node_color=[cmap(norm(v)) for v in shifts],
        edgecolors=["#1f2937" if g in candidates else "#9ca3af" for g in nodes.gene],
        linewidths=[1.6 if g in candidates else 0.7 for g in nodes.gene],
    )
    for gene in nodes.gene:
        x, y = positions[gene]
        candidate = gene in candidates
        ax.annotate(DRUG_NAME.get(gene, gene), (x, y), ha="center", va="center",
                    fontsize=8.4 if candidate else 7.4,
                    fontweight="bold" if candidate else "normal",
                    color=INK if candidate else INK2)

    if isolated:
        ax.text(-2.05, 0.45 - len(isolated) * 0.55, "no STRING edge\nabove threshold",
                ha="center", va="top", fontsize=8.6, style="italic", color=INK2)
    # See the module docstring: absence here means unmeasurable, not null.
    ax.scatter([-2.05], [-1.30], s=620, facecolors="#eef1f5", edgecolors="#9ca3af",
               linewidths=1.0, linestyle="--")
    ax.annotate("PD-L1", (-2.05, -1.30), ha="center", va="center",
                fontsize=8.0, color="#9aa8b8")
    ax.text(-2.05, -1.62, "not in T-cell\natlas", ha="center", va="top",
            fontsize=8.6, style="italic", color=INK2)

    ax.set_xlim(-2.55, 1.65)
    ax.set_ylim(-2.05, 1.45)
    ax.set_axis_off()
    ax.set_title(f"All {len(candidates)} candidates on the STRING map;\n"
                 "bold = panel b candidate,\ncolour = its SCLC→Normal shift",
                 fontsize=11, loc="left")
    ax.text(0.0, 1.14, "c", transform=ax.transAxes, fontsize=16,
            fontweight="bold", va="top")

    bar = ax.inset_axes([0.30, 0.02, 0.46, 0.035])
    colorbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap),
                            cax=bar, orientation="horizontal")
    colorbar.set_ticks([-bound, 0, bound])
    colorbar.set_ticklabels([f"{-bound:.3f}", "0", f"+{bound:.3f}"])
    colorbar.ax.tick_params(labelsize=8.5, pad=1.5)
    colorbar.set_label("Deletion shift, SCLC→Normal", fontsize=9)


def main() -> Path:
    ici = pd.read_csv(TABLES / "ici_target_perturbation.csv")
    cart = pd.read_csv(TABLES / "cart_engineering_perturbation.csv")
    nodes = pd.read_csv(TABLES / "network_node_perturbation.csv")
    edges = pd.read_csv(TABLES / "string_network_edges.csv")
    replicated = _replicated(cart)
    candidates = set(nodes[nodes.is_candidate].gene)

    fig, axes = plt.subplots(1, 3, figsize=FIGSIZE,
                             gridspec_kw={"width_ratios": [1.05, 1.0, 1.15],
                                          "wspace": 0.52})
    _panel_ici(axes[0], ici)
    _panel_rank(axes[1], cart, replicated)
    _panel_network(axes[2], nodes, edges, candidates)

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=11,
               markerfacecolor=BLUE, markeredgecolor="none", label="all donors agree"),
        Line2D([], [], marker="o", linestyle="none", markersize=10,
               markerfacecolor=BLUE_LT, markeredgecolor="none", label="majority agree"),
        Line2D([], [], marker="o", linestyle="none", markersize=11,
               markerfacecolor="none", markeredgecolor=ORANGE, markeredgewidth=2,
               label="single donor"),
        Line2D([], [], marker="X", linestyle="none", markersize=10,
               markerfacecolor="#9aa8b8", markeredgecolor="none", label="donors disagree"),
        Line2D([], [], marker="o", linestyle="none", markersize=6,
               markerfacecolor=GREY, markeredgecolor="none", label="not concordant"),
    ]
    fig.legend(handles=handles, loc="lower left", ncol=3, frameon=False,
               fontsize=10.5, bbox_to_anchor=(0.035, -0.01),
               title="Donor replication", title_fontsize=10.5)

    fig.subplots_adjust(left=0.075, right=0.995, top=0.78, bottom=0.30)
    out = FIGURES / "ici_cart_perturbation_network.png"
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


if __name__ == "__main__":
    print(f"wrote {main()}")
