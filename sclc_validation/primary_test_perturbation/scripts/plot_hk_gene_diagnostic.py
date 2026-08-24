#!/usr/bin/env python3
"""Figures for the HK/ubiquitous-gene enrichment diagnostic (T cells).

Two figures, same visual system as the other shift figures in this report:

  1. HK-family fraction among ISP hits vs the tested background, one bar pair
     per comparison, delete-vs-overexpress concordant hits and goal-vs-alt
     significant movers side by side.
  2. Detection-fraction gap scatter: for whole-genome HK-flagged concordant
     hits, source-state vs goal-state detection fraction (delete_n / total
     held-out cells). Near the diagonal means the gene is detected at the same
     rate in both T-cell states despite a "significant" reciprocal shift --
     the pattern that should raise suspicion. Ambient-flagged genes (already
     diagnosed as lineage-foreign, see METHODS_ambient_risk.md) are marked
     separately from HK-family genes since they are a different mechanism.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables/hk_gene_diagnostic"
FIGURES = HERE / "figures/hk_gene_diagnostic"

BLUE = "#17608f"
ORANGE = "#c2691a"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

DPI = 400
N_LABELED = 14

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
})

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]
STATE_LABEL = {"normal": "Normal", "sclc": "SCLC", "luad": "LUAD"}


def comparison_label(c: str) -> str:
    a, b = c.split("_to_")
    return f"{STATE_LABEL[a]}→{STATE_LABEL[b]}"


def plot_enrichment_bars() -> None:
    dvo = pd.read_csv(TABLES / "allgene_dvo_hk_enrichment.csv").set_index("comparison")
    ga = pd.read_csv(TABLES / "allgene_goal_alt_hk_enrichment.csv")
    ga_delete = ga[ga["arm"] == "delete"].set_index("comparison")
    ga_over = ga[ga["arm"] == "overexpress"].set_index("comparison")

    fig, ax = plt.subplots(figsize=(9, 5.5), dpi=DPI)
    x = np.arange(len(COMPARISONS))
    width = 0.2
    background = [dvo.loc[c, "pct_hk_background"] for c in COMPARISONS]

    ax.bar(x - 1.5 * width, background, width, color=GREY, label="background (all tested genes)")
    ax.bar(x - 0.5 * width, [dvo.loc[c, "pct_hk_among_hit"] for c in COMPARISONS],
           width, color=BLUE, label="delete↔overexpress concordant hits")
    ax.bar(x + 0.5 * width, [ga_delete.loc[c, "pct_hk_among_hit"] for c in COMPARISONS],
           width, color=ORANGE, label="delete-arm significant movers")
    ax.bar(x + 1.5 * width, [ga_over.loc[c, "pct_hk_among_hit"] for c in COMPARISONS],
           width, color="#8a4f9e", label="overexpress-arm significant movers")

    ax.set_xticks(x)
    ax.set_xticklabels([comparison_label(c) for c in COMPARISONS], fontsize=9, color=INK2)
    ax.set_ylabel("% of gene set flagged housekeeping/ubiquitous", fontsize=9.5, color=INK2)
    ax.tick_params(labelsize=8.5, colors=INK2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)
    ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, -0.14), ncol=2)
    ax.set_title("T cells — whole-genome screen", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    fig.suptitle("Housekeeping/ubiquitous-gene family membership: hits vs tested background",
                 fontsize=9.5, color=INK2, y=0.97, x=0.02, ha="left")
    fig.tight_layout()
    out = FIGURES / "hk_enrichment_bars.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent)}")


def plot_detection_gap() -> None:
    df = pd.read_csv(TABLES / "allgene_hk_concordant_detection_gap.csv")
    df = df.drop_duplicates(subset=["Gene_name", "comparison"])

    fig, ax = plt.subplots(figsize=(7, 7), dpi=DPI)
    ax.plot([0, 1], [0, 1], color=INK2, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
    colors = np.where(df["ambient_flag"], ORANGE, BLUE)
    ax.scatter(df["source_detect_frac"], df["goal_detect_frac"], s=26, c=colors,
               edgecolors="white", linewidths=0.3, alpha=0.85, zorder=2)

    top = df.reindex(df["detect_frac_gap"].sort_values(ascending=False).index[:N_LABELED])
    texts = [
        ax.text(row["source_detect_frac"], row["goal_detect_frac"], row["Gene_name"],
                 fontsize=7.5, style="italic", color=INK, zorder=4)
        for _, row in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))

    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_xlabel("detection fraction in source T-cell state", fontsize=9.5, color=INK2)
    ax.set_ylabel("detection fraction in goal T-cell state", fontsize=9.5, color=INK2)
    ax.tick_params(labelsize=8.5, colors=INK2, direction="out")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)

    handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=BLUE, markeredgecolor="none", label="HK family, not ambient-flagged"),
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=ORANGE, markeredgecolor="none", label="also ambient-flagged"),
    ]
    ax.legend(handles=handles, frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, -0.1))
    ax.set_title("T cells — whole-genome screen", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    fig.suptitle("HK-flagged, concordant hits: on the diagonal = equally detected in both states",
                 fontsize=9.5, color=INK2, y=1.0, x=0.02, ha="left")
    fig.tight_layout()
    out = FIGURES / "hk_detection_gap_scatter.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent)}")


def plot_targeted_panel_composition() -> None:
    df = pd.read_csv(TABLES / "targeted_panel_gene_flags.csv")
    df["flag"] = np.select(
        [df["is_hk"] & df["ambient_flag"], df["is_hk"], df["ambient_flag"]],
        ["HK + ambient", "HK family", "ambient"], default="neither",
    )
    order = ["neither", "ambient", "HK family", "HK + ambient"]
    colors = {"neither": GREY, "ambient": ORANGE, "HK family": BLUE, "HK + ambient": "#8a4f9e"}

    counts = df.groupby(["gene_source", "flag"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
    fig, ax = plt.subplots(figsize=(6, 5), dpi=DPI)
    bottom = np.zeros(len(counts))
    for flag in order:
        ax.bar(counts.index, counts[flag], bottom=bottom, color=colors[flag], label=flag, width=0.55)
        bottom += counts[flag].to_numpy()
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(["pre-registered\nimmune panel (21)", "prior top-driver\nhits (29)"], fontsize=9, color=INK2)
    ax.set_ylabel("genes", fontsize=9.5, color=INK2)
    ax.tick_params(labelsize=8.5, colors=INK2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left", bbox_to_anchor=(0.0, -0.12), ncol=2)
    ax.set_title("T cells — targeted 50-gene panel composition", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    fig.tight_layout()
    out = FIGURES / "targeted_panel_composition.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent)}")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    plot_enrichment_bars()
    plot_detection_gap()
    plot_targeted_panel_composition()


if __name__ == "__main__":
    main()
