#!/usr/bin/env python3
"""Delete-shift vs overexpress-shift, one high-resolution PNG per comparison,
targeted 50-gene panel.

Each point is a gene: x is its deletion shift toward the comparison's goal
state, y is its overexpression shift toward the same goal state. A gene whose
edit truly drives a cell toward/away from the goal state should show
opposite-signed shifts between the two arms (delete removes what overexpress
adds) -- the "concordant" flag already computed in the merged table,
quadrants II and IV. Marker area scales with the weaker of the two arms'
N_Detections (sqrt scale, shared across all six comparisons) and color marks
concordance. The most-displaced concordant genes are labeled.

Same visual style as `plot_goal_vs_alt_shift.py` in this folder, so Part 1 and
Part 2 of the shift report read as one system.

Reads the committed merged table, writes six PNGs. No GPU, no raw pickles.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
MERGED = HERE / "results" / "targeted_panel_delete_overexpress_merged.csv"
FIGURES = HERE / "figures" / "delete_vs_overexpress_shift"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

DPI = 400
N_LABELED = 10
STATE_LABEL = {"normal": "Normal T cells", "sclc": "SCLC T cells", "luad": "LUAD T cells"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
})


def size_scale(n: pd.Series, n_min: float, n_max: float, s_min: float = 10.0, s_max: float = 320.0) -> np.ndarray:
    n_clip = n.clip(lower=n_min, upper=n_max)
    frac = (np.sqrt(n_clip) - np.sqrt(n_min)) / (np.sqrt(n_max) - np.sqrt(n_min))
    return s_min + frac * (s_max - s_min)


def plot_one(comparison: str, df: pd.DataFrame, n_min: float, n_max: float, out: Path) -> None:
    weakest_n = df[["delete_n", "overexpress_n"]].min(axis=1)
    sizes = size_scale(weakest_n, n_min, n_max)
    colors = np.where(df["concordant"], BLUE, GREY)

    lo = min(df["delete_shift"].min(), df["overexpress_shift"].min())
    hi = max(df["delete_shift"].max(), df["overexpress_shift"].max())
    pad = 0.08 * (hi - lo)

    fig, ax = plt.subplots(figsize=(6.2, 6.2), dpi=DPI)
    ax.axhline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.axvline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.scatter(df["delete_shift"], df["overexpress_shift"], s=sizes, c=colors,
               edgecolors="white", linewidths=0.4, alpha=0.9, zorder=2)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)

    concordant = df[df["concordant"]].copy()
    concordant["distance"] = np.hypot(concordant["delete_shift"], concordant["overexpress_shift"])
    top = concordant.nlargest(N_LABELED, "distance")
    texts = [
        ax.text(row["delete_shift"], row["overexpress_shift"], row["Gene_name"],
                 fontsize=7.5, style="italic", color=INK, zorder=4)
        for _, row in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))

    source, goal = comparison.split("_to_")
    ax.set_title(f"{STATE_LABEL[source]} → {STATE_LABEL[goal]}", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    ax.set_xlabel("delete shift toward goal", fontsize=9.5, color=INK2)
    ax.set_ylabel("overexpress shift toward goal", fontsize=9.5, color=INK2)
    ax.tick_params(labelsize=8.5, colors=INK2, direction="out")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)

    color_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=BLUE, markeredgecolor="none", label="concordant"),
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=GREY, markeredgecolor="none", label="not concordant"),
    ]
    size_refs = sorted(set(int(round(v)) for v in np.geomspace(max(1, n_min), n_max, 4)))
    size_handles = [
        Line2D([], [], marker="o", linestyle="none", markerfacecolor="none",
               markeredgecolor=INK2, markeredgewidth=0.8,
               markersize=np.sqrt(size_scale(pd.Series([v]), n_min, n_max)[0]), label=f"N={v}")
        for v in size_refs
    ]
    leg1 = ax.legend(handles=color_handles, loc="upper left", bbox_to_anchor=(0.0, -0.11),
                      frameon=False, fontsize=8.5, ncol=1, handletextpad=0.4)
    ax.legend(handles=size_handles, loc="upper right", bbox_to_anchor=(1.0, -0.11),
              frameon=False, fontsize=8.5, ncol=1, title="min(N detections)", title_fontsize=8.5,
              handletextpad=0.8, labelspacing=1.1)
    ax.add_artist(leg1)
    fig.suptitle("T cells — targeted 50-gene panel", fontsize=9.5, color=INK2, y=1.0, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent.parent)}")


def main() -> None:
    df = pd.read_csv(MERGED)
    FIGURES.mkdir(parents=True, exist_ok=True)
    weakest_n = df[["delete_n", "overexpress_n"]].min(axis=1)
    n_min, n_max = max(1, weakest_n.min()), weakest_n.max()
    for comparison in COMPARISONS:
        plot_one(comparison, df[df["comparison"] == comparison], n_min, n_max, FIGURES / f"{comparison}.png")


if __name__ == "__main__":
    main()
