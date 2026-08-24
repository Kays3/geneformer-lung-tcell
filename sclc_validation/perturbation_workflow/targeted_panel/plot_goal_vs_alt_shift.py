#!/usr/bin/env python3
"""Shift toward the goal state vs shift toward the alternate state, targeted panel.

Every perturbed cell has two off-target-free readouts from the *same* edit: how
far it moved toward the comparison's goal state, and how far it moved toward the
third, unnamed "alt" state. A gene whose effect is specific to the goal sits
near the x-axis (large goal shift, ~0 alt shift); a gene sitting near the y=x
line is moving the cell toward both other states about equally -- not
specific to the goal at all.

One high-resolution PNG per arm (delete, overexpress) per comparison -- 12
files total. Marker area scales with N_Detections (sqrt scale, shared across
all six comparisons within an arm so the size legend is comparable) and color
marks whether the goal shift is FDR-significant (`Sig`, already computed
upstream as Goal_end_FDR < 0.05).

Reads the six per-comparison stat tables already committed under
`results/{arm}/`. No GPU, no raw pickles.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

N_LABELED = 10

HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]
ARMS = ("delete", "overexpress")

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

DPI = 400

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
})


def alt_column(df: pd.DataFrame) -> tuple[str, str]:
    col = next(c for c in df.columns if c.startswith("Shift_to_alt_end_"))
    return col, col.removeprefix("Shift_to_alt_end_")


def load_arm(arm: str) -> dict[str, pd.DataFrame]:
    return {c: pd.read_csv(RESULTS / arm / f"targeted_{arm}_{c}.csv") for c in COMPARISONS}


def size_scale(n: pd.Series, n_min: float, n_max: float, s_min: float = 10.0, s_max: float = 320.0) -> np.ndarray:
    n_clip = n.clip(lower=n_min, upper=n_max)
    frac = (np.sqrt(n_clip) - np.sqrt(n_min)) / (np.sqrt(n_max) - np.sqrt(n_min))
    return s_min + frac * (s_max - s_min)


def plot_one(arm: str, comparison: str, df: pd.DataFrame, n_min: float, n_max: float, out: Path) -> None:
    alt_col, alt_name = alt_column(df)
    sizes = size_scale(df["N_Detections"], n_min, n_max)
    colors = np.where(df["Sig"] == 1, BLUE, GREY)

    lo = min(df["Shift_to_goal_end"].min(), df[alt_col].min())
    hi = max(df["Shift_to_goal_end"].max(), df[alt_col].max())
    pad = 0.08 * (hi - lo)

    fig, ax = plt.subplots(figsize=(6.2, 6.2), dpi=DPI)
    ax.plot([lo - pad, hi + pad], [lo - pad, hi + pad], color=INK2, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
    ax.axhline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.axvline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.scatter(df["Shift_to_goal_end"], df[alt_col], s=sizes, c=colors,
               edgecolors="white", linewidths=0.4, alpha=0.9, zorder=2)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(lo - pad, hi + pad)

    significant = df[df["Sig"] == 1].copy()
    significant["distance"] = np.hypot(significant["Shift_to_goal_end"], significant[alt_col])
    top = significant.nlargest(N_LABELED, "distance")
    texts = [
        ax.text(row["Shift_to_goal_end"], row[alt_col], row["Gene_name"],
                 fontsize=7.5, style="italic", color=INK, zorder=4)
        for _, row in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))

    source, goal = comparison.split("_to_")
    ax.set_title(f"{source} → {goal}  (alt: {alt_name})", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    ax.set_xlabel("shift toward goal state", fontsize=9.5, color=INK2)
    ax.set_ylabel("shift toward alt state", fontsize=9.5, color=INK2)
    ax.tick_params(labelsize=8.5, colors=INK2, direction="out")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)

    color_handles = [
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=BLUE, markeredgecolor="none", label="goal shift FDR < 0.05"),
        Line2D([], [], marker="o", linestyle="none", markersize=6, markerfacecolor=GREY, markeredgecolor="none", label="not significant"),
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
              frameon=False, fontsize=8.5, ncol=1, title="N detections", title_fontsize=8.5,
              handletextpad=0.8, labelspacing=1.1)
    ax.add_artist(leg1)
    fig.suptitle(f"Targeted 50-gene panel — {arm}", fontsize=9.5, color=INK2, y=1.0, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent.parent)}")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for arm in ARMS:
        frames = load_arm(arm)
        all_n = pd.concat([f["N_Detections"] for f in frames.values()])
        n_min, n_max = max(1, all_n.min()), all_n.max()
        for comparison, df in frames.items():
            plot_one(arm, comparison, df, n_min, n_max, FIGURES / f"{arm}_{comparison}_goal_vs_alt_shift.png")


if __name__ == "__main__":
    main()
