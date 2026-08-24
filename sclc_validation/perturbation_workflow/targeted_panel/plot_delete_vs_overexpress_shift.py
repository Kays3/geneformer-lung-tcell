#!/usr/bin/env python3
"""Delete-shift vs overexpress-shift, one panel per comparison, targeted 50-gene panel.

Each point is a gene: x is its deletion shift toward the comparison's goal state,
y is its overexpression shift toward the same goal state. A gene whose edit truly
drives a cell toward/away from the goal state should show opposite-signed shifts
between the two arms (delete removes what overexpress adds) -- the "concordant"
flag already computed in the merged table, quadrants II and IV.

Reads the committed merged table, writes one PNG. No GPU, no raw pickles.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

HERE = Path(__file__).resolve().parent
MERGED = HERE / "results" / "targeted_panel_delete_overexpress_merged.csv"
OUT = HERE / "figures" / "delete_vs_overexpress_shift.png"

COMPARISONS = [
    "normal_to_sclc", "normal_to_luad",
    "sclc_to_normal", "sclc_to_luad",
    "luad_to_normal", "luad_to_sclc",
]

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"


def main() -> None:
    df = pd.read_csv(MERGED)
    OUT.parent.mkdir(parents=True, exist_ok=True)

    fig, axes = plt.subplots(2, 3, figsize=(12, 8), sharex=False, sharey=False)
    for ax, comparison in zip(axes.flat, COMPARISONS):
        sub = df[df["comparison"] == comparison]
        concordant = sub[sub["concordant"]]
        discordant = sub[~sub["concordant"]]

        ax.axhline(0, color=INK2, linewidth=0.8, zorder=1)
        ax.axvline(0, color=INK2, linewidth=0.8, zorder=1)
        ax.scatter(discordant["delete_shift"], discordant["overexpress_shift"],
                   s=22, color=GREY, edgecolors="none", zorder=2, label="not concordant")
        ax.scatter(concordant["delete_shift"], concordant["overexpress_shift"],
                   s=26, color=BLUE, edgecolors="none", zorder=3, label="concordant")

        ax.set_title(f"{comparison}  (n={len(sub)}, concordant={len(concordant)})",
                     fontsize=9, color=INK)
        ax.tick_params(labelsize=7, colors=INK2)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(INK2)

    for ax in axes[-1, :]:
        ax.set_xlabel("delete shift toward goal", fontsize=8, color=INK2)
    for ax in axes[:, 0]:
        ax.set_ylabel("overexpress shift toward goal", fontsize=8, color=INK2)

    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, frameon=False,
               fontsize=9, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle("Targeted 50-gene panel: deletion vs overexpression shift toward the goal state",
                 fontsize=11, color=INK, y=1.06)
    fig.tight_layout()
    fig.savefig(OUT, dpi=200, bbox_inches="tight")
    print(f"Wrote {OUT.relative_to(HERE.parent.parent.parent)}")


if __name__ == "__main__":
    main()
