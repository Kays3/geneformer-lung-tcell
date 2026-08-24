#!/usr/bin/env python3
"""Delete-shift vs overexpress-shift, one panel per comparison, whole-genome screen.

Same question as the targeted panel's version of this plot, at genome scale: for
every gene tested in both arms of the held-out all-gene SCLC/LUAD/normal
perturbation, does deleting it move a cell away from the goal state by roughly
the amount overexpressing it moves the cell toward that state? A concordant gene
(opposite-signed, FDR-significant, adequately detected in both arms -- the same
criterion `build_primary_report.py` uses for `concordant_primary`) sits in
quadrant II or IV.

Reads `stats/{arm}/heldout_allgene_{arm}_{comparison}.csv` under
SCLC_PERTURBATION_ROOT (default `~/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation`,
same resolution as build_primary_report.py) and writes one PNG plus the full
per-gene merged table (all genes, not just concordant hits) to `tables/`.

The response is heavy-tailed (a handful of genes carry |shift| > 0.3), so each
panel's axis limits are clipped to the 1st-99th percentile of that panel's data;
the print-out reports how many points fall outside the drawn frame.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_primary_report import ARMS, COMPARISONS, FDR, MIN_DETECTIONS, read_table  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIGURES = HERE / "figures"
OUT_TABLE = TABLES / "allgene_delete_overexpress_shift.csv"
OUT_FIGURE = FIGURES / "allgene_delete_vs_overexpress_shift.png"

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"


def merge_comparison(comparison: str) -> pd.DataFrame | None:
    delete, _ = read_table("delete", comparison)
    over, _ = read_table("overexpress", comparison)
    if delete is None or over is None:
        return None
    d = delete[["Gene_name", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(
        columns={"Shift_to_goal_end": "delete_shift", "Goal_end_FDR": "delete_fdr", "N_Detections": "delete_n"})
    o = over[["Gene_name", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(
        columns={"Shift_to_goal_end": "overexpress_shift", "Goal_end_FDR": "overexpress_fdr", "N_Detections": "overexpress_n"})
    merged = d.merge(o, on="Gene_name", how="inner")
    merged["comparison"] = comparison
    merged["concordant"] = (
        (merged["delete_fdr"] < FDR)
        & (merged["overexpress_fdr"] < FDR)
        & (merged["delete_shift"] * merged["overexpress_shift"] < 0)
        & (merged["delete_n"] >= MIN_DETECTIONS)
    )
    return merged


def build_table() -> pd.DataFrame:
    frames = [f for c in COMPARISONS if (f := merge_comparison(c)) is not None]
    table = pd.concat(frames, ignore_index=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_TABLE, index=False)
    return table


def plot(table: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(2, 3, figsize=(13, 8.5))
    for ax, comparison in zip(axes.flat, COMPARISONS):
        sub = table[table["comparison"] == comparison]
        if sub.empty:
            ax.set_visible(False)
            continue
        concordant = sub[sub["concordant"]]
        background = sub[~sub["concordant"]]

        xlo, xhi = sub["delete_shift"].quantile([0.01, 0.99])
        ylo, yhi = sub["overexpress_shift"].quantile([0.01, 0.99])
        outside = int(((sub["delete_shift"] < xlo) | (sub["delete_shift"] > xhi)
                        | (sub["overexpress_shift"] < ylo) | (sub["overexpress_shift"] > yhi)).sum())

        ax.axhline(0, color=INK2, linewidth=0.8, zorder=1)
        ax.axvline(0, color=INK2, linewidth=0.8, zorder=1)
        ax.scatter(background["delete_shift"], background["overexpress_shift"],
                   s=3, color=GREY, alpha=0.25, edgecolors="none", zorder=2,
                   label="not concordant", rasterized=True)
        ax.scatter(concordant["delete_shift"], concordant["overexpress_shift"],
                   s=5, color=BLUE, alpha=0.6, edgecolors="none", zorder=3,
                   label="concordant", rasterized=True)
        ax.set_xlim(xlo, xhi)
        ax.set_ylim(ylo, yhi)

        ax.set_title(f"{comparison}  (n={len(sub)}, concordant={len(concordant)}, "
                     f"{outside} outside frame)", fontsize=8.5, color=INK)
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
               fontsize=9, bbox_to_anchor=(0.5, 1.02), markerscale=4)
    fig.suptitle("Whole-genome screen: deletion vs overexpression shift toward the goal state\n"
                 "(axes clipped to the 1st-99th percentile per panel)",
                 fontsize=10.5, color=INK, y=1.07)
    fig.tight_layout()
    fig.savefig(OUT_FIGURE, dpi=220, bbox_inches="tight")


def main() -> None:
    table = build_table()
    plot(table)
    print(f"Wrote {OUT_TABLE.relative_to(HERE.parent.parent)}")
    print(f"Wrote {OUT_FIGURE.relative_to(HERE.parent.parent)}")
    print(table.groupby("comparison")["concordant"].agg(["size", "sum"]).rename(columns={"size": "n_genes", "sum": "n_concordant"}))


if __name__ == "__main__":
    main()
