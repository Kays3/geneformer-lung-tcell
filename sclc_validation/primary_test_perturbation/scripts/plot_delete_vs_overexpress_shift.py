#!/usr/bin/env python3
"""Delete-shift vs overexpress-shift, one high-resolution PNG per comparison,
whole-genome screen.

Same question as the targeted panel's version of this plot, at genome scale: for
every gene tested in both arms of the held-out all-gene SCLC/LUAD/normal
perturbation, does deleting it move a cell away from the goal state by roughly
the amount overexpressing it moves the cell toward that state? A concordant gene
(opposite-signed, FDR-significant, adequately detected in both arms -- the same
criterion `build_primary_report.py` uses for `concordant_primary`) sits in
quadrant II or IV. Marker area scales with the weaker of the two arms'
N_Detections; the most-displaced concordant genes are labeled.

Same visual style as `plot_goal_vs_alt_shift.py` in this folder, so Part 1 and
Part 2 of the shift report read as one system.

Reads `stats/{arm}/heldout_allgene_{arm}_{comparison}.csv` under
SCLC_PERTURBATION_ROOT (default `~/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation`,
same resolution as build_primary_report.py) and writes six PNGs plus the full
per-gene merged table (all genes, not just concordant hits) to `tables/`.

The response is heavy-tailed (a handful of genes carry |shift| > 0.3), so each
panel's axis limits are clipped to the 1st-99th percentile of that panel's data;
the panel title reports how many points fall outside the drawn frame.
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_primary_report import ARMS, COMPARISONS, FDR, MIN_DETECTIONS, read_table  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIGURES = HERE / "figures" / "delete_vs_overexpress_shift"
OUT_TABLE = TABLES / "allgene_delete_overexpress_shift.csv"

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

DPI = 400
N_LABELED = 12
STATE_LABEL = {"normal": "Normal T cells", "sclc": "SCLC T cells", "luad": "LUAD T cells"}

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
    "xtick.major.width": 0.75,
    "ytick.major.width": 0.75,
})


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


def size_scale(n: pd.Series, n_min: float, n_max: float, s_min: float = 2.0, s_max: float = 90.0) -> np.ndarray:
    n_clip = n.clip(lower=n_min, upper=n_max)
    frac = (np.sqrt(n_clip) - np.sqrt(n_min)) / (np.sqrt(n_max) - np.sqrt(n_min))
    return s_min + frac * (s_max - s_min)


def plot_one(comparison: str, df: pd.DataFrame, n_min: float, n_max: float, out: Path) -> None:
    weakest_n = df[["delete_n", "overexpress_n"]].min(axis=1)
    sizes = size_scale(weakest_n, n_min, n_max)
    concordant = df["concordant"]

    xlo, xhi = df["delete_shift"].quantile([0.01, 0.99])
    ylo, yhi = df["overexpress_shift"].quantile([0.01, 0.99])
    outside = int(((df["delete_shift"] < xlo) | (df["delete_shift"] > xhi)
                    | (df["overexpress_shift"] < ylo) | (df["overexpress_shift"] > yhi)).sum())

    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=DPI)
    ax.axhline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.axvline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.scatter(df.loc[~concordant, "delete_shift"], df.loc[~concordant, "overexpress_shift"],
               s=sizes[~concordant], c=GREY, edgecolors="none", alpha=0.35, zorder=2, rasterized=True)
    ax.scatter(df.loc[concordant, "delete_shift"], df.loc[concordant, "overexpress_shift"],
               s=sizes[concordant], c=BLUE, edgecolors="white", linewidths=0.25, alpha=0.6, zorder=3, rasterized=True)
    ax.set_xlim(xlo, xhi)
    ax.set_ylim(ylo, yhi)

    # Rank only among points inside the drawn frame.
    conc_df = df[concordant].copy()
    conc_df = conc_df[conc_df["delete_shift"].between(xlo, xhi) & conc_df["overexpress_shift"].between(ylo, yhi)]
    conc_df["distance"] = np.hypot(conc_df["delete_shift"], conc_df["overexpress_shift"])
    top = conc_df.nlargest(N_LABELED, "distance")
    texts = [
        ax.text(row["delete_shift"], row["overexpress_shift"], row["Gene_name"],
                 fontsize=7.5, style="italic", color=INK, zorder=4)
        for _, row in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))

    source, goal = comparison.split("_to_")
    ax.set_title(f"{STATE_LABEL[source]} → {STATE_LABEL[goal]}", fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    ax.text(1.0, 1.02, f"n={len(df)}, concordant={int(concordant.sum())}, {outside} outside frame",
            transform=ax.transAxes, fontsize=7.5, color=INK2, ha="right", va="bottom")
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
    fig.suptitle("T cells — whole-genome screen (axes clipped to 1st-99th pct)", fontsize=9.5, color=INK2, y=1.0, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent)}")


def plot(table: pd.DataFrame) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    weakest_n = table[["delete_n", "overexpress_n"]].min(axis=1)
    n_min, n_max = max(1, weakest_n.min()), weakest_n.max()
    for comparison in COMPARISONS:
        sub = table[table["comparison"] == comparison]
        if sub.empty:
            continue
        plot_one(comparison, sub, n_min, n_max, FIGURES / f"{comparison}.png")


def main() -> None:
    table = build_table()
    plot(table)
    print(f"Wrote {OUT_TABLE.relative_to(HERE.parent.parent)}")
    print(table.groupby("comparison")["concordant"].agg(["size", "sum"]).rename(columns={"size": "n_genes", "sum": "n_concordant"}))


if __name__ == "__main__":
    main()
