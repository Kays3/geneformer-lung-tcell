#!/usr/bin/env python3
"""Shift toward the goal state vs shift toward the alternate state, whole genome.

Genome-scale counterpart of the targeted panel's
`plot_goal_vs_alt_shift.py`: for every gene tested, how far the edit moves a
cell toward the comparison's goal state vs toward the third, "alt" state from
the same perturbation. Marker area scales with N_Detections (sqrt scale,
shared across all six comparisons within an arm) and color marks whether the
goal shift is FDR-significant (`Sig`). The top genes by movement magnitude
among significant hits are labeled.

One high-resolution PNG per arm per comparison -- 12 files total. Reads
`stats/{arm}/heldout_allgene_{arm}_{comparison}.csv` under
SCLC_PERTURBATION_ROOT (default `~/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation`,
same resolution as build_primary_report.py), and writes the full per-gene
table (all genes, both arms, all comparisons) to `tables/`.

The response is heavy-tailed, so each panel's axis limits are clipped to the
1st-99th percentile of that panel's own two columns; the title reports how
many points fall outside the drawn frame.
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
from build_primary_report import ARMS, COMPARISONS, STATS, path_for  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
FIGURES = HERE / "figures"
TABLES = HERE / "tables"
OUT_TABLE = TABLES / "allgene_goal_vs_alt_shift.csv"

BLUE = "#17608f"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"

DPI = 400
N_LABELED = 12
STATE_LABEL = {"normal": "Normal T cells", "sclc": "SCLC T cells", "luad": "LUAD T cells"}
ALT_LABEL = {"normal": "Normal T cells", "lung adenocarcinoma": "LUAD T cells", "small cell lung carcinoma": "SCLC T cells"}

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
    frames = {}
    for comparison in COMPARISONS:
        path = path_for(arm, comparison)
        if not path.exists():
            continue
        frames[comparison] = pd.read_csv(path).drop(columns=["Unnamed: 0"], errors="ignore")
    return frames


def size_scale(n: pd.Series, n_min: float, n_max: float, s_min: float = 2.0, s_max: float = 90.0) -> np.ndarray:
    n_clip = n.clip(lower=n_min, upper=n_max)
    frac = (np.sqrt(n_clip) - np.sqrt(n_min)) / (np.sqrt(n_max) - np.sqrt(n_min))
    return s_min + frac * (s_max - s_min)


def plot_one(arm: str, comparison: str, df: pd.DataFrame, n_min: float, n_max: float, out: Path) -> None:
    alt_col, alt_name = alt_column(df)
    sizes = size_scale(df["N_Detections"], n_min, n_max)
    significant = df["Sig"] == 1

    xlo, xhi = df["Shift_to_goal_end"].quantile([0.01, 0.99])
    ylo, yhi = df[alt_col].quantile([0.01, 0.99])
    lo, hi = min(xlo, ylo), max(xhi, yhi)
    outside = int(((df["Shift_to_goal_end"] < lo) | (df["Shift_to_goal_end"] > hi)
                    | (df[alt_col] < lo) | (df[alt_col] > hi)).sum())

    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=DPI)
    ax.plot([lo, hi], [lo, hi], color=INK2, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
    ax.axhline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.axvline(0, color=INK2, linewidth=0.7, zorder=1)
    ax.scatter(df.loc[~significant, "Shift_to_goal_end"], df.loc[~significant, alt_col],
               s=sizes[~significant], c=GREY, edgecolors="none", alpha=0.35, zorder=2, rasterized=True)
    ax.scatter(df.loc[significant, "Shift_to_goal_end"], df.loc[significant, alt_col],
               s=sizes[significant], c=BLUE, edgecolors="white", linewidths=0.25, alpha=0.6, zorder=3, rasterized=True)
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)

    # Rank only among points inside the drawn frame -- otherwise the biggest
    # movers are clipped outliers and nothing visible gets labeled.
    sig_df = df[significant].copy()
    sig_df = sig_df[sig_df["Shift_to_goal_end"].between(lo, hi) & sig_df[alt_col].between(lo, hi)]
    sig_df["distance"] = np.hypot(sig_df["Shift_to_goal_end"], sig_df[alt_col])
    top = sig_df.nlargest(N_LABELED, "distance")
    texts = [
        ax.text(row["Shift_to_goal_end"], row[alt_col], row["Gene_name"],
                 fontsize=7.5, style="italic", color=INK, zorder=4)
        for _, row in top.iterrows()
    ]
    if texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.5))

    source, goal = comparison.split("_to_")
    ax.set_title(f"{STATE_LABEL[source]} → {STATE_LABEL[goal]}  (alt: {ALT_LABEL.get(alt_name.lower(), alt_name)})",
                 fontsize=11, color=INK, loc="left", pad=10, fontweight="bold")
    ax.text(1.0, 1.02, f"n={len(df)}, {outside} outside frame", transform=ax.transAxes,
            fontsize=7.5, color=INK2, ha="right", va="bottom")
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
    fig.suptitle(f"T cells — whole-genome screen — {arm} (axes clipped to 1st-99th pct)",
                 fontsize=9.5, color=INK2, y=1.0, x=0.02, ha="left")
    fig.tight_layout()
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(HERE.parent.parent)}")


def build_table(all_frames: dict[str, dict[str, pd.DataFrame]]) -> None:
    rows = []
    for arm, frames in all_frames.items():
        for comparison, df in frames.items():
            alt_col, alt_name = alt_column(df)
            rows.append(pd.DataFrame({
                "arm": arm,
                "comparison": comparison,
                "Gene_name": df["Gene_name"],
                "alt_state": alt_name,
                "shift_to_goal": df["Shift_to_goal_end"],
                "shift_to_alt": df[alt_col],
                "n_detections": df["N_Detections"],
                "sig": df["Sig"],
            }))
    table = pd.concat(rows, ignore_index=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    table.to_csv(OUT_TABLE, index=False)
    print(f"Wrote {OUT_TABLE.relative_to(HERE.parent.parent)}")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    all_frames = {}
    for arm in ARMS:
        frames = load_arm(arm)
        if not frames:
            print(f"No data found for arm={arm} under {STATS}")
            continue
        all_frames[arm] = frames
        out_dir = FIGURES / "goal_vs_alt_shift" / arm
        out_dir.mkdir(parents=True, exist_ok=True)
        all_n = pd.concat([f["N_Detections"] for f in frames.values()])
        n_min, n_max = max(1, all_n.min()), all_n.max()
        for comparison, df in frames.items():
            plot_one(arm, comparison, df, n_min, n_max, out_dir / f"{comparison}.png")
    if all_frames:
        build_table(all_frames)


if __name__ == "__main__":
    main()
