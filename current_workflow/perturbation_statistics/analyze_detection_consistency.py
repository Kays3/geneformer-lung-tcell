#!/usr/bin/env python3
"""Rank repeatedly reported perturbations by cross-comparison cell detections."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUT = HERE / "source_tables" / "top_goal_shift_genes.csv"
OUTPUT = HERE / "detection_consistency"

COMPARISON_ORDER = [
    "LUSC → LUAD",
    "LUSC → NORMAL",
    "LUAD → LUSC",
    "LUAD → NORMAL",
    "NORMAL → LUAD",
    "NORMAL → LUSC",
]


def load_data() -> pd.DataFrame:
    frame = pd.read_csv(INPUT)
    required = {"Gene_name", "comparison_label", "N_Detections", "Shift_to_goal_end"}
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")
    return frame


def summarize(frame: pd.DataFrame) -> pd.DataFrame:
    summary = (
        frame.groupby("Gene_name", as_index=False)
        .agg(
            comparisons_present=("comparison_label", "nunique"),
            min_detected_cells=("N_Detections", "min"),
            median_detected_cells=("N_Detections", "median"),
            mean_detected_cells=("N_Detections", "mean"),
            max_detected_cells=("N_Detections", "max"),
            median_goal_shift=("Shift_to_goal_end", "median"),
        )
    )
    summary["comparison_coverage"] = summary["comparisons_present"] / len(COMPARISON_ORDER)
    summary = summary.sort_values(
        ["comparisons_present", "min_detected_cells", "median_detected_cells"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    summary.insert(0, "consistency_rank", np.arange(1, len(summary) + 1))
    return summary


def make_bubble_plot(frame: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    repeated = summary.loc[summary["comparisons_present"] >= 2].copy()
    genes = repeated["Gene_name"].tolist()
    plot = frame.loc[frame["Gene_name"].isin(genes)].copy()
    plot["comparison_label"] = pd.Categorical(
        plot["comparison_label"], categories=COMPARISON_ORDER, ordered=True
    )
    plot["Gene_name"] = pd.Categorical(plot["Gene_name"], categories=genes[::-1], ordered=True)

    fig, ax = plt.subplots(figsize=(10.8, max(4.8, 0.43 * len(genes) + 2.0)))
    sizes = 45 + 455 * np.sqrt(plot["N_Detections"] / plot["N_Detections"].max())
    points = ax.scatter(
        plot["comparison_label"],
        plot["Gene_name"],
        s=sizes,
        c=plot["Shift_to_goal_end"],
        cmap="cividis",
        edgecolor="#263238",
        linewidth=0.7,
        alpha=0.9,
    )
    ax.set_title("Repeated perturbations across directional comparisons", loc="left", weight="bold")
    ax.set_xlabel("Directional comparison")
    ax.set_ylabel("Perturbed gene (ordered by consistency rank)")
    ax.grid(axis="both", color="#d9dde2", linewidth=0.7, alpha=0.65)
    ax.set_axisbelow(True)
    ax.tick_params(axis="x", rotation=28)
    colorbar = fig.colorbar(points, ax=ax, pad=0.02)
    colorbar.set_label("Shift toward goal endpoint")

    legend_counts = sorted({int(plot["N_Detections"].min()), int(plot["N_Detections"].median()), int(plot["N_Detections"].max())})
    handles = [
        ax.scatter([], [], s=45 + 455 * np.sqrt(v / plot["N_Detections"].max()),
                   facecolor="#7895b2", edgecolor="#263238", linewidth=0.7)
        for v in legend_counts
    ]
    ax.legend(handles, [f"{v} cells" for v in legend_counts], title="Bubble size",
              frameon=False, loc="upper left", bbox_to_anchor=(1.12, 1.0))
    fig.text(
        0.01,
        0.01,
        "Scope: published top 15 qualified toward-goal shifts per comparison; absence is not evidence of zero detections.",
        fontsize=8.5,
        color="#4f5963",
    )
    fig.tight_layout(rect=(0, 0.05, 0.9, 1))
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def make_consistency_plot(summary: pd.DataFrame, path: Path) -> None:
    repeated = summary.loc[summary["comparisons_present"] >= 2].copy()
    fig, ax = plt.subplots(figsize=(9.2, 6.4))
    sizes = 120 * repeated["comparisons_present"] ** 2
    points = ax.scatter(
        repeated["median_detected_cells"],
        repeated["min_detected_cells"],
        s=sizes,
        c=repeated["median_goal_shift"],
        cmap="cividis",
        edgecolor="#263238",
        linewidth=0.8,
        alpha=0.9,
    )
    for row in repeated.itertuples(index=False):
        ax.annotate(row.Gene_name, (row.median_detected_cells, row.min_detected_cells),
                    xytext=(5, 4), textcoords="offset points", fontsize=8.5)
    ax.set_title("Detection support among recurring perturbations", loc="left", weight="bold")
    ax.set_xlabel("Median detected cells across appearances")
    ax.set_ylabel("Minimum detected cells across appearances")
    ax.grid(color="#d9dde2", linewidth=0.7, alpha=0.65)
    ax.set_axisbelow(True)
    colorbar = fig.colorbar(points, ax=ax, pad=0.02)
    colorbar.set_label("Median shift toward goal endpoint")
    handles = [
        ax.scatter([], [], s=120 * n**2, facecolor="#7895b2", edgecolor="#263238", linewidth=0.7)
        for n in (2, 3)
    ]
    ax.legend(handles, ["2 comparisons", "3 comparisons"], title="Bubble size",
              frameon=False, loc="lower right")
    fig.text(0.01, 0.01, "Upper-right indicates stronger detected-cell support in both typical and weakest appearances.",
             fontsize=8.5, color="#4f5963")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_report(frame: pd.DataFrame, summary: pd.DataFrame, path: Path) -> None:
    repeated = summary.loc[summary["comparisons_present"] >= 2].copy()
    top = repeated.head(10)
    rows = []
    for row in top.itertuples(index=False):
        rows.append(
            f"| {row.consistency_rank} | {row.Gene_name} | {row.comparisons_present}/6 | "
            f"{row.min_detected_cells:.0f} | {row.median_detected_cells:.1f} | "
            f"{row.max_detected_cells:.0f} | {row.median_goal_shift:.4f} |"
        )

    leaders = ", ".join(top.head(3)["Gene_name"])
    text = f"""# Cross-comparison perturbation detection consistency

## Technical summary

Within the published top-15 qualified shifts for each of six directional comparisons,
**{len(repeated)} genes recur in at least two comparisons**. The highest-coverage genes are
**{leaders}**. Ranking prioritizes the number of comparisons represented, then the minimum
detected-cell count (a conservative consistency criterion), then the median count.

![Bubble plot](bubble_plot.png)

Bubble area represents `N_Detections`; color represents the positive shift toward the goal
endpoint. Empty cells mean the gene was not in that comparison's published top 15, not that
the gene was untested or had zero detections.

![Consistency summary](consistency_summary.png)

This second view compares typical support (median detected cells) with worst-observed support
(minimum detected cells). Bubble area is comparison coverage; upper-right genes have stronger
support by both measures.

## Highest recurring detection counts

| Rank | Gene | Comparisons | Minimum cells | Median cells | Maximum cells | Median goal shift |
|---:|---|---:|---:|---:|---:|---:|
{chr(10).join(rows)}

## Scope and metric definitions

- Source: `source_tables/top_goal_shift_genes.csv` (90 rows: top 15 qualified shifts in each comparison).
- Qualified shifts were previously filtered to `Goal_end_FDR < 0.05`, positive goal shift, and
  `N_Detections >= 25`.
- `comparisons_present` counts directional comparisons in which a gene appears in the published top 15.
- `minimum cells` is the lowest `N_Detections` among those appearances and is used to reward stable support.

## Limitations and next step

This is a consistency check of the **published top-ranked subset**, not the complete perturbation
matrix. A definitive six-comparison consistency analysis requires the six full
`heldout_allgene_<comparison>.csv` files. Once available, the same script should be pointed to those
tables so every tested gene can be distinguished from a gene absent only because of top-15 truncation.
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    frame = load_data()
    summary = summarize(frame)
    summary.to_csv(OUTPUT / "detection_consistency_ranking.csv", index=False)
    make_bubble_plot(frame, summary, OUTPUT / "bubble_plot.png")
    make_consistency_plot(summary, OUTPUT / "consistency_summary.png")
    write_report(frame, summary, OUTPUT / "REPORT.md")
    print(f"Wrote {OUTPUT / 'REPORT.md'}")


if __name__ == "__main__":
    main()
