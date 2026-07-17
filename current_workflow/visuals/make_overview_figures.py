#!/usr/bin/env python3
"""Reproduce the compact figures used by the repository README."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = Path(__file__).resolve().parent
BLUE = "#2867B2"
ORANGE = "#D9822B"
INK = "#25313C"
GRID = "#D9E0E6"


def style() -> None:
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.titleweight": "bold",
            "axes.labelcolor": INK,
            "text.color": INK,
            "xtick.color": INK,
            "ytick.color": INK,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def performance_context() -> None:
    prior = pd.read_csv(
        ROOT / "archive/prior_nsclc_workflow/tables/classifier_metrics.csv"
    )
    manifest = json.loads((ROOT / "current_workflow/experiment_manifest.json").read_text())
    final = manifest["fine_tuning"]
    data = pd.DataFrame(
        [
            {
                "model": "Prior: cell type\n9 classes",
                "accuracy": prior.loc[prior.stage == "Stage 1", "accuracy"].iloc[0],
                "macro_f1": prior.loc[prior.stage == "Stage 1", "macro_f1"].iloc[0],
            },
            {
                "model": "Prior: disease\n4 classes",
                "accuracy": prior.loc[prior.stage == "Stage 2", "accuracy"].iloc[0],
                "macro_f1": prior.loc[prior.stage == "Stage 2", "macro_f1"].iloc[0],
            },
            {
                "model": "Current: T-cell disease\n3 classes",
                "accuracy": final["test_accuracy"],
                "macro_f1": final["test_macro_f1"],
            },
        ]
    )
    data.to_csv(OUTPUT / "model_performance_context.csv", index=False)

    y = np.arange(len(data))
    height = 0.28
    fig, ax = plt.subplots(figsize=(9.4, 4.4))
    ax.barh(y - height / 2, data.accuracy, height, color=BLUE, label="Accuracy")
    ax.barh(y + height / 2, data.macro_f1, height, color=ORANGE, label="Macro F1")
    for row, (_, values) in enumerate(data.iterrows()):
        ax.text(values.accuracy + 0.012, row - height / 2, f"{values.accuracy:.3f}", va="center", weight="bold")
        ax.text(values.macro_f1 + 0.012, row + height / 2, f"{values.macro_f1:.3f}", va="center", weight="bold")
    ax.set_yticks(y, data.model)
    ax.invert_yaxis()
    ax.set_xlim(0, 1.04)
    ax.set_xlabel("Held-out score")
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.set_title("Classifier performance context", loc="left", pad=26)
    ax.text(
        0,
        1.04,
        "Different cohorts and tasks; values are contextual, not a model ranking.",
        transform=ax.transAxes,
        fontsize=9,
        color="#596773",
    )
    ax.legend(frameon=False, ncol=1, loc="center right")
    fig.tight_layout()
    fig.savefig(OUTPUT / "model_performance_context.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


def confusion_matrix() -> None:
    matrix_df = pd.read_csv(
        ROOT / "current_workflow/source_tables/final_test_confusion_matrix.csv", index_col=0
    )
    matrix = matrix_df.to_numpy()
    row_percent = matrix / matrix.sum(axis=1, keepdims=True) * 100
    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    image = ax.imshow(row_percent, cmap="Blues", vmin=0, vmax=100)
    for row in range(matrix.shape[0]):
        for column in range(matrix.shape[1]):
            color = "white" if row_percent[row, column] > 52 else INK
            ax.text(
                column,
                row,
                f"{matrix[row, column]:,}\n{row_percent[row, column]:.1f}%",
                ha="center",
                va="center",
                color=color,
                weight="bold",
            )
    ax.set_xticks(range(3), matrix_df.columns)
    ax.set_yticks(range(3), matrix_df.index)
    ax.set_xlabel("Predicted disease")
    ax.set_ylabel("Actual disease")
    ax.set_title("Current T-cell classifier: held-out confusion matrix", loc="left", pad=24)
    ax.text(
        0,
        1.03,
        "Labels show cell count and percentage within each actual class (n=3,379).",
        transform=ax.transAxes,
        fontsize=9,
        color="#596773",
    )
    colorbar = fig.colorbar(image, ax=ax, fraction=0.046, pad=0.04)
    colorbar.set_label("Row percentage")
    for spine in ax.spines.values():
        spine.set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT / "final_tcell_confusion_matrix.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    style()
    performance_context()
    confusion_matrix()
