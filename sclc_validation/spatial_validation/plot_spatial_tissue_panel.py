#!/usr/bin/env python3
"""Render GSE263196 tissue images with spatial validation score overlays."""

from __future__ import annotations

import gzip
import io
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.gridspec import GridSpec
import numpy as np
import pandas as pd
from PIL import Image


HERE = Path(__file__).resolve().parent
RAW_DIR = HERE.parent / "audit" / "source_metadata" / "GSE263196_RAW"
RESULTS_DIR = HERE / "results"
FIGURE_DIR = HERE / "figures"

SAMPLES = [
    ("GSM8187469", "SCLC3"),
    ("GSM8187470", "SCLC4"),
    ("GSM8187471", "SCLC8"),
    ("GSM8187472", "SCLC9"),
    ("GSM8187473", "SCLC12"),
]

INK = "#20262E"
MUTED = "#66717E"
GRID = "#D7DCE2"
BLUE = "#2F6690"
ORANGE = "#D97745"


def load_tissue_bundle(gsm: str, label: str, scores: pd.DataFrame) -> tuple[Image.Image, pd.DataFrame]:
    prefix = f"{gsm}_{label}"
    image_path = RAW_DIR / f"{prefix}_tissue_lowres_image.png.gz"
    positions_path = RAW_DIR / f"{prefix}_tissue_positions_list.csv.gz"
    scale_path = RAW_DIR / f"{prefix}_scalefactors_json.json.gz"
    missing = [path for path in (image_path, positions_path, scale_path) if not path.exists()]
    if missing:
        missing_text = "\n".join(str(path) for path in missing)
        raise FileNotFoundError(
            "Missing GEO tissue assets. Extract the GSE263196_RAW.tar image, "
            f"position, and scale-factor files into {RAW_DIR}:\n{missing_text}"
        )

    with gzip.open(image_path, "rb") as handle:
        image = Image.open(io.BytesIO(handle.read())).convert("RGB")
    with gzip.open(scale_path, "rt") as handle:
        scale = json.load(handle)["tissue_lowres_scalef"]

    positions = pd.read_csv(
        positions_path,
        header=None,
        names=[
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ],
    )
    sample_scores = scores.loc[scores["sample_gsm"].eq(gsm)].copy()
    spatial = sample_scores.merge(positions, on="barcode", how="left", validate="one_to_one")
    if spatial["pxl_row_in_fullres"].isna().any():
        raise ValueError(f"Missing spatial coordinates after barcode join for {gsm}")
    if not spatial["in_tissue"].eq(1).all():
        raise ValueError(f"Non-tissue spots found in scored output for {gsm}")

    spatial["image_x"] = spatial["pxl_col_in_fullres"] * scale
    spatial["image_y"] = spatial["pxl_row_in_fullres"] * scale
    return image, spatial


def format_p_value(value: float) -> str:
    if value < 1e-99:
        return "p<10⁻⁹⁹"
    if value < 0.001:
        exponent = int(np.floor(np.log10(value)))
        coefficient = value / (10**exponent)
        superscript = str(exponent).translate(str.maketrans("0123456789-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁻"))
        return f"p={coefficient:.1f}×10{superscript}"
    return f"p={value:.3f}"


def draw_tissue_axis(
    axis: plt.Axes,
    image: Image.Image,
    spatial: pd.DataFrame,
    score_column: str,
    normalization: Normalize,
    colormap: str,
) -> None:
    axis.imshow(image)
    axis.scatter(
        spatial["image_x"],
        spatial["image_y"],
        c=spatial[score_column],
        cmap=colormap,
        norm=normalization,
        s=8,
        alpha=0.76,
        linewidths=0,
        rasterized=True,
    )
    axis.set_xlim(0, image.width)
    axis.set_ylim(image.height, 0)
    axis.set_aspect("equal")
    axis.set_axis_off()


def draw_forest(axis: plt.Axes, per_sample: pd.DataFrame, pooled: pd.Series) -> None:
    rows = per_sample.to_dict(orient="records")
    labels = [
        f"{row['sample_label']}  (n={int(row['n_spots']):,})"
        for row in rows
    ] + [f"Pooled  ({len(rows)} samples; n={int(per_sample['n_spots'].sum()):,})"]
    estimates = [row["rho"] for row in rows] + [pooled["rho_pooled"]]
    lows = [row["ci_low"] for row in rows] + [pooled["ci_low"]]
    highs = [row["ci_high"] for row in rows] + [pooled["ci_high"]]
    p_values = [row["p_value"] for row in rows] + [pooled["p_value"]]
    y_positions = np.arange(len(labels))

    axis.axvline(0, color=MUTED, linewidth=1, linestyle="--", zorder=1)
    for index, (estimate, low, high) in enumerate(zip(estimates, lows, highs)):
        color = ORANGE if index == len(labels) - 1 else BLUE
        marker = "D" if index == len(labels) - 1 else "o"
        axis.plot([low, high], [index, index], color=color, linewidth=2.2, zorder=2)
        axis.scatter(
            estimate,
            index,
            color=color,
            marker=marker,
            s=62,
            edgecolor="white",
            linewidth=0.8,
            zorder=3,
        )
        axis.text(
            1.02,
            index,
            f"ρ={estimate:.3f}  {format_p_value(p_values[index])}",
            va="center",
            ha="left",
            color=INK,
            fontsize=8.5,
            transform=axis.get_yaxis_transform(),
        )

    axis.set_yticks(y_positions, labels)
    axis.invert_yaxis()
    axis.set_xlim(-0.06, 0.45)
    axis.set_xlabel("Spearman ρ: T-cell score vs. dysfunction score", color=INK)
    axis.set_title(
        "Per-sample effects and fixed-effect Fisher z pooled estimate",
        loc="left",
        fontsize=11,
        fontweight="bold",
        color=INK,
        pad=8,
    )
    axis.grid(axis="x", color=GRID, linewidth=0.8)
    axis.set_axisbelow(True)
    axis.tick_params(axis="both", colors=INK, labelsize=9)
    axis.spines[["top", "right", "left"]].set_visible(False)
    axis.spines["bottom"].set_color(MUTED)


def main() -> None:
    scores = pd.read_csv(RESULTS_DIR / "spatial_spot_level_scores.csv.gz")
    per_sample = pd.read_csv(RESULTS_DIR / "spatial_tcell_dysfunction_correlation_by_sample.csv")
    pooled = pd.read_csv(RESULTS_DIR / "spatial_tcell_dysfunction_correlation_pooled.csv").iloc[0]

    bundles = {
        gsm: load_tissue_bundle(gsm, label, scores)
        for gsm, label in SAMPLES
    }
    tcell_limits = scores["tcell_score"].quantile([0.02, 0.98]).to_numpy()
    dysfunction_limits = scores["dysfunction_score"].quantile([0.02, 0.98]).to_numpy()
    tcell_norm = Normalize(*tcell_limits, clip=True)
    dysfunction_norm = Normalize(*dysfunction_limits, clip=True)

    figure = plt.figure(figsize=(17, 11.5), facecolor="white")
    grid = GridSpec(
        3,
        5,
        figure=figure,
        height_ratios=[1, 1, 0.72],
        left=0.055,
        right=0.94,
        top=0.87,
        bottom=0.075,
        hspace=0.12,
        wspace=0.035,
    )

    figure.suptitle(
        "Spatial T-cell and dysfunction scores in SCLC tissue",
        x=0.055,
        y=0.965,
        ha="left",
        fontsize=20,
        fontweight="bold",
        color=INK,
    )
    figure.text(
        0.055,
        0.925,
        "GSE263196 · actual Visium H&E tissue images · 15,632 in-tissue spots · "
        "score colors clipped to pooled 2nd–98th percentiles",
        ha="left",
        fontsize=10.5,
        color=MUTED,
    )

    tcell_mappable = None
    dysfunction_mappable = None
    for column, (gsm, label) in enumerate(SAMPLES):
        image, spatial = bundles[gsm]
        statistics = per_sample.loc[per_sample["sample_gsm"].eq(gsm)].iloc[0]

        tcell_axis = figure.add_subplot(grid[0, column])
        draw_tissue_axis(tcell_axis, image, spatial, "tcell_score", tcell_norm, "Blues")
        tcell_axis.set_title(
            f"{label} · n={int(statistics['n_spots']):,}\n"
            f"ρ={statistics['rho']:.3f} [{statistics['ci_low']:.3f}, {statistics['ci_high']:.3f}]",
            fontsize=9.5,
            color=INK,
            pad=5,
        )

        dysfunction_axis = figure.add_subplot(grid[1, column])
        draw_tissue_axis(
            dysfunction_axis,
            image,
            spatial,
            "dysfunction_score",
            dysfunction_norm,
            "magma",
        )
        tcell_mappable = tcell_axis.collections[-1]
        dysfunction_mappable = dysfunction_axis.collections[-1]

    figure.text(
        0.018,
        0.72,
        "A  T-cell identity score",
        rotation=90,
        va="center",
        ha="center",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    figure.text(
        0.018,
        0.445,
        "B  Dysfunction score",
        rotation=90,
        va="center",
        ha="center",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )
    figure.text(
        0.018,
        0.185,
        "C  Validation effect",
        rotation=90,
        va="center",
        ha="center",
        fontsize=11,
        fontweight="bold",
        color=INK,
    )

    tcell_colorbar_axis = figure.add_axes([0.952, 0.625, 0.012, 0.19])
    tcell_colorbar = figure.colorbar(tcell_mappable, cax=tcell_colorbar_axis)
    tcell_colorbar.set_label("T-cell score", fontsize=9, color=INK)
    tcell_colorbar.ax.tick_params(labelsize=8, colors=INK)
    tcell_colorbar.outline.set_visible(False)

    dysfunction_colorbar_axis = figure.add_axes([0.952, 0.35, 0.012, 0.19])
    dysfunction_colorbar = figure.colorbar(dysfunction_mappable, cax=dysfunction_colorbar_axis)
    dysfunction_colorbar.set_label("Dysfunction score", fontsize=9, color=INK)
    dysfunction_colorbar.ax.tick_params(labelsize=8, colors=INK)
    dysfunction_colorbar.outline.set_visible(False)

    forest_axis = figure.add_subplot(grid[2, :])
    forest_axis.set_position([0.13, 0.075, 0.71, 0.20])
    draw_forest(forest_axis, per_sample, pooled)

    figure.text(
        0.94,
        0.025,
        "Scores: log-normalized Visium expression; T-cell markers and pre-registered "
        "dysfunction markers defined in spatial_validation.py. Tissue images: NCBI GEO GSE263196.",
        ha="right",
        fontsize=8.5,
        color=MUTED,
    )

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    png_path = FIGURE_DIR / "spatial_tissue_validation_panel.png"
    pdf_path = FIGURE_DIR / "spatial_tissue_validation_panel.pdf"
    figure.savefig(png_path, dpi=220, facecolor="white")
    figure.savefig(pdf_path, dpi=220, facecolor="white")
    plt.close(figure)
    print(f"Wrote {png_path}")
    print(f"Wrote {pdf_path}")


if __name__ == "__main__":
    main()
