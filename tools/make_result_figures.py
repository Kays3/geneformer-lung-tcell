#!/usr/bin/env python3
"""Regenerate the poster's data-derived figures from source CSVs/JSON.

Run whenever a source result table changes:

    python tools/make_result_figures.py

Each figure here is captioned or tabled elsewhere in the poster with the
SAME numbers, so keeping the figure regeneration script next to those
numbers (rather than hand-editing a saved PNG) is what keeps them from
drifting apart. See poster/README.md, "Confusion-matrix fix", for the bug
this discipline was written to prevent: an old figure from a different
cohort sat next to a table with the correct cohort's numbers for several
poster drafts before anyone noticed the mismatch.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parent.parent
SCLC = REPO / "sclc_validation"
OUT_PERTURB = SCLC / "perturbation_workflow" / "figures"
OUT_PERTURB.mkdir(parents=True, exist_ok=True)


def confusion_matrix_figure():
    """Fig. 1 — held-out confusion matrix, SCLC/Normal/LUAD.

    Source: sclc_validation/perturbation_workflow/results/test_confusion_matrix.csv
    Must stay in this class order (SCLC, Normal, LUAD) to match the table
    rendered next to it in poster_template.html.
    """
    cm = pd.read_csv(
        SCLC / "perturbation_workflow/results/test_confusion_matrix.csv", index_col=0
    )
    order = ["small cell lung carcinoma", "normal", "lung adenocarcinoma"]
    cm = cm.loc[order, order]
    labels = ["SCLC", "Normal", "LUAD"]
    counts = cm.values
    row_pct = counts / counts.sum(axis=1, keepdims=True) * 100

    fig, ax = plt.subplots(figsize=(6.2, 5.6), dpi=200)
    im = ax.imshow(row_pct, cmap="Blues", vmin=0, vmax=100)
    for i in range(3):
        for j in range(3):
            val = row_pct[i, j]
            ax.text(j, i, f"{counts[i, j]:,}\n{val:.1f}%", ha="center", va="center",
                    fontsize=13, fontweight="bold",
                    color="white" if val > 50 else "#1a1a1a")
    ax.set_xticks(range(3)); ax.set_xticklabels(labels, fontsize=12)
    ax.set_yticks(range(3)); ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlabel("Predicted disease", fontsize=12)
    ax.set_ylabel("Actual disease", fontsize=12)
    n_total = counts.sum()
    ax.set_title(f"Current T-cell classifier: held-out confusion matrix\n"
                 f"Labels show cell count and percentage within each actual class (n={n_total:,}).",
                 fontsize=12, fontweight="bold", loc="left")
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Row percentage", fontsize=11)
    fig.tight_layout()
    out = OUT_PERTURB / "sclc_confusion_matrix.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def bidirectional_scatter_figure():
    """Fig. 4 — delete vs. overexpress shift, 21-gene panel, all 6 comparisons.

    Source: perturbation_workflow/targeted_panel/results/
            targeted_panel_delete_overexpress_merged.csv
    A concordant hit (filled circle) falls in the upper-left or lower-right
    quadrant: deletion and overexpression move the cell in opposite
    directions along the same disease axis. Symlog axes are necessary - two
    positive-control genes (ASCL1, NEUROD1) have shifts an order of
    magnitude larger than the rest of the panel and collapse everything else
    to the origin on a linear scale.
    """
    merged = pd.read_csv(
        SCLC / "perturbation_workflow/targeted_panel/results/"
        "targeted_panel_delete_overexpress_merged.csv"
    )
    panel = merged[merged["gene_set"] == "panel"].copy()
    comparisons = sorted(panel["comparison"].unique())
    colors = plt.cm.tab10(np.linspace(0, 1, len(comparisons)))

    fig, ax = plt.subplots(figsize=(6.5, 6.8), dpi=200)
    for c, col in zip(comparisons, colors):
        sub = panel[panel["comparison"] == c]
        concordant = sub[sub["concordant"] == True]  # noqa: E712
        discordant = sub[sub["concordant"] == False]  # noqa: E712
        ax.scatter(discordant["delete_shift"], discordant["overexpress_shift"],
                   color=col, alpha=0.5, s=34, marker="x", linewidths=1.6)
        ax.scatter(concordant["delete_shift"], concordant["overexpress_shift"],
                   color=col, alpha=0.9, s=52, marker="o", label=c.replace("_", " "),
                   edgecolors="white", linewidths=0.4)

    lin_thresh = 0.02
    ax.set_xscale("symlog", linthresh=lin_thresh)
    ax.set_yscale("symlog", linthresh=lin_thresh)
    lim = 0.5
    ax.axhline(0, color="#999", lw=0.6, zorder=0)
    ax.axvline(0, color="#999", lw=0.6, zorder=0)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Delete: mean goal-state shift (symlog)", fontsize=12)
    ax.set_ylabel("Overexpress: mean goal-state shift (symlog)", fontsize=12)
    ax.set_title("Bidirectional perturbation: delete vs. overexpress shift\n"
                 "21-gene panel, all 6 comparisons  \u00b7  "
                 "filled circle = concordant, x = discordant",
                 fontsize=11, fontweight="bold", loc="left")
    ax.legend(fontsize=8, loc="lower left", framealpha=0.92, ncol=1)
    fig.tight_layout()
    out = OUT_PERTURB / "bidirectional_scatter.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def concordance_by_comparison_figure():
    """Fig. 5 — concordance rate by comparison and gene set.

    Source: targeted_panel_concordance_summary_by_comparison.csv
    n_concordant summed across all rows must equal the poster's headline
    "123 concordant hits" figure (51 from the 21-gene panel + 72 from the
    29-gene top-driver set) - this is asserted below rather than assumed.
    """
    conc = pd.read_csv(
        SCLC / "perturbation_workflow/targeted_panel/results/"
        "targeted_panel_concordance_summary_by_comparison.csv"
    )
    conc["pct"] = conc["n_concordant"] / conc["n_genes"] * 100
    total = int(conc["n_concordant"].sum())
    assert total == 123, (
        f"concordance summary totals {total}, expected 123 - the poster's "
        f"headline number and this figure have drifted apart; check "
        f"targeted_panel_concordance_summary_by_comparison.csv")

    comparisons_order = sorted(conc["comparison"].unique())
    gene_sets = conc["gene_set"].unique()
    set_labels = {"panel": "21-gene exhaustion panel",
                  "top_driver_luad_lusc_normal": "29-gene top-driver set"}
    set_colors = {"panel": "#2176ae", "top_driver_luad_lusc_normal": "#c77d2e"}

    fig, ax = plt.subplots(figsize=(6.8, 5.2), dpi=200)
    x = np.arange(len(comparisons_order))
    width = 0.35
    for i, gs in enumerate(gene_sets):
        sub = conc[conc["gene_set"] == gs].set_index("comparison").loc[comparisons_order]
        bars = ax.bar(x + (i - 0.5) * width, sub["pct"], width,
                       label=set_labels.get(gs, gs), color=set_colors.get(gs, f"C{i}"))
        for b, n_c, n_g in zip(bars, sub["n_concordant"], sub["n_genes"]):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 1.5,
                    f"{n_c}/{n_g}", ha="center", va="bottom", fontsize=8)

    ax.set_xticks(x)
    ax.set_xticklabels([c.replace("_", "\n") for c in comparisons_order], fontsize=9)
    ax.set_ylabel("Concordant genes (%)", fontsize=11)
    ax.set_ylim(0, 75)
    ax.set_title("Bidirectional concordance rate by comparison and gene set",
                 fontsize=11.5, fontweight="bold", loc="left")
    ax.legend(fontsize=9, loc="upper right", framealpha=0.9)
    ax.grid(axis="y", alpha=0.25)
    fig.tight_layout()
    out = OUT_PERTURB / "concordance_by_comparison.png"
    fig.savefig(out, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return out


def main():
    for fn in (confusion_matrix_figure, bidirectional_scatter_figure,
               concordance_by_comparison_figure):
        out = fn()
        print(f"wrote {out.relative_to(REPO)}  ({out.stat().st_size / 1e3:.0f} KB)")


if __name__ == "__main__":
    main()
