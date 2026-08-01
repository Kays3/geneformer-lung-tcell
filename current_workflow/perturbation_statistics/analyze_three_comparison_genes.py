#!/usr/bin/env python3
"""Visualize genes with the maximum observed cross-comparison coverage."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
INPUT = HERE / "source_tables" / "top_goal_shift_genes.csv"
OUTPUT = HERE / "three_comparison_genes"

COLORS = {
    "LUAD → LUSC": "#2f6690",
    "NORMAL → LUAD": "#d39b2a",
    "NORMAL → LUSC": "#8a5a9e",
}

KEGG_PATHWAYS = [
    ("03040", "Spliceosome", "Genetic information"),
    ("04010", "MAPK signaling", "Signal transduction"),
    ("04141", "Protein processing in ER", "Proteostasis"),
    ("04144", "Endocytosis", "Transport / catabolism"),
    ("04213", "Longevity regulation", "Stress response"),
    ("04612", "Antigen processing & presentation", "Immune system"),
    ("04915", "Estrogen signaling", "Endocrine signaling"),
    ("05020", "Prion disease", "Disease pathway"),
    ("05134", "Legionellosis", "Infectious disease"),
    ("05145", "Toxoplasmosis", "Infectious disease"),
    ("05162", "Measles", "Infectious disease"),
    ("05417", "Lipid and atherosclerosis", "Cardiovascular disease"),
]


def focused_data() -> pd.DataFrame:
    frame = pd.read_csv(INPUT)
    coverage = frame.groupby("Gene_name")["comparison_label"].nunique()
    maximum = int(coverage.max())
    genes = coverage.index[coverage == maximum]
    focused = frame.loc[frame["Gene_name"].isin(genes)].copy()
    focused["minus_log10_fdr"] = -np.log10(focused["Goal_end_FDR"].clip(lower=np.finfo(float).tiny))
    focused["max_comparison_coverage"] = maximum
    return focused.sort_values(["Gene_name", "N_Detections"])


def plot_detection_effect(frame: pd.DataFrame, path: Path) -> None:
    genes = sorted(frame["Gene_name"].unique())
    fig, axes = plt.subplots(1, len(genes), figsize=(11.4, 5.4), sharey=False)
    if len(genes) == 1:
        axes = [axes]
    size_min, size_max = 130, 520
    fdr = frame["minus_log10_fdr"]
    scaled = size_min + (size_max - size_min) * (fdr - fdr.min()) / (fdr.max() - fdr.min())
    frame = frame.assign(marker_size=scaled)

    for ax, gene in zip(axes, genes):
        subset = frame.loc[frame["Gene_name"] == gene]
        for row in subset.itertuples(index=False):
            ax.scatter(row.N_Detections, row.Shift_to_goal_end, s=row.marker_size,
                       color=COLORS[row.comparison_label], edgecolor="#263238",
                       linewidth=0.9, alpha=0.9)
            ax.annotate(row.comparison_label, (row.N_Detections, row.Shift_to_goal_end),
                        xytext=(6, 5), textcoords="offset points", fontsize=8.5)
        ax.set_title(gene, weight="bold")
        ax.set_xlabel("Detected cells")
        ax.grid(color="#d9dde2", linewidth=0.7, alpha=0.7)
        ax.set_axisbelow(True)
        ax.set_xlim(0, subset["N_Detections"].max() * 1.22)
        ax.set_ylim(0, subset["Shift_to_goal_end"].max() * 1.18)
        ax.set_ylabel("Shift toward goal endpoint")
    fig.suptitle("Detection support versus goal-shift effect in three comparisons",
                 x=0.06, ha="left", weight="bold", fontsize=15)
    fig.text(0.06, 0.01,
             "Bubble area encodes −log10(FDR); all six points pass the published qualified-shift filters.",
             fontsize=9, color="#4f5963")
    fig.tight_layout(rect=(0, 0.05, 1, 0.93))
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_kegg_membership(path: Path) -> None:
    pathways = pd.DataFrame(KEGG_PATHWAYS, columns=["kegg_id", "pathway", "category"])
    category_order = list(dict.fromkeys(pathways["category"]))
    palette = dict(zip(category_order, plt.cm.cividis(np.linspace(0.12, 0.88, len(category_order)))))

    fig, (ax, note_ax) = plt.subplots(1, 2, figsize=(12.4, 7.2),
                                      gridspec_kw={"width_ratios": [3.3, 1.35]})
    y = np.arange(len(pathways))[::-1]
    ax.hlines(y, 0, 1, color="#d9dde2", linewidth=1.0)
    ax.scatter(np.ones(len(pathways)), y, s=165,
               c=[palette[c] for c in pathways["category"]],
               edgecolor="#263238", linewidth=0.8)
    for yi, row in zip(y, pathways.itertuples(index=False)):
        ax.text(0.97, yi, f"hsa{row.kegg_id}  {row.pathway}", ha="right", va="center", fontsize=9.5)
    ax.text(1, len(pathways) + 0.15, "HSPA1B", ha="center", weight="bold", fontsize=11)
    ax.set_xlim(-0.02, 1.12)
    ax.set_ylim(-0.8, len(pathways) + 0.8)
    ax.axis("off")

    note_ax.axis("off")
    note_ax.text(0.02, 0.92, "SFTPC", fontsize=15, weight="bold", transform=note_ax.transAxes)
    note_ax.text(0.02, 0.82, "KEGG KO: K26068", fontsize=10.5, transform=note_ax.transAxes)
    note_ax.text(0.02, 0.72, "No canonical KEGG pathway\nassignment", fontsize=12,
                 color="#8a5a2b", weight="bold", transform=note_ax.transAxes)
    note_ax.text(0.02, 0.57,
                 "KEGG classifies SFTPC as a structural\nprotein and links it to pulmonary\ndisease entries rather than a pathway.",
                 fontsize=10, linespacing=1.5, transform=note_ax.transAxes)
    note_ax.add_patch(plt.Rectangle((0, 0.49), 0.96, 0.49, fill=False,
                                    edgecolor="#aab2ba", linewidth=1.1,
                                    transform=note_ax.transAxes))
    fig.suptitle("KEGG pathway membership of maximum-coverage genes",
                 x=0.04, ha="left", weight="bold", fontsize=16)
    fig.text(0.04, 0.02,
             "Memberships are from KEGG human gene entries hsa:3304 (HSPA1B) and hsa:6440 (SFTPC).",
             fontsize=9, color="#4f5963")
    fig.tight_layout(rect=(0, 0.05, 1, 0.94))
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_normal_goal_comparison(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(INPUT)
    normal = frame.loc[frame["comparison"].isin(["normal_to_luad", "normal_to_lusc"])].copy()
    shifts = normal.pivot(index=["Gene_name", "Ensembl_ID"], columns="comparison",
                          values="Shift_to_goal_end").dropna()
    cells = normal.pivot(index=["Gene_name", "Ensembl_ID"], columns="comparison",
                         values="N_Detections").dropna()
    wide = shifts.join(cells, lsuffix="_shift", rsuffix="_cells").reset_index()
    wide = wide.rename(columns={
        "normal_to_luad_shift": "normal_to_luad_shift",
        "normal_to_lusc_shift": "normal_to_lusc_shift",
        "normal_to_luad_cells": "normal_to_luad_cells",
        "normal_to_lusc_cells": "normal_to_lusc_cells",
    })
    wide["detected_cells"] = wide[["normal_to_luad_cells", "normal_to_lusc_cells"]].min(axis=1)
    wide["maximum_coverage_gene"] = wide["Gene_name"].isin(["HSPA1B", "SFTPC"])

    fig, ax = plt.subplots(figsize=(9.2, 7.2))
    sizes = 110 + 920 * np.sqrt(wide["detected_cells"] / wide["detected_cells"].max())
    colors = np.where(wide["maximum_coverage_gene"], "#d39b2a", "#4f7da3")
    ax.scatter(wide["normal_to_lusc_shift"], wide["normal_to_luad_shift"], s=sizes,
               c=colors, edgecolor="#263238", linewidth=1.0, alpha=0.9)
    offsets = {
        "FBLN1": (8, -18),
        "HSPA1B": (12, -37),
        "MGP": (8, 12),
        "SFTPA1": (14, 14),
        "SFTPA2": (-66, -22),
        "SFTPC": (10, 10),
    }
    for row in wide.itertuples(index=False):
        ax.annotate(f"{row.Gene_name}\n{int(row.detected_cells)} cells",
                    (row.normal_to_lusc_shift, row.normal_to_luad_shift),
                    xytext=offsets.get(row.Gene_name, (7, 5)), textcoords="offset points",
                    fontsize=9, weight="bold" if row.maximum_coverage_gene else "normal",
                    arrowprops={"arrowstyle": "-", "color": "#68737d", "linewidth": 0.7})
    diagonal_max = max(wide["normal_to_lusc_shift"].max(), wide["normal_to_luad_shift"].max()) * 1.12
    ax.plot([0, diagonal_max], [0, diagonal_max], color="#68737d", linewidth=1.0,
            linestyle="--", label="Equal LUAD and LUSC goal shift")
    ax.set_xlim(0, wide["normal_to_lusc_shift"].max() * 1.42)
    ax.set_ylim(0, wide["normal_to_luad_shift"].max() * 1.18)
    ax.set_xlabel("NORMAL → LUSC shift toward goal endpoint")
    ax.set_ylabel("NORMAL → LUAD shift toward goal endpoint")
    ax.set_title("Shared NORMAL-source perturbations: LUAD versus LUSC goal shifts",
                 loc="left", weight="bold")
    ax.grid(color="#d9dde2", linewidth=0.7, alpha=0.7)
    ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left")
    fig.text(0.01, 0.01,
             "Bubble area and labels show detected cells; gold marks genes also recurring in a third comparison.",
             fontsize=9, color="#4f5963")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    return wide


def write_report(frame: pd.DataFrame, path: Path) -> None:
    rows = []
    for row in frame.sort_values(["Gene_name", "comparison_label"]).itertuples(index=False):
        rows.append(
            f"| {row.Gene_name} | {row.comparison_label} | {row.N_Detections} | "
            f"{row.Shift_to_goal_end:.4f} | {row.Goal_end_FDR:.2e} |"
        )
    text = f"""# Maximum-coverage perturbations: SFTPC and HSPA1B

## Technical summary

`SFTPC` and `HSPA1B` are the only genes appearing among the published top 15 qualified
toward-goal shifts in three comparisons—the maximum coverage observed. `HSPA1B` has the
largest detection counts (270–641 cells), while `SFTPC` has the larger goal-shift effects,
especially NORMAL → LUAD (0.1644).

## Detection support and effect size

![Detection versus effect plot](detection_effect_2d.png)

| Gene | Comparison | Detected cells | Goal shift | FDR |
|---|---|---:|---:|---:|
{chr(10).join(rows)}

## NORMAL-source shifts are consistently larger toward LUAD than toward LUSC

![NORMAL to LUAD versus NORMAL to LUSC](normal_luad_vs_lusc.png)

Six genes occur in both published NORMAL-source top-15 lists. Every point lies above the
equal-shift diagonal, meaning its modeled deletion shift toward LUAD is larger than its shift
toward LUSC. `SFTPC` has the largest shift on both axes (465 detected cells), while `HSPA1B`
has the largest detected-cell support (641 cells). Bubble area and labels encode detected
cells. Gold identifies `SFTPC` and `HSPA1B`, which also recur in LUAD → LUSC.

## Biological interpretation

### HSPA1B points to stress/proteostasis, but not a disease-specific mechanism by itself

HSPA1B encodes an inducible HSP70-family chaperone that stabilizes proteins against
aggregation and assists folding of newly translated proteins. KEGG maps it to protein
processing in the endoplasmic reticulum, MAPK signaling, endocytosis, antigen processing
and presentation, and several stress/disease pathways. Its high cell coverage therefore
supports a broadly shared stress/proteostasis axis. Because this is a ubiquitous stress
response gene, deletion effects should not automatically be interpreted as NSCLC-specific
T-cell biology.

### SFTPC is a lung epithelial signal and a contamination-sensitive result in T cells

SFTPC encodes pulmonary surfactant protein C, is strongly lung-restricted, and is essential
for alveolar surfactant function and lung homeostasis. In a nominal T-cell cohort, its repeated
appearance is therefore more consistent with an alveolar epithelial transcript burden,
ambient RNA, or epithelial–immune doublets than with a canonical T-cell-intrinsic pathway.
This is an inference from tissue specificity and should be tested directly with per-cell
SFTPC counts, epithelial marker burden, and doublet/decontamination sensitivity analyses.

## KEGG pathway visualization

![KEGG pathway membership](kegg_pathway_membership.png)

KEGG assigns HSPA1B to 12 human pathways. KEGG assigns SFTPC to orthology K26068 and pulmonary
disease entries but explicitly places it outside canonical pathway/BRITE categories; the figure
retains that absence rather than inventing a pathway connection.

## Scope and limitations

The analysis uses the 90-row `top_goal_shift_genes.csv`, which contains only the top 15
qualified shifts per comparison. It does not establish that these genes are absent from the
other three comparisons. Perturbation shifts are model-derived associations, not evidence of
causality or therapeutic tractability.

## Recommended next steps

1. Re-run this analysis on all six full perturbation tables.
2. For SFTPC, stratify cells by epithelial-marker burden and repeat after ambient-RNA correction
   and doublet removal.
3. For HSPA1B, test donor consistency and correlate the signal with broader heat-shock and
   unfolded-protein-response modules.

## Sources

- KEGG HSPA1B: https://www.kegg.jp/entry/hsa:3304
- KEGG SFTPC: https://www.kegg.jp/entry/hsa:6440
- KEGG protein processing in ER: https://www.kegg.jp/pathway/hsa04141
- KEGG antigen processing and presentation: https://www.kegg.jp/pathway/hsa04612
- NCBI HSPA1B: https://www.ncbi.nlm.nih.gov/gene/3304
- NCBI SFTPC: https://www.ncbi.nlm.nih.gov/gene/6440
"""
    path.write_text(text, encoding="utf-8")


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    frame = focused_data()
    frame.to_csv(OUTPUT / "three_comparison_genes.csv", index=False)
    plot_detection_effect(frame, OUTPUT / "detection_effect_2d.png")
    normal_wide = plot_normal_goal_comparison(OUTPUT / "normal_luad_vs_lusc.png")
    normal_wide.to_csv(OUTPUT / "normal_luad_vs_lusc.csv", index=False)
    plot_kegg_membership(OUTPUT / "kegg_pathway_membership.png")
    write_report(frame, OUTPUT / "REPORT.md")
    print(f"Wrote {OUTPUT / 'REPORT.md'}")


if __name__ == "__main__":
    main()
