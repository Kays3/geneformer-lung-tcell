#!/usr/bin/env python3
"""Build figures, tables, and a reproducible notebook for perturbation exploration."""

from pathlib import Path
import json
import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
DEFAULT_ROOT = Path("/home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation")
ROOT = Path(os.environ.get("SCLC_PERTURBATION_ROOT", DEFAULT_ROOT))
STATS = ROOT / "stats"
FIGURES = HERE / "figures" / "exploration"
TABLES = HERE / "tables"
FIGURES.mkdir(parents=True, exist_ok=True)
TABLES.mkdir(parents=True, exist_ok=True)

COMPARISONS = ["sclc_to_luad", "sclc_to_normal", "luad_to_sclc", "luad_to_normal", "normal_to_sclc", "normal_to_luad"]


def load_balanced() -> pd.DataFrame:
    frames = []
    for comparison in COMPARISONS:
        delete = pd.read_csv(STATS / "delete" / f"heldout_allgene_delete_{comparison}.csv").drop(columns=["Unnamed: 0"], errors="ignore")
        over = pd.read_csv(STATS / "overexpress" / f"heldout_allgene_overexpress_{comparison}.csv").drop(columns=["Unnamed: 0"], errors="ignore")
        delete = delete[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(columns={"Shift_to_goal_end": "delete_shift", "Goal_end_FDR": "delete_fdr", "N_Detections": "delete_n"})
        over = over[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(columns={"Shift_to_goal_end": "overexpress_shift", "Goal_end_FDR": "overexpress_fdr", "N_Detections": "overexpress_n"})
        frame = delete.merge(over, on=["Gene_name", "Ensembl_ID"], how="inner")
        frame["comparison"] = comparison
        frame = frame[(frame.delete_fdr < 0.05) & (frame.overexpress_fdr < 0.05) & (frame.delete_shift * frame.overexpress_shift < 0) & (frame.delete_n >= 25) & (frame.overexpress_n >= 25)].copy()
        frame["min_detection"] = frame[["delete_n", "overexpress_n"]].min(axis=1)
        frame["min_abs_shift"] = frame[["delete_shift", "overexpress_shift"]].abs().min(axis=1)
        frame["mean_abs_shift"] = frame[["delete_shift", "overexpress_shift"]].abs().mean(axis=1)
        frame["detection_norm"] = frame.min_detection / frame.min_detection.max()
        shift_cap = frame.min_abs_shift.quantile(0.95)
        frame["shift_norm"] = (frame.min_abs_shift / shift_cap).clip(upper=1)
        frame["balanced_score"] = np.sqrt(frame.detection_norm * frame.shift_norm)
        frames.append(frame)
    return pd.concat(frames, ignore_index=True).sort_values(["comparison", "balanced_score", "min_detection"], ascending=[True, False, False])


def functional_label(gene: str) -> str:
    if gene in {"B2M", "HLA-A", "HLA-B", "HLA-C", "CD74", "HLA-DRA", "HLA-DRB1"} or gene.startswith("HLA-"):
        return "Antigen presentation / immune";
    if gene in {"JUN", "JUNB", "FOS", "KLF2", "KLF6", "MAP2K2", "TSC22D3", "CXCR4", "RORA"}:
        return "Signaling / transcription";
    if gene in {"S100A4", "S100A6", "S100A10", "VIM", "CD81", "TMSB10", "TMSB4X"}:
        return "Cytoskeleton / motility";
    if gene in {"HSPA1A", "HSPA1B", "HSPA8", "HSPB1", "DNAJB1", "MIF"}:
        return "Stress / proteostasis";
    if gene.startswith("MT-"):
        return "Mitochondrial";
    if gene.startswith("RPL") or gene.startswith("RPS"):
        return "Translation / ribosome";
    if gene in {"SCGB1A1", "WFDC2", "MMP7", "TFF3", "KRT17", "SFTPB", "SFTPC", "NAPSA"}:
        return "Epithelial / ambient-RNA caution";
    if gene in {"NBEAL1", "S100A4", "S100A10", "PPP1R14B", "H4C3"}:
        return "Other prioritized signal";
    return "Other / unannotated here"


def build_artifacts() -> tuple[pd.DataFrame, pd.DataFrame]:
    balanced = load_balanced()
    balanced["functional_category"] = balanced.Gene_name.map(functional_label)
    balanced.to_csv(TABLES / "exploration_balanced_rankings.csv", index=False)
    top = balanced.sort_values(["balanced_score", "min_detection"], ascending=False).head(100)
    functional = top.groupby("functional_category", as_index=False).agg(
        genes=("Gene_name", "nunique"),
        comparison_rows=("Gene_name", "size"),
        median_score=("balanced_score", "median"),
        median_min_detection=("min_detection", "median"),
    ).sort_values(["genes", "median_score"], ascending=False)
    functional.to_csv(TABLES / "exploration_functional_summary.csv", index=False)
    balanced.groupby("comparison", group_keys=False).head(15).to_csv(TABLES / "exploration_top_candidates_by_comparison.csv", index=False)

    fig, ax = plt.subplots(figsize=(10, 6))
    for comparison, group in balanced.groupby("comparison"):
        ax.scatter(group.min_detection, group.min_abs_shift, s=10 + 70 * group.balanced_score, alpha=0.35, label=comparison)
    labels = balanced.sort_values("balanced_score", ascending=False).head(10)
    ax.scatter(labels.min_detection, labels.min_abs_shift, facecolors="none", edgecolors="#111827", linewidths=1.2, s=100)
    for _, row in labels.iterrows():
        ax.annotate(row.Gene_name, (row.min_detection, row.min_abs_shift), xytext=(4, 4), textcoords="offset points", fontsize=8)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Minimum detections across delete and overexpression (log scale)")
    ax.set_ylabel("Minimum absolute cosine shift (log scale)")
    ax.set_title("Balanced bidirectional perturbation candidates")
    ax.grid(alpha=0.2)
    ax.legend(fontsize=7, ncol=2, frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURES / "balanced_detection_shift_scatter.png", dpi=180)
    plt.close(fig)

    top_per_comparison = balanced.groupby("comparison", group_keys=False).head(8)
    pivot = top_per_comparison.pivot_table(index="Gene_name", columns="comparison", values="balanced_score", aggfunc="max", fill_value=0)
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).head(25).index]
    fig, ax = plt.subplots(figsize=(12, 9))
    image = ax.imshow(pivot.values, cmap="viridis", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), pivot.columns, rotation=40, ha="right")
    ax.set_yticks(range(len(pivot.index)), pivot.index)
    ax.set_title("Balanced score across disease comparisons")
    fig.colorbar(image, ax=ax, label="Balanced score")
    fig.tight_layout()
    fig.savefig(FIGURES / "balanced_score_heatmap.png", dpi=180)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5))
    functional.sort_values("genes").plot.barh(x="functional_category", y="genes", ax=ax, color="#3b6ea8", legend=False)
    ax.set_xlabel("Rows among overall top 100 balanced candidates")
    ax.set_ylabel("")
    ax.set_title("Functional categories among prioritized candidates")
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURES / "functional_categories.png", dpi=180)
    plt.close(fig)

    # Opposite quadrants show bidirectional shifts; dot size is minimum detection.
    fig, axes = plt.subplots(2, 3, figsize=(18, 11), squeeze=False)
    for ax, comparison in zip(axes.flat, COMPARISONS):
        group = balanced[balanced.comparison == comparison]
        sizes = 12 + 260 * (group.min_detection / group.min_detection.max()) ** 0.7
        ax.scatter(group.delete_shift, group.overexpress_shift, s=sizes, c=group.balanced_score, cmap="viridis", alpha=0.42, linewidths=0)
        top = group.head(6)
        ax.scatter(top.delete_shift, top.overexpress_shift, s=50 + 260 * (top.min_detection / group.min_detection.max()) ** 0.7, facecolors="none", edgecolors="#111827", linewidths=1.1)
        for _, row in top.iterrows():
            ax.annotate(row.Gene_name, (row.delete_shift, row.overexpress_shift), xytext=(3, 3), textcoords="offset points", fontsize=7)
        lim = max(abs(group.delete_shift).max(), abs(group.overexpress_shift).max()) * 1.08
        ax.axhline(0, color="#9ca3af", linewidth=0.7); ax.axvline(0, color="#9ca3af", linewidth=0.7)
        ax.plot([-lim, lim], [lim, -lim], color="#d1d5db", linewidth=0.7, linestyle="--")
        ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
        ax.set_title(comparison.replace("_to_", " → "))
        ax.set_xlabel("Delete cosine shift"); ax.set_ylabel("Overexpression cosine shift")
        ax.grid(alpha=0.15)
    fig.suptitle("Bidirectional candidates by comparison\n(dot size = minimum detection; color = balanced score)", y=1.01, fontsize=14)
    fig.tight_layout()
    fig.savefig(FIGURES / "bidirectional_candidates_by_comparison.png", dpi=180, bbox_inches="tight")
    plt.close(fig)
    return balanced, functional


def md(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(True)}


def code(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source.splitlines(True)}


def build_notebook() -> Path:
    notebook = {
        "cells": [
            md("# Exploratory bidirectional perturbation analysis\n\nThis notebook explores the completed donor-held-out Geneformer delete and overexpression results. It ranks genes by a balanced combination of cosine-shift strength and detection coverage, then provides functional labels for the prioritized candidates. These are model-level sensitivity signals, not experimental causal effects."),
            md("## TL;DR\n\nThe balanced score uses the weaker perturbation arm for both effect size and coverage. This prevents highly detected housekeeping genes with tiny shifts from dominating, while avoiding unstable large shifts supported by only a few cells."),
            code("from pathlib import Path\nimport os\nimport numpy as np\nimport pandas as pd\nfrom IPython.display import Image, display\n\nHERE = Path.cwd()\nROOT = Path(os.environ.get('SCLC_PERTURBATION_ROOT', '/home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation'))\nSTATS = ROOT / 'stats'\nTABLES = HERE / 'tables'\nFIGURES = HERE / 'figures' / 'exploration'\nrankings = pd.read_csv(TABLES / 'exploration_balanced_rankings.csv')\nfunctional_summary = pd.read_csv(TABLES / 'exploration_functional_summary.csv')\ntop_by_comparison = pd.read_csv(TABLES / 'exploration_top_candidates_by_comparison.csv')\nrankings.head(10)"),
            md("## Methods and quality checks\n\nA candidate must be significant in both arms (`FDR < 0.05`), have opposite-signed shifts, and have at least 25 detections in each arm. `min_abs_shift` is the smaller absolute shift across the two arms; `min_detection` is the smaller detection count. Both are normalized within each comparison, with shift values capped at the 95th percentile. The final score is `sqrt(normalized_detection * normalized_shift)`."),
            code("assert rankings['balanced_score'].between(0, 1).all()\nassert (rankings['min_detection'] >= 25).all()\nassert (rankings['delete_shift'] * rankings['overexpress_shift'] < 0).all()\nprint('Qualified comparison rows:', len(rankings))\nprint('Unique genes:', rankings['Ensembl_ID'].nunique())\nrankings.groupby('comparison').size().rename('qualified_rows').to_frame()"),
            md("## Balanced ranking\n\nUse `comparison` to focus on one disease transition. The table below shows the highest-ranked candidates overall; the saved CSV contains all qualified rows."),
            code("rankings.sort_values(['balanced_score', 'min_detection'], ascending=False)[['comparison', 'Gene_name', 'delete_shift', 'overexpress_shift', 'min_abs_shift', 'min_detection', 'balanced_score']].head(30)"),
            md("## Figures\n\nThe scatter plot shows the two objectives directly. Point size encodes the balanced score; labels identify the overall top candidates."),
            code("display(Image(filename=str(FIGURES / 'balanced_detection_shift_scatter.png')))\ndisplay(Image(filename=str(FIGURES / 'balanced_score_heatmap.png')))\ndisplay(Image(filename=str(FIGURES / 'bidirectional_candidates_by_comparison.png')))\ndisplay(Image(filename=str(FIGURES / 'functional_categories.png')))"),
            md("## Candidates by comparison\n\nEach row below is one of the top 15 balanced candidates for its disease transition. In the six-panel figure, marker size is the minimum detection count across deletion and overexpression."),
            code("top_by_comparison[['comparison', 'Gene_name', 'delete_shift', 'overexpress_shift', 'min_detection', 'balanced_score', 'functional_category']]") ,
            md("## Gene functionality\n\nFunctional labels are transparent, coarse categories for exploration. They are based on curated symbol sets and naming patterns in this notebook; they are not a replacement for GO/Reactome enrichment or literature review. The epithelial/ambient-RNA category is intentionally flagged as a contamination-sensitive interpretation."),
            code("functional_summary"),
            code("FUNCTIONAL_LABELS = {\n    'B2M': 'Antigen presentation / immune', 'HLA-C': 'Antigen presentation / immune',\n    'CD74': 'Antigen presentation / immune', 'JUN': 'Signaling / transcription',\n    'KLF6': 'Signaling / transcription', 'MAP2K2': 'Signaling / transcription',\n    'CXCR4': 'Signaling / transcription', 'S100A4': 'Cytoskeleton / motility',\n    'VIM': 'Cytoskeleton / motility', 'CD81': 'Cytoskeleton / motility',\n    'HSPA1A': 'Stress / proteostasis', 'HSPA1B': 'Stress / proteostasis',\n    'DNAJB1': 'Stress / proteostasis', 'NBEAL1': 'Other prioritized signal',\n}\ndef gene_lookup(gene):\n    cols = ['comparison', 'Gene_name', 'delete_shift', 'overexpress_shift', 'min_detection', 'balanced_score', 'functional_category']\n    return rankings[rankings.Gene_name.eq(gene)].sort_values('balanced_score', ascending=False)[cols]\n\ngene_lookup('HLA-C')"),
            md("## Takeaways and caveats\n\n- Prefer candidates that recur across comparisons with high `min_detection` and a substantial `min_abs_shift`.\n- `HLA-C`, `B2M`, `NBEAL1`, `S100A4`, `CXCR4`, `JUN`, and `MAP2K2` are useful starting points for review.\n- Ribosomal, mitochondrial, heat-shock, and housekeeping genes can be statistically strong but may reflect global cell state or technical effects.\n- Epithelial markers require ambient-RNA, doublet, and epithelial-burden sensitivity checks.\n- Donor-stratified stability and independent biological validation are required before causal interpretation."),
        ],
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"}},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    output = HERE / "perturbation_results_exploration.ipynb"
    output.write_text(json.dumps(notebook, indent=1) + "\n")
    return output


if __name__ == "__main__":
    balanced, functional = build_artifacts()
    output = build_notebook()
    print(f"Wrote {output}")
    print(f"Qualified rows: {len(balanced)}; unique genes: {balanced.Ensembl_ID.nunique()}")
    print(f"Functional categories: {len(functional)}")
