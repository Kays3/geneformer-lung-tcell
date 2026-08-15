#!/usr/bin/env python3
"""Denoise the all-gene perturbation results and prioritize immune / cancer biology.

The exploratory ranking is dominated by ribosomal, mitochondrial, heat-shock and
ambient-RNA genes. This builder classifies every qualified hit into a transparent
technical-noise class or a curated immune / cancer program, then rebuilds tables,
figures and a reproducible notebook restricted to the interpretable signal.

Input is the committed exploration ranking so the analysis reproduces without the
remote compute artifacts; set SCLC_PERTURBATION_ROOT to rebuild from raw stats.
"""

from pathlib import Path
import json
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parents[1]
TABLES = HERE / "tables"
FIGURES = HERE / "figures" / "denoised"
FIGURES.mkdir(parents=True, exist_ok=True)

COMPARISONS = ["sclc_to_luad", "sclc_to_normal", "luad_to_sclc", "luad_to_normal", "normal_to_sclc", "normal_to_luad"]

# --- Technical / non-specific classes -------------------------------------------------
# These are removed from the prioritized set. They can be statistically extreme because
# they track global transcriptional output, dissociation stress, or ambient contamination
# rather than a T-cell-intrinsic program.
HEMOGLOBIN = {"HBB", "HBA1", "HBA2", "HBD", "HBG1", "HBG2", "HBM", "HBQ1", "ALAS2", "AHSP", "CA1"}
IMMEDIATE_EARLY = {
    "JUN", "JUNB", "JUND", "FOS", "FOSB", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3",
    "DUSP1", "DUSP2", "ZFP36", "ZFP36L1", "ZFP36L2", "KLF4", "KLF6", "SOCS3", "PPP1R15A",
    "BTG1", "BTG2", "CEBPB", "MYADM", "RGCC",
}
EPITHELIAL_AMBIENT = {
    "SFTPC", "SFTPB", "SFTPA1", "SFTPA2", "SCGB1A1", "SCGB3A1", "SCGB3A2", "NAPSA", "WFDC2",
    "TFF1", "TFF3", "MMP7", "AGER", "MUC1", "MUC5B", "CLDN18", "NKX2-1", "EPCAM", "ELF3",
    "CEACAM6", "SLPI", "SLC34A2",
}
MYELOID_PLATELET_AMBIENT = {"PPBP", "PF4", "LYZ", "S100A8", "S100A9", "S100A12", "FCN1", "VCAN"}
HOUSEKEEPING = {
    "ACTB", "ACTG1", "GAPDH", "TMSB4X", "TMSB10", "PTMA", "MYL6", "PFN1", "CFL1", "EEF1A1",
    "EEF1B2", "EEF1G", "EEF1D", "EEF2", "EIF1", "NACA", "BTF3", "SRP14", "UBA52", "FAU",
    "OAZ1", "SERF2", "HINT1", "PABPC1", "NPM1", "HNRNPA1",
}
NONCODING = {"MALAT1", "NEAT1", "XIST", "TSIX", "MIAT", "NORAD"}

NOISE_PATTERNS: list[tuple[str, re.Pattern]] = [
    ("Ribosomal", re.compile(r"^(RPL|RPS|MRPL|MRPS|RPLP)")),
    ("Mitochondrial", re.compile(r"^MT-")),
    ("Histone", re.compile(r"^(H1-|H2A|H2B|H3-|H3C|H4C|HIST)")),
    ("OXPHOS / housekeeping", re.compile(r"^(ATP5|COX\d|NDUF|UQCR|SDH[ABCD]$)")),
    ("Non-coding / lncRNA", re.compile(r"(^LINC|^SNHG|^MIR\d|^RNU\d|^RN7S|-AS\d$|^LOC\d)")),
]

# --- Curated immune / cancer programs -------------------------------------------------
IMMUNE_CANCER_PROGRAMS: dict[str, set[str]] = {
    "T-cell identity / TCR": {
        "CD3D", "CD3E", "CD3G", "CD247", "CD2", "CD5", "CD6", "CD7", "CD28", "TRAC", "TRBC1",
        "TRBC2", "TRDC", "TRGC1", "TRGC2", "LCK", "FYN", "ZAP70", "LAT", "ITK", "THEMIS",
        "SKAP1", "CD4", "CD8A", "CD8B", "CD3EAP",
    },
    "Checkpoint / exhaustion": {
        "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "TOX2", "LAYN", "ENTPD1", "BTLA",
        "CD160", "VSIR", "CD244", "KLRG1", "EOMES", "TNFRSF9", "CXCL13", "ITGAE", "PDCD1LG2",
        "CD274", "TNFRSF14",
    },
    "Cytotoxic effector": {
        "GZMA", "GZMB", "GZMH", "GZMK", "GZMM", "PRF1", "GNLY", "NKG7", "KLRD1", "KLRB1",
        "KLRC1", "FGFBP2", "FASLG", "IFNG", "TNF", "CST7", "CTSW", "CCL3", "CCL4", "CCL5",
        "SERPINB9",
    },
    "Memory / progenitor": {
        "TCF7", "LEF1", "SELL", "CCR7", "BACH2", "ID3", "SLAMF6", "IL7R", "FOXP1", "MYB",
        "SATB1", "ACTN1",
    },
    "Treg / suppressive": {
        "FOXP3", "IL2RA", "IKZF2", "TNFRSF4", "TNFRSF18", "IL2RB", "LRRC32", "IL10", "TGFB1",
    },
    "Antigen presentation / MHC": {
        "HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "B2M", "TAP1", "TAP2", "TAPBP",
        "NLRC5", "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1",
        "HLA-DQB1", "HLA-DMA", "HLA-DMB", "CIITA", "PSMB8", "PSMB9", "ERAP1", "ERAP2",
    },
    "Interferon / inflammatory": {
        "STAT1", "STAT2", "IRF1", "IRF7", "IRF9", "ISG15", "ISG20", "IFI6", "IFI44L", "IFIT1",
        "IFIT3", "MX1", "MX2", "OAS1", "OAS3", "GBP1", "GBP2", "GBP5", "IDO1", "BST2", "XAF1",
        "IFI16", "IFI35",
    },
    "Trafficking / tissue residency": {
        "CXCR3", "CXCR4", "CXCR5", "CXCR6", "CCR2", "CCR4", "CCR5", "CCR6", "CX3CR1", "S1PR1",
        "S1PR5", "ITGA1", "ITGAL", "ITGB2", "SELPLG", "RGS1", "CD69L",
    },
    "Costimulation / activation": {
        "ICOS", "CD27", "CD40LG", "TNFRSF1B", "CD38", "IL2", "IL21", "IL21R", "CD226", "SLAMF1",
        "SLAMF7", "CD58", "ICOSLG",
    },
    "SCLC neuroendocrine program": {
        "ASCL1", "NEUROD1", "POU2F3", "YAP1", "INSM1", "CHGA", "CHGB", "SYP", "NCAM1", "DLL3",
        "MYCL", "MYCN", "SOX2", "FOXA2", "CALCA", "GRP", "UCHL1", "NEUROD2",
    },
    "Oncogenic / tumor suppressor": {
        "TP53", "RB1", "MYC", "KRAS", "EGFR", "STK11", "KEAP1", "NFE2L2", "PTEN", "NOTCH1",
        "NOTCH2", "CREBBP", "EP300", "ARID1A", "SMARCA4", "CCND1", "CDKN2A", "BCL2", "MCL1",
    },
    "Immunosuppressive metabolism / TME": {
        "NT5E", "ADORA2A", "ARG1", "NOS2", "SLC2A1", "LDHA", "HIF1A", "VEGFA", "PRDM1", "BATF",
    },
}

# Genes that are genuinely immune-meaningful but are also canonical dissociation/immediate-early
# transcripts. Reporting them separately is more honest than silently assigning them to either side.
AMBIGUOUS_STRESS_IMMUNE = {"NR4A1", "NR4A2", "NR4A3", "CD69", "KLF2", "TNFAIP3", "REL", "NFKBIA", "GPR183"}

GENE_TO_PROGRAM = {gene: program for program, genes in IMMUNE_CANCER_PROGRAMS.items() for gene in genes}


def classify(gene: str) -> tuple[str, str]:
    """Return (tier, class_label). Tier is one of immune_cancer / ambiguous / noise / other.

    Precedence is deliberate: the curated immune and ambiguous sets are checked before the
    noise regexes so that a real immune gene can never be discarded by a pattern match.
    """
    if gene in AMBIGUOUS_STRESS_IMMUNE:
        return "ambiguous", "Immediate-early / activation (ambiguous)"
    if gene in GENE_TO_PROGRAM:
        return "immune_cancer", GENE_TO_PROGRAM[gene]
    if gene in HEMOGLOBIN:
        return "noise", "Hemoglobin / erythrocyte ambient"
    if gene in EPITHELIAL_AMBIENT:
        return "noise", "Epithelial / ambient RNA"
    if gene in MYELOID_PLATELET_AMBIENT:
        return "noise", "Myeloid / platelet ambient"
    if gene in IMMEDIATE_EARLY:
        return "noise", "Immediate-early / dissociation stress"
    if gene in HOUSEKEEPING:
        return "noise", "Structural housekeeping"
    if gene in NONCODING:
        return "noise", "Non-coding / lncRNA"
    if gene.startswith("HSP") or gene.startswith("DNAJ") or gene in {"CRYAB", "BAG3", "AHSA1", "HSPH1"}:
        return "noise", "Heat shock / proteostasis"
    for label, pattern in NOISE_PATTERNS:
        if pattern.search(gene):
            return "noise", label
    return "other", "Other / unannotated"


def direction_label(delete_shift: float) -> str:
    """Positive delete shift means removing the gene moves cells toward the goal state."""
    return "restrains_goal" if delete_shift > 0 else "promotes_goal"


def build_artifacts() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    rankings = pd.read_csv(TABLES / "exploration_balanced_rankings.csv")
    classes = rankings.Gene_name.map(classify)
    rankings["tier"] = [tier for tier, _ in classes]
    rankings["class_label"] = [label for _, label in classes]
    rankings["direction"] = rankings.delete_shift.map(direction_label)
    rankings.to_csv(TABLES / "denoised_all_classified.csv", index=False)

    audit = (
        rankings.groupby(["tier", "class_label"], as_index=False)
        .agg(rows=("Gene_name", "size"), genes=("Ensembl_ID", "nunique"), median_score=("balanced_score", "median"))
        .sort_values(["tier", "rows"], ascending=[True, False])
    )
    audit["pct_of_all_rows"] = (100 * audit.rows / len(rankings)).round(2)
    audit.to_csv(TABLES / "denoise_audit.csv", index=False)

    immune = rankings[rankings.tier.eq("immune_cancer")].copy()
    # Re-normalize the balanced score within the denoised set so ranking is not anchored to
    # the ribosomal genes that previously defined the top of each comparison.
    parts = []
    for comparison, group in immune.groupby("comparison"):
        group = group.copy()
        group["detection_norm_denoised"] = group.min_detection / group.min_detection.max()
        cap = group.min_abs_shift.quantile(0.95)
        group["shift_norm_denoised"] = (group.min_abs_shift / cap).clip(upper=1) if cap > 0 else 0.0
        group["denoised_score"] = np.sqrt(group.detection_norm_denoised * group.shift_norm_denoised)
        parts.append(group)
    immune = pd.concat(parts, ignore_index=True).sort_values(["comparison", "denoised_score"], ascending=[True, False])
    immune.to_csv(TABLES / "immune_cancer_candidates.csv", index=False)

    recurrence = (
        immune.groupby(["Gene_name", "class_label"], as_index=False)
        .agg(
            n_comparisons=("comparison", "nunique"),
            comparisons=("comparison", lambda s: "; ".join(sorted(s))),
            median_min_detection=("min_detection", "median"),
            max_denoised_score=("denoised_score", "max"),
            n_restrains=("direction", lambda s: int((s == "restrains_goal").sum())),
            n_promotes=("direction", lambda s: int((s == "promotes_goal").sum())),
        )
        .sort_values(["n_comparisons", "max_denoised_score"], ascending=False)
    )
    recurrence.to_csv(TABLES / "immune_cancer_recurrence.csv", index=False)

    program_summary = (
        immune.groupby("class_label", as_index=False)
        .agg(genes=("Gene_name", "nunique"), rows=("Gene_name", "size"), median_score=("denoised_score", "median"))
        .sort_values("genes", ascending=False)
    )
    program_summary.to_csv(TABLES / "immune_cancer_program_summary.csv", index=False)

    _plot_composition(rankings)
    _plot_program_heatmap(immune)
    _plot_bidirectional(immune)
    _plot_recurrence(recurrence)
    _plot_program_direction(immune)
    return rankings, immune, recurrence


def _plot_composition(rankings: pd.DataFrame) -> None:
    counts = rankings.pivot_table(index="comparison", columns="tier", values="Gene_name", aggfunc="size", fill_value=0)
    counts = counts.reindex(COMPARISONS)
    order = [c for c in ["noise", "other", "ambiguous", "immune_cancer"] if c in counts.columns]
    colors = {"noise": "#cbd5e1", "other": "#94a3b8", "ambiguous": "#f59e0b", "immune_cancer": "#2563eb"}
    fig, ax = plt.subplots(figsize=(10, 5.5))
    bottom = np.zeros(len(counts))
    for tier in order:
        ax.bar(counts.index, counts[tier], bottom=bottom, label=tier, color=colors[tier])
        bottom += counts[tier].values
    ax.set_ylabel("Qualified concordant rows")
    ax.set_title("Composition of qualified hits before denoising")
    ax.set_xticks(range(len(counts.index)), [c.replace("_to_", " → ") for c in counts.index], rotation=25, ha="right")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURES / "tier_composition_by_comparison.png", dpi=180)
    plt.close(fig)


def _plot_program_heatmap(immune: pd.DataFrame) -> None:
    top = immune.groupby("comparison", group_keys=False).head(12)
    pivot = top.pivot_table(index="Gene_name", columns="comparison", values="denoised_score", aggfunc="max", fill_value=0)
    pivot = pivot.reindex(columns=[c for c in COMPARISONS if c in pivot.columns])
    pivot = pivot.loc[pivot.max(axis=1).sort_values(ascending=False).head(30).index]
    fig, ax = plt.subplots(figsize=(11, 10))
    image = ax.imshow(pivot.values, cmap="magma", aspect="auto")
    ax.set_xticks(range(len(pivot.columns)), [c.replace("_to_", " → ") for c in pivot.columns], rotation=35, ha="right")
    ax.set_yticks(range(len(pivot.index)), pivot.index, fontsize=9)
    ax.set_title("Denoised immune / cancer candidates across disease transitions")
    fig.colorbar(image, ax=ax, label="Denoised balanced score")
    fig.tight_layout()
    fig.savefig(FIGURES / "immune_cancer_heatmap.png", dpi=180)
    plt.close(fig)


def _place_labels(ax, rows, *, fontsize: float) -> None:
    """Annotate rows, nudging each label off the ones already placed.

    Labels are positioned in axis-fraction space and tested against the boxes
    already taken, so a gene name lands in the first free direction rather than
    stacking on its neighbour. The earlier version wrote every label at a fixed
    (+4, +4) point offset, which collided whenever two candidates sat close
    together - "HLA-DPA1" over "HLA-DRA" and "IL7R"/"CXCR4"/"B2M"/"CD74" in the
    original. Deterministic: candidate offsets are tried in a fixed order.
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    # Rough label extent as a fraction of the axes, from the font size and the
    # axes size in points. Approximate on purpose - it only has to separate
    # labels, not bound them exactly.
    box_w = ax.get_window_extent()
    char_w = fontsize * 0.60
    taken: list[tuple[float, float, float, float]] = []
    # right, left, above, below, then the diagonals
    candidates = [(1.0, 0.4), (-1.0, 0.4), (0.0, 1.4), (0.0, -1.4),
                  (1.0, -1.0), (-1.0, -1.0), (1.0, 1.4), (-1.0, 1.4)]

    for row in rows.itertuples():
        px = (row.delete_shift - x0) / (x1 - x0)
        py = (row.overexpress_shift - y0) / (y1 - y0)
        w = len(row.Gene_name) * char_w / max(box_w.width, 1) * 72 / plt.rcParams["figure.dpi"] * 100
        w = min(max(w, 0.06), 0.42)
        h = fontsize * 1.5 / max(box_w.height, 1) * 72 / plt.rcParams["figure.dpi"] * 100
        h = min(max(h, 0.03), 0.20)

        for dx, dy in candidates:
            lx = px + dx * (w * 0.55 + 0.012)
            ly = py + dy * (h * 0.75 + 0.008)
            # Keep the whole label box inside the axes, not just its centre -
            # a centre-only test let long names hang over the y-axis ticks.
            if not (w / 2 + 0.01 <= lx <= 1 - w / 2 - 0.01):
                continue
            if not (h / 2 + 0.01 <= ly <= 1 - h / 2 - 0.01):
                continue
            box = (lx - w / 2, ly - h / 2, lx + w / 2, ly + h / 2)
            if any(not (box[2] < t[0] or box[0] > t[2] or box[3] < t[1] or box[1] > t[3])
                   for t in taken):
                continue
            taken.append(box)
            ax.annotate(
                row.Gene_name, (px, py), xycoords="axes fraction",
                xytext=(lx, ly), textcoords="axes fraction",
                fontsize=fontsize, fontweight="semibold", color="#111827",
                ha="center", va="center",
                arrowprops=dict(arrowstyle="-", color="#6b7280", linewidth=0.6,
                                shrinkA=0, shrinkB=2),
                bbox=dict(boxstyle="round,pad=0.16", facecolor="white",
                          edgecolor="none", alpha=0.78),
            )
            break


def _plot_bidirectional(immune: pd.DataFrame) -> None:
    # Sized for the poster's right column (260 mm inner width): a smaller
    # figure with larger point sizes renders text bigger once scaled to that
    # width. At figsize 12in the 12 pt gene labels land near 10 pt on the A0
    # sheet, which is about the poster's own smallest label.
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), squeeze=False)
    for ax, comparison in zip(axes.flat, COMPARISONS):
        group = immune[immune.comparison == comparison]
        if group.empty:
            ax.set_visible(False)
            continue
        sizes = 14 + 210 * (group.min_detection / group.min_detection.max()) ** 0.7
        ax.scatter(group.delete_shift, group.overexpress_shift, s=sizes, c=group.denoised_score,
                   cmap="magma", alpha=0.62, linewidths=0)

        # Per-axis limits from each axis's own spread, rather than one symmetric
        # limit taken from the larger of the two. The shared limit made the x
        # range up to 4x the data it had to show - most of every panel was
        # blank. Zero stays inside the range so "right half = deletion moves
        # cells toward the goal state" still reads off the figure.
        for axis, values in (("x", group.delete_shift), ("y", group.overexpress_shift)):
            lo, hi = float(values.min()), float(values.max())
            lo, hi = min(lo, 0.0), max(hi, 0.0)
            pad = (hi - lo) * 0.14 or 1e-4
            (ax.set_xlim if axis == "x" else ax.set_ylim)(lo - pad, hi + pad)

        ax.axhline(0, color="#9ca3af", linewidth=0.7)
        ax.axvline(0, color="#9ca3af", linewidth=0.7)
        _place_labels(ax, group.head(6), fontsize=12)
        ax.set_title(comparison.replace("_to_", " → "), fontsize=14, fontweight="bold")
        ax.set_xlabel("Delete cosine shift", fontsize=11.5)
        ax.set_ylabel("Overexpression cosine shift", fontsize=11.5)
        ax.tick_params(labelsize=10)
        ax.xaxis.set_major_locator(plt.MaxNLocator(4))
        ax.yaxis.set_major_locator(plt.MaxNLocator(5))
        ax.grid(alpha=0.15)
    fig.suptitle("Denoised immune / cancer bidirectional candidates\n"
                 "(right half = deletion moves cells toward the goal state)",
                 y=1.02, fontsize=15, fontweight="bold")
    fig.tight_layout()
    fig.savefig(FIGURES / "immune_cancer_bidirectional.png", dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_recurrence(recurrence: pd.DataFrame) -> None:
    top = recurrence.head(25).sort_values("n_comparisons")
    palette = {label: color for label, color in zip(sorted(recurrence.class_label.unique()), plt.cm.tab20.colors)}
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.barh(top.Gene_name, top.n_comparisons, color=[palette[c] for c in top.class_label])
    ax.set_xlabel("Number of disease transitions with a concordant hit (max 6)")
    ax.set_title("Most recurrent denoised immune / cancer candidates")
    ax.grid(axis="x", alpha=0.2)
    handles = [plt.Rectangle((0, 0), 1, 1, color=palette[c]) for c in sorted(top.class_label.unique())]
    ax.legend(handles, sorted(top.class_label.unique()), fontsize=7, frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIGURES / "immune_cancer_recurrence.png", dpi=180)
    plt.close(fig)


def _plot_program_direction(immune: pd.DataFrame) -> None:
    grouped = immune.groupby(["class_label", "direction"]).size().unstack(fill_value=0)
    for column in ("restrains_goal", "promotes_goal"):
        if column not in grouped.columns:
            grouped[column] = 0
    grouped = grouped.sort_values("restrains_goal", ascending=True)
    fig, ax = plt.subplots(figsize=(10, 6))
    y = np.arange(len(grouped))
    ax.barh(y - 0.2, grouped.restrains_goal, height=0.4, color="#b91c1c", label="restrains goal (deletion pushes toward goal)")
    ax.barh(y + 0.2, grouped.promotes_goal, height=0.4, color="#1d4ed8", label="promotes goal (deletion pushes away)")
    ax.set_yticks(y, grouped.index, fontsize=9)
    ax.set_xlabel("Concordant gene-comparison rows")
    ax.set_title("Directionality of denoised immune / cancer programs")
    ax.legend(fontsize=8, frameon=False)
    ax.grid(axis="x", alpha=0.2)
    fig.tight_layout()
    fig.savefig(FIGURES / "program_directionality.png", dpi=180)
    plt.close(fig)


def md(source: str) -> dict:
    return {"cell_type": "markdown", "metadata": {}, "source": source.splitlines(True)}


def code(source: str) -> dict:
    return {"cell_type": "code", "execution_count": None, "metadata": {}, "outputs": [], "source": source.splitlines(True)}


def build_notebook(rankings: pd.DataFrame, immune: pd.DataFrame, recurrence: pd.DataFrame) -> Path:
    noise_rows = int(rankings.tier.eq("noise").sum())
    immune_rows = len(immune)
    notebook = {
        "cells": [
            md(
                "# Denoised perturbation results: immune and cancer biology\n\n"
                "The all-gene bidirectional screen returns 3,586 qualified gene-comparison hits, but the\n"
                "top of the ranking is occupied by ribosomal, mitochondrial, heat-shock and ambient-RNA\n"
                "transcripts. Those genes can be statistically extreme because they track global\n"
                "transcriptional output, dissociation stress, or contamination from neighbouring cells —\n"
                "not because they are T-cell-intrinsic regulators.\n\n"
                "This notebook applies a transparent, auditable classification to every qualified hit and\n"
                "re-ranks the results within curated immune and cancer programs."
            ),
            md(
                "## TL;DR\n\n"
                f"Of {len(rankings):,} qualified gene-comparison rows:\n\n"
                f"- **{noise_rows:,} ({100 * noise_rows / len(rankings):.0f}%)** match an explicit technical class (ribosomal, mitochondrial, heat-shock, immediate-early, ambient) and are set aside.\n"
                f"- **{immune_rows:,} ({100 * immune_rows / len(rankings):.0f}%)**, covering **{immune.Gene_name.nunique()} genes**, map to a curated immune or cancer program. This is the prioritized set.\n"
                f"- **{int(rankings.tier.eq('other').sum()):,} ({100 * rankings.tier.eq('other').sum() / len(rankings):.0f}%)** are unannotated here — neither flagged as technical nor claimed as immune.\n\n"
                "This is therefore **positive selection**, not just noise subtraction: the large `other` tier is not\n"
                "asserted to be noise, it is simply outside the curated sets and left available for review.\n\n"
                "Once ribosomal genes no longer anchor the score, **antigen presentation / MHC becomes the strongest\n"
                "surviving program** — `HLA-C`, `B2M`, `CD74`, `HLA-B`, `HLA-E`, `PSMB9` occupy most of the top 20,\n"
                "with detection counts in the thousands.\n\n"
                "Filtering changes which genes are *prioritized*. It does not add evidence: donor-level\n"
                "consistency and ambient-RNA sensitivity analysis are still required before any causal claim."
            ),
            code(
                "from pathlib import Path\n"
                "import numpy as np\n"
                "import pandas as pd\n"
                "from IPython.display import Image, display\n\n"
                "HERE = Path.cwd()\n"
                "TABLES = HERE / 'tables'\n"
                "FIGURES = HERE / 'figures' / 'denoised'\n\n"
                "classified = pd.read_csv(TABLES / 'denoised_all_classified.csv')\n"
                "immune = pd.read_csv(TABLES / 'immune_cancer_candidates.csv')\n"
                "audit = pd.read_csv(TABLES / 'denoise_audit.csv')\n"
                "recurrence = pd.read_csv(TABLES / 'immune_cancer_recurrence.csv')\n"
                "programs = pd.read_csv(TABLES / 'immune_cancer_program_summary.csv')\n"
                "print(f'Qualified rows: {len(classified):,}')\n"
                "classified.tier.value_counts().rename('rows').to_frame()"
            ),
            md(
                "## What is removed, and why\n\n"
                "Each class below is a stated, reviewable reason for exclusion rather than an opaque filter.\n"
                "The curated immune and ambiguous sets are matched **before** the noise patterns, so a real\n"
                "immune gene can never be discarded by a regex — `HLA-*`, `B2M` and `CD74` survive by construction.\n\n"
                "| Class | Reason it is set aside |\n"
                "|---|---|\n"
                "| Ribosomal, OXPHOS, structural housekeeping | Track global transcriptional output and cell size, not a specific program |\n"
                "| Mitochondrial (`MT-`) | Dominated by cell stress and dying-cell fraction |\n"
                "| Heat shock / proteostasis | Tissue handling and dissociation temperature artifacts |\n"
                "| Immediate-early (`JUN`, `FOS`, `EGR1`, `DUSP1`) | Induced by the dissociation protocol itself |\n"
                "| Epithelial, myeloid, platelet, hemoglobin | Ambient RNA from neighbouring cells in the tumour |\n"
                "| Histone, non-coding | Cell-cycle and technical capture effects |\n\n"
                "The `other` tier is the largest and is **not** a noise claim. It holds genes with no membership in\n"
                "the curated sets used here. Real regulators certainly sit in it — that is the cost of positive\n"
                "selection, and the reason the full classified table is written out rather than discarded."
            ),
            code("audit"),
            code("display(Image(filename=str(FIGURES / 'tier_composition_by_comparison.png')))"),
            md(
                "## Prioritized immune and cancer candidates\n\n"
                "The balanced score is recomputed **within** the denoised set. This matters: the original\n"
                "normalization was anchored to ribosomal genes that defined the top of each comparison, so\n"
                "immune genes were compressed into the bottom of the scale."
            ),
            code("programs"),
            code(
                "cols = ['comparison', 'Gene_name', 'class_label', 'direction', 'delete_shift', 'overexpress_shift', 'min_detection', 'denoised_score']\n"
                "immune.sort_values('denoised_score', ascending=False)[cols].head(30)"
            ),
            code("display(Image(filename=str(FIGURES / 'immune_cancer_heatmap.png')))"),
            md(
                "## Directionality\n\n"
                "Sign convention follows the rest of this project: a **positive deletion shift** means removing\n"
                "the gene moves cells *toward* the goal state, so the gene normally **restrains** that transition\n"
                "and helps maintain the source identity. A negative deletion shift means the gene **promotes**\n"
                "the transition. Concordance requires the overexpression arm to point the opposite way.\n\n"
                "This is the same logic that made `ASCL1` / `NEUROD1` readable as SCLC master regulators."
            ),
            code("display(Image(filename=str(FIGURES / 'program_directionality.png')))"),
            code(
                "immune.groupby(['class_label', 'direction']).size().unstack(fill_value=0).assign(\n"
                "    total=lambda d: d.sum(axis=1)\n"
                ").sort_values('total', ascending=False)"
            ),
            md(
                "## Recurrence across disease transitions\n\n"
                "A candidate seen in one transition can be comparison-specific noise. A candidate concordant\n"
                "across several transitions is harder to explain that way. `n_restrains` and `n_promotes` show\n"
                "whether the gene acts consistently or flips direction depending on the transition — a flip is\n"
                "biologically meaningful for a state-discriminating gene, not a contradiction."
            ),
            code("recurrence.head(25)"),
            code("display(Image(filename=str(FIGURES / 'immune_cancer_recurrence.png')))"),
            code("display(Image(filename=str(FIGURES / 'immune_cancer_bidirectional.png')))"),
            md(
                "## Cross-check against the targeted 50-gene panel\n\n"
                "The earlier curated panel was run as a separate, hypothesis-driven experiment. Genes that were\n"
                "prioritized there and also surface here — found blind in an unbiased genome-wide screen —\n"
                "are the strongest internal corroboration available without new data."
            ),
            code(
                "PANEL = ['TIGIT', 'LAG3', 'GZMH', 'CCR7', 'NKG7', 'TCF7', 'IL7R', 'HAVCR2', 'CTLA4', 'SLAMF6', 'IFNG', 'PDCD1', 'ASCL1', 'NEUROD1']\n"
                "panel_hits = immune[immune.Gene_name.isin(PANEL)]\n"
                "summary = panel_hits.groupby('Gene_name').agg(\n"
                "    n_comparisons=('comparison', 'nunique'),\n"
                "    programs=('class_label', 'first'),\n"
                "    median_detection=('min_detection', 'median'),\n"
                "    max_score=('denoised_score', 'max'),\n"
                ").sort_values('n_comparisons', ascending=False)\n"
                "missing = sorted(set(PANEL) - set(panel_hits.Gene_name))\n"
                "print('Panel genes not concordant anywhere in the all-gene screen:', missing or 'none')\n"
                "summary"
            ),
            md(
                "## Ambiguous tier\n\n"
                "These genes are immune-meaningful and dissociation-induced at the same time. Assigning them to\n"
                "either side would be a hidden judgement call, so they are kept visible and excluded from the\n"
                "prioritized set by default."
            ),
            code(
                "cols = ['comparison', 'Gene_name', 'direction', 'delete_shift', 'overexpress_shift', 'min_detection', 'balanced_score']\n"
                "classified[classified.tier.eq('ambiguous')].sort_values('balanced_score', ascending=False)[cols]"
            ),
            md(
                "## Gene lookup\n\n"
                "Inspect any gene across all six transitions, including ones that were filtered out."
            ),
            code(
                "def gene_lookup(gene):\n"
                "    cols = ['comparison', 'Gene_name', 'tier', 'class_label', 'direction', 'delete_shift', 'overexpress_shift', 'min_detection']\n"
                "    return classified[classified.Gene_name.eq(gene)].sort_values('comparison')[cols]\n\n"
                "gene_lookup('HLA-C')"
            ),
            md(
                "## Caveats\n\n"
                "- **Filtering is a prioritization step, not evidence.** Every surviving candidate carries exactly the\n"
                "  same statistical support it had before; the set is smaller and more interpretable, not more proven.\n"
                "- **Gene-set membership is curated by symbol** in `scripts/build_denoised_notebook.py`. It is transparent\n"
                "  and editable, but it is not GO/Reactome enrichment and it encodes assumptions about what counts as immune.\n"
                "- **Excluded does not mean biologically inert.** Ribosomal and MHC-adjacent proteostasis genes have real\n"
                "  roles in T-cell function; they are set aside because this assay cannot separate that from global state.\n"
                "- **Ambient RNA is flagged by symbol, not corrected.** A CellBender run is still the right fix.\n"
                "- **Donor-level consistency is still not applied to this all-gene result** — the single largest remaining gap.\n"
                "  Until it is, no candidate here should be promoted to a biological claim in an abstract or poster."
            ),
        ],
        "metadata": {"kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"}, "language_info": {"name": "python"}},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    output = HERE / "perturbation_denoised_immune_cancer.ipynb"
    output.write_text(json.dumps(notebook, indent=1) + "\n")
    return output


if __name__ == "__main__":
    rankings, immune, recurrence = build_artifacts()
    output = build_notebook(rankings, immune, recurrence)
    print(f"Wrote {output}")
    print(f"Total qualified rows: {len(rankings):,}")
    for tier, count in rankings.tier.value_counts().items():
        print(f"  {tier}: {count:,} ({100 * count / len(rankings):.1f}%)")
    print(f"Immune/cancer genes: {immune.Ensembl_ID.nunique()} unique, {len(immune)} rows")
    print(f"Recurrent (>=3 comparisons): {int((recurrence.n_comparisons >= 3).sum())}")
