#!/usr/bin/env python3
"""Build reproducible perturbation summaries and pathway-enrichment inputs."""

from __future__ import annotations

import json
import gzip
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests


HERE = Path(__file__).resolve().parent
DEFAULT_STATS = Path(
    "/home/petadimensionlab/workspace/Geneformer/KD/"
    "tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/stats"
)
STATS_DIR = Path(os.environ.get("PERTURBATION_STATS_DIR", DEFAULT_STATS))
TABLE_DIR = HERE / "source_tables"
FIGURE_DIR = HERE / "figures"
RAW_ENRICHMENT_DIR = HERE / "enrichment_raw"

FDR_THRESHOLD = 0.05
MIN_DETECTIONS = 25
GPROFILER_URL = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
GPROFILER_SOURCES = ("REAC", "KEGG", "GO:BP")

COMPARISON_ORDER = (
    "lusc_to_luad",
    "lusc_to_normal",
    "luad_to_lusc",
    "luad_to_normal",
    "normal_to_luad",
    "normal_to_lusc",
)
DISPLAY = {
    "lusc_to_luad": "LUSC → LUAD",
    "lusc_to_normal": "LUSC → NORMAL",
    "luad_to_lusc": "LUAD → LUSC",
    "luad_to_normal": "LUAD → NORMAL",
    "normal_to_luad": "NORMAL → LUAD",
    "normal_to_lusc": "NORMAL → LUSC",
}


def load_comparisons() -> dict[str, pd.DataFrame]:
    comparisons: dict[str, pd.DataFrame] = {}
    for comparison in COMPARISON_ORDER:
        path = STATS_DIR / f"heldout_allgene_{comparison}.csv"
        frame = pd.read_csv(path).drop(columns=["Unnamed: 0"], errors="ignore")
        frame["comparison"] = comparison
        frame["comparison_label"] = DISPLAY[comparison]
        frame["direction"] = np.where(frame["Shift_to_goal_end"] > 0, "toward goal", "away from goal")
        frame["qualified_goal_shift"] = (
            (frame["Goal_end_FDR"] < FDR_THRESHOLD)
            & (frame["Shift_to_goal_end"] > 0)
            & (frame["N_Detections"] >= MIN_DETECTIONS)
        )
        comparisons[comparison] = frame
    return comparisons


def summarize(comparisons: dict[str, pd.DataFrame]) -> pd.DataFrame:
    rows = []
    for comparison, frame in comparisons.items():
        significant = frame[frame["Goal_end_FDR"] < FDR_THRESHOLD]
        qualified = frame[frame["qualified_goal_shift"]]
        rows.append(
            {
                "comparison": comparison,
                "comparison_label": DISPLAY[comparison],
                "genes_tested": len(frame),
                "fdr_significant": len(significant),
                "significant_toward_goal": int((significant["Shift_to_goal_end"] > 0).sum()),
                "significant_away_from_goal": int((significant["Shift_to_goal_end"] < 0).sum()),
                "qualified_toward_goal": len(qualified),
                "qualified_fraction": len(qualified) / len(frame),
                "max_positive_shift": qualified["Shift_to_goal_end"].max() if len(qualified) else np.nan,
                "median_qualified_shift": qualified["Shift_to_goal_end"].median() if len(qualified) else np.nan,
                "median_qualified_detections": qualified["N_Detections"].median() if len(qualified) else np.nan,
            }
        )
    return pd.DataFrame(rows)


def top_gene_table(comparisons: dict[str, pd.DataFrame], limit: int = 15) -> pd.DataFrame:
    frames = []
    for comparison, frame in comparisons.items():
        qualified = frame[frame["qualified_goal_shift"]].copy()
        qualified = qualified.sort_values(
            ["Shift_to_goal_end", "Goal_end_FDR"], ascending=[False, True]
        ).head(limit)
        qualified["rank"] = np.arange(1, len(qualified) + 1)
        frames.append(
            qualified[
                [
                    "comparison",
                    "comparison_label",
                    "rank",
                    "Gene_name",
                    "Ensembl_ID",
                    "Shift_to_goal_end",
                    "Goal_end_FDR",
                    "N_Detections",
                ]
            ]
        )
    return pd.concat(frames, ignore_index=True)


def gprofiler_enrichment(
    comparison: str, query_genes: list[str], background_genes: list[str]
) -> tuple[pd.DataFrame, dict]:
    payload = {
        "organism": "hsapiens",
        "query": query_genes,
        "sources": list(GPROFILER_SOURCES),
        "user_threshold": 0.05,
        "significance_threshold_method": "fdr",
        "domain_scope": "custom",
        "background": background_genes,
        "no_evidences": False,
    }
    response = requests.post(GPROFILER_URL, json=payload, timeout=120)
    response.raise_for_status()
    raw = response.json()
    rows = []
    for result in raw.get("result", []):
        evidence = result.get("intersections", [])
        intersection_genes = [
            gene for gene, annotations in zip(query_genes, evidence) if annotations
        ]
        rows.append(
            {
                "comparison": comparison,
                "comparison_label": DISPLAY[comparison],
                "source": result["source"],
                "term_id": result["native"],
                "term_name": result["name"],
                "adjusted_p_value": result["p_value"],
                "minus_log10_p": -np.log10(max(result["p_value"], np.finfo(float).tiny)),
                "term_size": result["term_size"],
                "query_size": result["query_size"],
                "intersection_size": result["intersection_size"],
                "precision": result["precision"],
                "recall": result["recall"],
                "effective_domain_size": result["effective_domain_size"],
                "intersection_genes": ";".join(intersection_genes),
            }
        )
    return pd.DataFrame(rows), raw


def pathway_analysis(comparisons: dict[str, pd.DataFrame]) -> pd.DataFrame:
    frames = []
    RAW_ENRICHMENT_DIR.mkdir(parents=True, exist_ok=True)
    for comparison, frame in comparisons.items():
        qualified = frame[frame["qualified_goal_shift"]]
        query = sorted(qualified["Gene_name"].dropna().astype(str).unique())
        background = sorted(frame["Gene_name"].dropna().astype(str).unique())
        enriched, raw = gprofiler_enrichment(comparison, query, background)
        with gzip.open(
            RAW_ENRICHMENT_DIR / f"{comparison}.json.gz", "wt", encoding="utf-8"
        ) as handle:
            json.dump(raw, handle, separators=(",", ":"))
        if not enriched.empty:
            frames.append(enriched)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True).sort_values(
        ["comparison", "source", "adjusted_p_value"]
    )


PATHWAY_THEMES = {
    "Translation / ribosome": r"translation|ribosom|peptide chain|protein biosynthetic|selenocysteine|srp-dependent",
    "Immune / cytokine": r"immune|cytokine|interferon|lymphocyte|t cell|innate|defense response",
    "Oxidative phosphorylation": r"oxidative phosphorylation|electron transport|respiratory chain|atp synthesis",
    "IL-17 / inflammation": r"il-17|interleukin-17|inflamm|neutrophil",
    "Antigen processing": r"antigen processing|antigen presentation|major histocompatibility|mhc",
    "Cellular stress": r"response to stress|cellular stress|heat shock",
}


def summarize_pathway_themes(pathways: pd.DataFrame) -> pd.DataFrame:
    """Consolidate redundant enrichment terms into transparent, predefined themes."""
    rows = []
    for comparison in COMPARISON_ORDER:
        subset = pathways[pathways["comparison"] == comparison]
        for theme, pattern in PATHWAY_THEMES.items():
            matches = subset[
                subset["term_name"].str.lower().str.contains(pattern, regex=True, na=False)
            ]
            if matches.empty:
                continue
            best = matches.sort_values("adjusted_p_value").iloc[0]
            rows.append(
                {
                    "comparison": comparison,
                    "comparison_label": DISPLAY[comparison],
                    "theme": theme,
                    "chart_label": f"{DISPLAY[comparison]} · {theme}",
                    "best_adjusted_p_value": best["adjusted_p_value"],
                    "minus_log10_p": best["minus_log10_p"],
                    "matching_term_count": len(matches),
                    "representative_source": best["source"],
                    "representative_term": best["term_name"],
                    "intersection_size": best["intersection_size"],
                }
            )
    return pd.DataFrame(rows).sort_values(["comparison", "minus_log10_p"], ascending=[True, False])


def plot_goal_shifts(top_genes: pd.DataFrame) -> None:
    fig, axes = plt.subplots(3, 2, figsize=(13, 15), constrained_layout=True)
    for axis, comparison in zip(axes.flat, COMPARISON_ORDER):
        plot = top_genes[top_genes["comparison"] == comparison].sort_values(
            "Shift_to_goal_end"
        )
        axis.barh(plot["Gene_name"], plot["Shift_to_goal_end"], color="#2867B2")
        axis.axvline(0, color="#2F3941", linewidth=0.8)
        axis.set_title(DISPLAY[comparison], loc="left", fontweight="bold")
        axis.set_xlabel("Mean cosine shift toward goal centroid")
        axis.grid(axis="x", color="#E3E8EC", linewidth=0.8)
        axis.spines[["top", "right"]].set_visible(False)
    fig.suptitle(
        "Coverage-qualified positive goal shifts",
        x=0.01,
        ha="left",
        fontsize=18,
        fontweight="bold",
    )
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / "goal_shift_top_genes.png", dpi=180, facecolor="white")
    plt.close(fig)


def plot_pathways(pathways: pd.DataFrame) -> None:
    if pathways.empty:
        return
    top = (
        summarize_pathway_themes(pathways)
        .sort_values("best_adjusted_p_value")
        .groupby("comparison", sort=False)
        .head(5)
        .copy()
    )
    top["label"] = top["theme"] + " · " + top["comparison_label"]
    top = top.sort_values("minus_log10_p")
    fig, axis = plt.subplots(figsize=(13, 12), constrained_layout=True)
    axis.barh(top["label"], top["minus_log10_p"], color="#B08900")
    axis.set_xlabel("−log10 adjusted enrichment p-value")
    axis.set_title("Enriched biological themes across directional goal shifts", loc="left", fontweight="bold")
    axis.grid(axis="x", color="#E3E8EC", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_DIR / "pathway_enrichment.png", dpi=180, facecolor="white")
    plt.close(fig)


def main() -> None:
    TABLE_DIR.mkdir(parents=True, exist_ok=True)
    comparisons = load_comparisons()
    summary = summarize(comparisons)
    top_genes = top_gene_table(comparisons)
    pathways = pathway_analysis(comparisons)
    pathway_themes = summarize_pathway_themes(pathways)

    summary.to_csv(TABLE_DIR / "comparison_summary.csv", index=False)
    top_genes.to_csv(TABLE_DIR / "top_goal_shift_genes.csv", index=False)
    pathways.to_csv(TABLE_DIR / "pathway_enrichment.csv", index=False)
    pathway_themes.to_csv(TABLE_DIR / "pathway_theme_summary.csv", index=False)
    plot_goal_shifts(top_genes)
    plot_pathways(pathways)

    metadata = {
        "stats_directory": str(STATS_DIR),
        "comparison_order": list(COMPARISON_ORDER),
        "goal_shift_filter": {
            "goal_end_fdr_lt": FDR_THRESHOLD,
            "shift_to_goal_end_gt": 0,
            "n_detections_gte": MIN_DETECTIONS,
        },
        "pathway_method": {
            "service": "g:Profiler g:GOSt",
            "organism": "hsapiens",
            "sources": list(GPROFILER_SOURCES),
            "domain_scope": "custom",
            "background": "all genes tested in the corresponding source-state table",
            "multiple_testing": "Benjamini-Hochberg FDR",
        },
    }
    (HERE / "analysis_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    print(summary.to_string(index=False))
    print(f"\nPathway rows: {len(pathways):,}")


if __name__ == "__main__":
    main()
