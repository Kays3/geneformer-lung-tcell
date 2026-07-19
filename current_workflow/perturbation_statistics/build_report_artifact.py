#!/usr/bin/env python3
"""Build the canonical portable-report artifact from reviewed analysis tables."""

from __future__ import annotations

import json
from datetime import datetime
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
TABLES = HERE / "source_tables"
GENERATED_AT = datetime.now().astimezone().isoformat(timespec="seconds")


def records(frame: pd.DataFrame) -> list[dict]:
    return json.loads(frame.to_json(orient="records"))


summary = pd.read_csv(TABLES / "comparison_summary.csv")
genes = pd.read_csv(TABLES / "top_goal_shift_genes.csv")
themes = pd.read_csv(TABLES / "pathway_theme_summary.csv")

summary["qualified_percent"] = summary["qualified_fraction"] * 100
summary["fdr_significant_percent"] = summary["fdr_significant"] / summary["genes_tested"] * 100
overview = pd.DataFrame(
    [
        {
            "completed_comparisons": 6,
            "total_comparisons": 6,
            "cells_completed": 3379,
            "cells_total": 3379,
            "deletions_completed": 348313 + 1150097 + 1439366,
            "qualified_gene_shifts": int(summary["qualified_toward_goal"].sum()),
            "fdr_significant_gene_tests": int(summary["fdr_significant"].sum()),
        }
    ]
)
overview.to_csv(TABLES / "overview_summary.csv", index=False)

sources = [
    {
        "id": "overview",
        "label": "Completed run overview",
        "path": "source_tables/overview_summary.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('source_tables/overview_summary.csv');",
            "description": "Run totals consolidated from the completed live progress snapshot and six statistical tables.",
            "metric_definitions": [
                "Deletion observations = 348,313 LUSC + 1,150,097 LUAD + 1,439,366 normal deletion observations.",
                "Qualified toward-goal tests = sum across six tables after FDR, positive-shift, and detection filters.",
            ],
        },
    },
    {
        "id": "comparison-summary",
        "label": "Directional perturbation comparison summary",
        "path": "source_tables/comparison_summary.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('source_tables/comparison_summary.csv');",
            "description": "Per-comparison counts from six Geneformer goal-state-shift result tables.",
            "filters": ["Goal_end_FDR < 0.05", "Shift_to_goal_end > 0", "N_Detections >= 25 for qualified shifts"],
            "metric_definitions": [
                "FDR-significant = Goal_end_FDR < 0.05, regardless of shift direction.",
                "Qualified toward-goal = Goal_end_FDR < 0.05 AND Shift_to_goal_end > 0 AND N_Detections >= 25.",
            ],
        },
    },
    {
        "id": "top-genes",
        "label": "Coverage-qualified goal-shift gene rankings",
        "path": "source_tables/top_goal_shift_genes.csv",
        "query": {
            "engine": "DuckDB",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('source_tables/top_goal_shift_genes.csv');",
            "description": "Top 15 positive goal shifts per comparison after the prespecified FDR, direction, and coverage filter.",
            "filters": ["Goal_end_FDR < 0.05", "Shift_to_goal_end > 0", "N_Detections >= 25"],
        },
    },
    {
        "id": "pathways",
        "label": "g:Profiler pathway enrichment and thematic consolidation",
        "path": "source_tables/pathway_theme_summary.csv",
        "href": "https://biit.cs.ut.ee/gprofiler/page/apis",
        "query": {
            "engine": "g:Profiler g:GOSt",
            "language": "sql",
            "sql": "SELECT * FROM read_csv_auto('source_tables/pathway_theme_summary.csv');",
            "description": "Over-representation analysis for each qualified gene set against all genes tested in that comparison; predefined regex themes consolidate redundant significant terms.",
            "filters": ["Organism: hsapiens", "Sources: Reactome, KEGG, GO Biological Process", "FDR < 0.05"],
            "metric_definitions": ["Bar length is -log10 of the best adjusted enrichment p-value within each predefined theme."],
        },
    },
    {
        "id": "geneformer-method",
        "label": "Geneformer transfer-learning and perturbation framework",
        "href": "https://www.nature.com/articles/s41586-023-06139-9",
    },
]

cards = [
    {
        "id": "cells-card",
        "dataset": "overview",
        "description": "All planned source cells completed perturbation.",
        "sourceId": "overview",
        "metrics": [{"label": "Cells completed", "field": "cells_completed", "format": "number"}],
    },
    {
        "id": "deletions-card",
        "dataset": "overview",
        "description": "Marker-deletion observations written across LUSC, LUAD, and normal sources.",
        "sourceId": "overview",
        "metrics": [{"label": "Deletion observations", "field": "deletions_completed", "format": "compact"}],
    },
    {
        "id": "comparisons-card",
        "dataset": "overview",
        "description": "Every directional source-to-goal statistical table is present.",
        "sourceId": "overview",
        "metrics": [{"label": "Comparisons complete", "field": "completed_comparisons", "format": "number"}],
    },
    {
        "id": "qualified-card",
        "dataset": "overview",
        "description": "Sum across comparisons; a gene can appear in more than one direction.",
        "sourceId": "overview",
        "metrics": [{"label": "Qualified toward-goal tests", "field": "qualified_gene_shifts", "format": "number"}],
    },
]

charts = [
    {
        "id": "direction-chart",
        "title": "FDR-significant shifts split almost evenly by direction",
        "subtitle": "Counts at Goal_end_FDR < 0.05; no coverage threshold",
        "intent": "composition",
        "question": "How many significant shifts point toward versus away from each goal centroid?",
        "rationale": "A stacked horizontal bar prevents the raw FDR count from being mistaken for positive goal shifts.",
        "type": "horizontalStackedBar",
        "dataset": "comparison_summary",
        "sourceId": "comparison-summary",
        "encodings": {
            "x": {"field": "comparison_label", "type": "nominal", "label": "Comparison"},
            "y": {"fields": ["significant_toward_goal", "significant_away_from_goal"], "type": "quantitative", "label": "Gene tests"},
            "tooltip": [
                {"field": "genes_tested", "type": "quantitative", "label": "Genes tested"},
                {"field": "qualified_toward_goal", "type": "quantitative", "label": "Coverage-qualified toward goal"},
            ],
        },
        "layout": "full",
        "maxRows": 6,
    },
    {
        "id": "pathway-chart",
        "title": "Translation dominates most directions; immune and respiratory programs add context",
        "subtitle": "Best adjusted enrichment p-value per predefined biological theme",
        "intent": "comparison",
        "question": "Which broad biological themes recur among qualified toward-goal deletions?",
        "rationale": "Theme consolidation reduces term redundancy while preserving the best term, source, intersection size, and number of matching terms for audit.",
        "type": "horizontalBar",
        "dataset": "pathway_themes",
        "sourceId": "pathways",
        "encodings": {
            "x": {"field": "chart_label", "type": "nominal", "label": "Comparison · theme"},
            "y": {"field": "minus_log10_p", "type": "quantitative", "label": "−log10 adjusted p-value"},
            "color": {"field": "theme", "type": "nominal", "label": "Theme"},
            "tooltip": [
                {"field": "representative_term", "type": "text", "label": "Representative term"},
                {"field": "representative_source", "type": "text", "label": "Database"},
                {"field": "intersection_size", "type": "quantitative", "label": "Intersection genes"},
                {"field": "matching_term_count", "type": "quantitative", "label": "Matching enriched terms"},
            ],
        },
        "layout": "full",
        "maxRows": 40,
    },
]

blocks = [
    {"id": "title", "type": "markdown", "body": "# NSCLC T-cell perturbation statistics and biological interpretation"},
    {
        "id": "technical-summary",
        "type": "markdown",
        "body": "## Technical summary\n\n**The computational screen and all six directional statistical analyses are complete.** Across 3,379 cells and 2,937,776 deletion observations, 3,698 gene-comparison tests met FDR < 0.05. Only 1,604 also shifted positively toward the goal centroid with at least 25 detections; these are the report's primary ranking set.\n\n**The result is useful for prioritization, not yet for causal or therapeutic claims.** Directional shifts are close to evenly divided between toward-goal and away-from-goal effects. Translation/ribosome, immune/inflammatory, antigen-processing, and—especially for LUAD → LUSC—oxidative-phosphorylation programs are enriched. Several leading genes are canonical lung epithelial, alveolar, stromal, or myeloid-associated transcripts despite the T-cell selection, making ambient RNA, doublets, or tissue-mixture effects a major alternative explanation.",
    },
    {"id": "metrics", "type": "metric-strip", "cardIds": ["cells-card", "deletions-card", "comparisons-card", "qualified-card"]},
    {
        "id": "status",
        "type": "markdown",
        "sourceId": "comparison-summary",
        "body": "## Completion and statistical status\n\nAll result tables are non-empty and contain one row per tested gene with no missing values or duplicate gene identifiers. Tested genes range from 11,242 in LUSC-source comparisons to 14,923 in normal-source comparisons. The coverage-qualified toward-goal rate ranges from **1.25% (LUSC → LUAD)** to **2.94% (LUAD → NORMAL)** of tested genes.",
    },
    {
        "id": "direction-interpretation",
        "type": "markdown",
        "sourceId": "comparison-summary",
        "body": "### Statistical significance is not the same as movement toward the goal\n\nThe original significance flag uses only FDR. Of 3,698 FDR-significant tests, 1,780 have a positive shift and 1,918 have a negative shift. The plot below keeps both directions visible; downstream biological ranking additionally requires a positive shift and at least 25 detections.",
    },
    {"id": "direction-chart-block", "type": "chart", "chartId": "direction-chart"},
]

interpretation = {
    "lusc_to_luad": "**LUSC → LUAD:** MMP12, KRT17, S100A2, S100A7, S100A8/S100A9, and CXCL8 lead the positive shifts. The coherent squamous/inflammatory pattern suggests that deleting source-associated inflammatory or epithelial-like features moves embeddings away from the LUSC centroid; it does not show that these genes are LUAD-promoting.",
    "lusc_to_normal": "**LUSC → NORMAL:** the leading effects are dominated by ribosomal, translational, and housekeeping genes (including RPS27, RPL21, UQCR11, and PTMA), with relatively small absolute shifts. This pattern is more compatible with a broad state/fitness axis than a specific normal-T-cell program.",
    "luad_to_lusc": "**LUAD → LUSC:** HBB, SFTPB, SFTPC, NAPSA, SCGB3A1, and stress/trafficking genes rank highly. Loss of alveolar/secretory transcripts can separate cells from a LUAD-associated centroid, but their appearance in a T-cell analysis is also a strong contamination or doublet warning.",
    "luad_to_normal": "**LUAD → NORMAL:** PIGR, MUC1, AGR2, WFDC2, SLC34A2, and KRT8/18/19 form a conspicuous epithelial/secretory program. These shifts may reflect removal of tumor-tissue RNA signatures rather than restoration of an intrinsic normal T-cell state.",
    "normal_to_luad": "**NORMAL → LUAD:** SFTPC, SFTPA1, FBLN1, DCN, MGP, IGFBP6, and ACKR1 combine alveolar, stromal, and vascular-associated signals. The mixture argues for donor/tissue-composition sensitivity and warrants cell-level decontamination review before target nomination.",
    "normal_to_lusc": "**NORMAL → LUSC:** SFTPC, SFTPA1/2, FBLN1, MGP, HSPA1B, and RAB4B again emphasize lung-environment and stress features. The direction is statistically reproducible at the cell-observation level, but cellular provenance is uncertain.",
}

for comparison, text in interpretation.items():
    label = summary.loc[summary["comparison"] == comparison, "comparison_label"].iloc[0]
    chart_id = f"goal-{comparison}"
    charts.append(
        {
            "id": chart_id,
            "title": f"{label}: top coverage-qualified goal shifts",
            "subtitle": "Mean cosine shift after gene deletion; top 10 of the qualified ranking",
            "intent": "comparison",
            "question": f"Which deletions most strongly move the source embedding toward the {label.split('→')[1].strip()} centroid?",
            "rationale": "Horizontal ranking preserves gene labels and shows the magnitude concentration within a comparison.",
            "type": "horizontalBar",
            "dataset": f"genes_{comparison}",
            "sourceId": "top-genes",
            "encodings": {
                "x": {"field": "Gene_name", "type": "nominal", "label": "Deleted gene"},
                "y": {"field": "Shift_to_goal_end", "type": "quantitative", "label": "Mean cosine shift toward goal"},
                "tooltip": [
                    {"field": "Goal_end_FDR", "type": "quantitative", "label": "Goal-end FDR"},
                    {"field": "N_Detections", "type": "quantitative", "label": "Detections"},
                    {"field": "Ensembl_ID", "type": "text", "label": "Ensembl ID"},
                ],
            },
            "layout": "full",
            "maxRows": 10,
        }
    )
    blocks.extend(
        [
            {"id": f"interpret-{comparison}", "type": "markdown", "body": text, "sourceId": "top-genes"},
            {"id": f"chart-{comparison}", "type": "chart", "chartId": chart_id},
        ]
    )

blocks.extend(
    [
        {
            "id": "pathway-findings",
            "type": "markdown",
            "body": "## Pathway analysis: recurring programs, with substantial redundancy\n\nEnrichment used each comparison's qualified gene set against its own tested-gene background. Translation/ribosome terms dominate four directions and remain present in the other two. Immune/cytokine programs occur in every direction; oxidative phosphorylation is especially prominent for LUAD → LUSC and appears more modestly for NORMAL → LUSC. Antigen-processing and inflammatory themes recur but often share genes with broader immune terms. **These are overlapping annotations, not independent pathway discoveries**, and enrichment significance does not measure effect size or causal importance.",
            "sourceId": "pathways",
        },
        {"id": "pathway-chart-block", "type": "chart", "chartId": "pathway-chart"},
        {
            "id": "scope-methods",
            "type": "markdown",
            "body": "## Scope, definitions, and methodology\n\n**Unit and comparison.** The screen deletes one expressed gene at a time in held-out cells and measures the change in cosine position relative to source and goal disease centroids. A positive `Shift_to_goal_end` means the perturbed embedding moved toward the specified goal centroid. It does **not** mean the gene is more highly expressed in the goal state, that deletion converts cell identity, or that the gene is a beneficial intervention target.\n\n**Statistics.** Geneformer goal-state statistics compare each gene's shift distribution with a pooled random all-gene perturbation distribution using a Wilcoxon rank-sum procedure, followed by Benjamini–Hochberg correction within each result table. This report defines a conservative ranking set as FDR < 0.05, positive goal shift, and at least 25 detections. FDR is not jointly controlled across all six tables.\n\n**Pathways.** g:Profiler g:GOSt tests Reactome, KEGG, and GO Biological Process terms using all genes tested in the same source-state table as the custom background and FDR correction. The displayed themes are a transparent regex consolidation of redundant significant terms; the underlying 2,762 enriched term rows remain in `source_tables/pathway_enrichment.csv`. See the [Geneformer framework](https://www.nature.com/articles/s41586-023-06139-9), [g:Profiler methodology](https://pmc.ncbi.nlm.nih.gov/articles/PMC10320099/), and [GO enrichment guidance](https://www.geneontology.org/docs/go-enrichment-analysis/).",
        },
        {
            "id": "limitations",
            "type": "markdown",
            "body": "## Limitations and robustness assessment\n\n**Overall confidence: share with caveats; hypothesis-generating.**\n\n- **Cell observations are not independent donor replicates.** The held-out cohort includes 19 LUAD, 12 normal, and only 5 LUSC donors. Donor-stratified direction consistency and donor-blocked inference are not yet available.\n- **Ambient RNA/doublets are a plausible dominant confounder.** Lung epithelial/alveolar genes (for example SFTPC, SFTPB, NAPSA, MUC1, PIGR, and keratins) and stromal/vascular genes appear among top T-cell shifts. Ambient RNA is known to create cell-type-inappropriate expression and false pathway signals; review [DecontX](https://pmc.ncbi.nlm.nih.gov/articles/PMC7059395/) and [SoupX](https://pmc.ncbi.nlm.nih.gov/articles/PMC7763177/).\n- **The null is pooled and cell-weighted.** Unequal expression frequency, donor composition, and cell-state abundance can influence detection counts and shift distributions.\n- **Centroid movement is representation-dependent.** A small cosine shift can be statistically precise without being biologically large, and the unusually large LUSC → LUAD MMP12 shift should be checked for outliers and donor concentration.\n- **Pathway terms are correlated.** Thousands of overlapping GO/Reactome/KEGG terms amplify apparent breadth; theme consolidation is descriptive and does not add a second inferential test.",
        },
        {
            "id": "next-steps",
            "type": "markdown",
            "body": "## Recommended next steps\n\n1. Recompute top-gene direction and effect size per donor, then require consistent sign in a prespecified majority of donors; report leave-one-donor-out stability, especially for the five LUSC donors.\n2. Quantify epithelial/alveolar marker burden per T cell, doublet score, and ambient-RNA correction sensitivity. Repeat rankings after excluding high-burden cells.\n3. Compare each disease-direction ranking with source-versus-goal differential expression. Distinguish deletion of a source marker from genuine movement toward a goal-specific T-cell program.\n4. Replace raw term lists with redundancy-aware pathway modules and confirm modules by donor-stratified gene-set scores.\n5. Nominate targets only after effect size, donor stability, T-cell expression specificity, and an orthogonal perturbation or validation dataset agree.",
        },
        {
            "id": "questions",
            "type": "markdown",
            "body": "## Further questions\n\n- Do the leading shifts remain after ambient-RNA correction and doublet exclusion?\n- Are effects distributed across donors, or driven by one tissue specimen or sequencing batch?\n- Which shifts are reciprocal (source → goal and goal → source) versus disease-direction specific?\n- Do qualified deletions move canonical T-cell functional programs, or primarily remove tissue-of-origin transcripts?\n- Which candidates reproduce in a donor-balanced external NSCLC T-cell cohort?",
        },
    ]
)

tables = [
    {
        "id": "comparison-table",
        "title": "Exact per-comparison statistical counts",
        "dataset": "comparison_summary",
        "sourceId": "comparison-summary",
        "density": "dense",
        "layout": "full",
        "columns": [
            {"field": "comparison_label", "label": "Comparison", "type": "text"},
            {"field": "genes_tested", "label": "Genes tested", "format": "number"},
            {"field": "fdr_significant", "label": "FDR significant", "format": "number"},
            {"field": "significant_toward_goal", "label": "Toward goal", "format": "number"},
            {"field": "significant_away_from_goal", "label": "Away from goal", "format": "number"},
            {"field": "qualified_toward_goal", "label": "Qualified", "format": "number"},
            {"field": "qualified_fraction", "label": "Qualified rate", "format": "percent"},
        ],
    }
]
blocks.insert(6, {"id": "comparison-table-block", "type": "table", "tableId": "comparison-table"})

datasets = {
    "overview": records(overview),
    "comparison_summary": records(summary),
    "pathway_themes": records(themes),
}
for comparison in summary["comparison"]:
    datasets[f"genes_{comparison}"] = records(
        genes[genes["comparison"] == comparison].sort_values("rank").head(10)
    )

artifact = {
    "surface": "report",
    "manifest": {
        "version": 1,
        "surface": "report",
        "title": "NSCLC T-cell perturbation statistics and biological interpretation",
        "description": "Completed directional Geneformer goal-state statistics, goal-shift rankings, pathway enrichment, and robustness caveats.",
        "generatedAt": GENERATED_AT,
        "cards": cards,
        "charts": charts,
        "tables": tables,
        "sources": sources,
        "blocks": blocks,
    },
    "snapshot": {"version": 1, "generatedAt": GENERATED_AT, "status": "ready", "datasets": datasets},
    "sources": sources,
}

(HERE / "artifact.json").write_text(json.dumps(artifact, indent=2, ensure_ascii=False) + "\n")
print(f"Wrote {HERE / 'artifact.json'} with {len(blocks)} blocks and {len(charts)} charts")
