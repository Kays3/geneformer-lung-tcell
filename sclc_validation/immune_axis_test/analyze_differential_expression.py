#!/usr/bin/env python3
"""Analyze T5's genome-wide DE output (PLAN.md T5a/T5b). Run after pulling
back `differential_expression.py`'s output for both POPULATION values.

T5a (complete atlas, results/*_complete.csv):
  1. Where the 7 pre-registered exhaustion genes rank among all genome-wide
     SCLC-vs-LUAD DE genes.
  2. A data-driven dysfunction score (top-K LUAD-up genes by the data itself,
     not a curated list) vs. the curated exhaustion program, both computed on
     the same complete-atlas pooled expression table for a fair comparison.

T5b (held-out test only, results/*_test_only.csv):
  3. What fraction of whole-genome ISP hits (delete-vs-overexpress concordant,
     goal-vs-alt significant) are also DE-significant in the matching
     comparison, against the tested-background rate (Fisher's exact test).

Every number is read from tables `differential_expression.py` and the
existing whole-genome ISP tables already wrote; this script only summarizes.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
ISP_DVO = ROOT / "sclc_validation/primary_test_perturbation/tables/allgene_delete_overexpress_shift.csv"
ISP_GOAL_ALT = ROOT / "sclc_validation/primary_test_perturbation/tables/allgene_goal_vs_alt_shift.csv"

BLUE = "#17608f"
ORANGE = "#c2691a"
GREY = "#b7c2ce"
INK = "#16202c"
INK2 = "#4a5768"
DPI = 400
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
    "font.size": 9,
    "axes.linewidth": 0.75,
})

EXHAUSTION_GENES = ["PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT", "TOX", "LAYN"]
TOP_K = 50
FDR = 0.05

# Undirected DE pair -> the ISP tables' two directed comparisons that involve it.
PAIR_TO_ISP_COMPARISONS = {
    "sclc_vs_luad": ["sclc_to_luad", "luad_to_sclc"],
    "sclc_vs_normal": ["sclc_to_normal", "normal_to_sclc"],
    "luad_vs_normal": ["luad_to_normal", "normal_to_luad"],
}


def load(population: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    de = pd.read_csv(RESULTS / f"differential_expression_{population}.csv")
    pooled = pd.read_csv(RESULTS / f"state_pooled_expression_{population}.csv")
    return de, pooled


def goal1_exhaustion_rank(de_complete: pd.DataFrame) -> pd.DataFrame:
    sub = de_complete[de_complete["comparison"] == "sclc_vs_luad"].copy()
    sub = sub.sort_values(["padj", "pval"]).reset_index(drop=True)
    sub["rank"] = np.arange(1, len(sub) + 1)
    sub["percentile"] = 100 * sub["rank"] / len(sub)
    result = sub[sub["gene_symbol"].isin(EXHAUSTION_GENES)][
        ["gene_symbol", "log2fc", "pval", "padj", "rank", "percentile"]
    ].rename(columns={"gene_symbol": "gene"})
    result.to_csv(RESULTS / "t5a_exhaustion_gene_rank.csv", index=False)
    return result.sort_values("percentile")


def goal2_dysfunction_score(de_complete: pd.DataFrame, pooled_complete: pd.DataFrame) -> pd.DataFrame:
    sub = de_complete[de_complete["comparison"] == "sclc_vs_luad"]
    # group=sclc, reference=luad: negative log2fc means higher in LUAD.
    luad_up = sub[sub["log2fc"] < 0].sort_values("padj").head(TOP_K)
    data_driven_genes = luad_up["gene_symbol"].tolist()

    def program_score(genes: list[str], label: str) -> pd.DataFrame:
        rows = []
        for state in ("normal", "sclc", "luad"):
            vals = pooled_complete[(pooled_complete["gene_symbol"].isin(genes)) & (pooled_complete["state"] == state)]
            rows.append({"program": label, "state": state, "n_genes": len(genes),
                         "mean_detect_rate": vals["detect_rate"].mean(),
                         "mean_expression": vals["mean_log1p_cp10k"].mean()})
        return pd.DataFrame(rows)

    comparison = pd.concat([
        program_score(EXHAUSTION_GENES, "curated_exhaustion"),
        program_score(data_driven_genes, f"data_driven_top{TOP_K}_luad_up"),
    ], ignore_index=True)
    comparison.to_csv(RESULTS / "t5a_dysfunction_score_comparison.csv", index=False)
    pd.DataFrame({"gene": data_driven_genes}).to_csv(RESULTS / "t5a_data_driven_gene_set.csv", index=False)
    return comparison


def goal3_isp_overlap(de_test_only: pd.DataFrame) -> pd.DataFrame:
    dvo = pd.read_csv(ISP_DVO)
    goal_alt = pd.read_csv(ISP_GOAL_ALT)
    rows = []
    for pair, isp_comparisons in PAIR_TO_ISP_COMPARISONS.items():
        de_tested_genes = set(de_test_only[de_test_only["comparison"] == pair]["gene_symbol"])
        de_sig_genes_genome = set(
            de_test_only[(de_test_only["comparison"] == pair) & (de_test_only["padj"] < FDR)]["gene_symbol"]
        )

        concordant_genes = set(dvo[dvo["comparison"].isin(isp_comparisons) & dvo["concordant"]]["Gene_name"])
        # Universe for the Fisher test is genes tested in BOTH ISP (this pair)
        # and DE (this pair) -- ISP is detection-limited so its gene set is a
        # strict subset of DE's genome-wide one; every cell of the 2x2 table
        # must draw from the same universe or the test is not well-posed.
        tested_isp_genes = set(dvo[dvo["comparison"].isin(isp_comparisons)]["Gene_name"]) & de_tested_genes
        concordant_genes &= tested_isp_genes
        de_sig_genes = de_sig_genes_genome & tested_isp_genes

        for label, hit_set in [("delete_overexpress_concordant", concordant_genes)]:
            if not tested_isp_genes:
                continue
            in_de_and_hit = len(hit_set & de_sig_genes)
            in_de_not_hit = len(de_sig_genes - hit_set)
            hit_not_de = len(hit_set - de_sig_genes)
            neither = len(tested_isp_genes - hit_set - de_sig_genes)
            table = [[in_de_and_hit, hit_not_de], [in_de_not_hit, neither]]
            odds_ratio, p_value = fisher_exact(table)
            rows.append({
                "pair": pair, "isp_hit_type": label, "n_isp_hits": len(hit_set),
                "n_isp_hits_also_de_sig": in_de_and_hit,
                "pct_isp_hits_also_de_sig": round(100 * in_de_and_hit / len(hit_set), 1) if hit_set else float("nan"),
                "background_pct_de_sig_isp_universe": round(100 * len(de_sig_genes) / len(tested_isp_genes), 2),
                "background_pct_de_sig_genomewide": round(100 * len(de_sig_genes_genome) / len(de_tested_genes), 2),
                "odds_ratio": round(odds_ratio, 3), "fisher_p": p_value,
            })
    result = pd.DataFrame(rows)
    result.to_csv(RESULTS / "t5b_isp_overlap.csv", index=False)
    return result


def plot_exhaustion_rank(rank_table: pd.DataFrame, n_total: int) -> None:
    fig, ax = plt.subplots(figsize=(6, 4.5), dpi=DPI)
    ax.barh(rank_table["gene"], 100 - rank_table["percentile"], color=BLUE)
    ax.set_xlabel("percentile of genome-wide DE significance (SCLC vs LUAD)", fontsize=9, color=INK2)
    ax.set_xlim(0, 100)
    ax.axvline(95, color=ORANGE, linestyle=(0, (3, 2)), linewidth=1, label="top 5%")
    ax.set_title(f"T cells — where do exhaustion genes rank among {n_total:,} genes?",
                 fontsize=10, color=INK, loc="left", pad=8, fontweight="bold")
    ax.legend(frameon=False, fontsize=8)
    ax.tick_params(labelsize=8, colors=INK2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)
    fig.tight_layout()
    out = FIGURES / "t5a_exhaustion_gene_rank.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(ROOT)}")


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    de_complete, pooled_complete = load("complete")
    de_test_only, _ = load("test_only")

    rank_table = goal1_exhaustion_rank(de_complete)
    print("Goal 1 -- exhaustion gene rank, SCLC vs LUAD:")
    print(rank_table.to_string(index=False))
    n_total = len(de_complete[de_complete["comparison"] == "sclc_vs_luad"])
    plot_exhaustion_rank(rank_table, n_total)

    score_comparison = goal2_dysfunction_score(de_complete, pooled_complete)
    print("\nGoal 2 -- curated vs data-driven dysfunction score:")
    print(score_comparison.pivot(index="program", columns="state", values="mean_expression")[["normal", "sclc", "luad"]].round(4).to_string())

    overlap = goal3_isp_overlap(de_test_only)
    print("\nGoal 3 -- ISP hit / DE overlap:")
    print(overlap.to_string(index=False))


if __name__ == "__main__":
    main()
