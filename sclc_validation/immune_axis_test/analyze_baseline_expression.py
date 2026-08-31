#!/usr/bin/env python3
"""Analyze T2's baseline-expression measurements against PLAN.md's own
pre-registered criteria (section 8), and cross-validate the HK-gene review
report's detection-fraction proxy against this ground-truth measurement.

Reads the two tables `measure_baseline_expression.py` produced on the compute
host and pulled back (`results/baseline_expression_{pooled,per_donor}.csv`),
plus the whole-genome HK diagnostic's own detection-fraction-from-N_Detections
proxy, and writes a per-gene/per-program summary table plus two figures. No
new measurement; this only summarizes what T2 already wrote.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from adjustText import adjust_text

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
HK_GAP_TABLE = ROOT / "sclc_validation/primary_test_perturbation/tables/hk_gene_diagnostic/allgene_hk_concordant_detection_gap.csv"

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

PROGRAM_ORDER = ["exhaustion", "cytotoxicity", "progenitor", "sclc_subtype_tf"]
STATE_ORDER = ["normal", "sclc", "luad"]
STATE_LABEL = {"normal": "Normal", "sclc": "SCLC", "luad": "LUAD"}


def program_summary(pooled: pd.DataFrame) -> pd.DataFrame:
    programs = pooled[pooled["program"].isin(PROGRAM_ORDER)]
    summary = programs.groupby(["program", "state"])[["detect_rate", "mean_log1p_cp10k"]].mean().reset_index()
    return summary


def check_t2_criterion(summary: pd.DataFrame) -> str:
    exh = summary[summary["program"] == "exhaustion"].set_index("state")
    sclc_vs_luad_detect = exh.loc["sclc", "detect_rate"] - exh.loc["luad", "detect_rate"]
    sclc_vs_luad_expr = exh.loc["sclc", "mean_log1p_cp10k"] - exh.loc["luad", "mean_log1p_cp10k"]
    overturned = sclc_vs_luad_detect >= 0 or sclc_vs_luad_expr >= 0
    lines = [
        "PLAN.md section 8 falsification check: "
        "'T2 shows SCLC T cells at equal or higher exhaustion expression than LUAD T cells'",
        f"  SCLC - LUAD, exhaustion detect_rate:     {sclc_vs_luad_detect:+.4f}",
        f"  SCLC - LUAD, exhaustion mean_log1p_cp10k: {sclc_vs_luad_expr:+.4f}",
        f"  => criterion {'TRIGGERED (contradiction sharpens)' if overturned else 'NOT triggered (SCLC measurably below LUAD on both metrics)'}",
    ]
    return "\n".join(lines)


def plot_program_by_state(summary: pd.DataFrame, per_donor: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 4, figsize=(15, 4.2), dpi=DPI, sharey=False)
    for ax, program in zip(axes, PROGRAM_ORDER):
        sub = per_donor[per_donor["program"] == program]
        genes = sorted(sub["gene"].unique())
        x = np.arange(len(STATE_ORDER))
        for gene in genes:
            gsub = sub[sub["gene"] == gene]
            means = [gsub[gsub["state"] == s]["detect_rate"].mean() for s in STATE_ORDER]
            ax.plot(x, means, marker="o", markersize=4, linewidth=1, alpha=0.75, label=gene)
        # Normal is a single donor -- mark its point distinctly on every line.
        ax.axvspan(-0.4, 0.4, color=GREY, alpha=0.15, zorder=0)
        ax.set_xticks(x)
        ax.set_xticklabels([STATE_LABEL[s] for s in STATE_ORDER], fontsize=9, color=INK2)
        ax.set_title(program, fontsize=10, color=INK, loc="left", fontweight="bold")
        ax.tick_params(labelsize=8, colors=INK2)
        for spine in ("top", "right"):
            ax.spines[spine].set_visible(False)
        for spine in ("left", "bottom"):
            ax.spines[spine].set_color(INK2)
        ax.legend(fontsize=6.5, frameon=False, loc="upper left", bbox_to_anchor=(0.0, -0.08), ncol=2)
    axes[0].set_ylabel("detection rate (held-out test cells)", fontsize=9, color=INK2)
    fig.suptitle("T cells — T2: measured detection rate by program and state\n"
                 "(shaded band = Normal, single donor)", fontsize=10.5, color=INK, y=1.08, x=0.02, ha="left", fontweight="bold")
    fig.tight_layout()
    out = FIGURES / "t2_program_detection_rate.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(ROOT)}")


def plot_proxy_validation(pooled: pd.DataFrame) -> None:
    if not HK_GAP_TABLE.exists():
        print(f"[skip] {HK_GAP_TABLE} not found")
        return
    proxy = pd.read_csv(HK_GAP_TABLE).drop_duplicates(subset=["Gene_name", "comparison"])
    hk = pooled[pooled["program"] == "hk_concordant_whole_genome"].set_index(["gene", "state"])["detect_rate"]

    rows = []
    for _, row in proxy.iterrows():
        try:
            measured_source = hk.loc[(row["Gene_name"], row["source_state"])]
            measured_goal = hk.loc[(row["Gene_name"], row["goal_state"])]
        except KeyError:
            continue
        rows.append({
            "gene": row["Gene_name"], "comparison": row["comparison"],
            "proxy_source": row["source_detect_frac"], "measured_source": measured_source,
            "proxy_goal": row["goal_detect_frac"], "measured_goal": measured_goal,
        })
    val = pd.DataFrame(rows)
    val_long = pd.concat([
        val[["gene", "proxy_source", "measured_source"]].rename(columns={"proxy_source": "proxy", "measured_source": "measured"}),
        val[["gene", "proxy_goal", "measured_goal"]].rename(columns={"proxy_goal": "proxy", "measured_goal": "measured"}),
    ], ignore_index=True).drop_duplicates()

    r = val_long["proxy"].corr(val_long["measured"])
    fig, ax = plt.subplots(figsize=(5.5, 5.5), dpi=DPI)
    ax.plot([0, 1], [0, 1], color=INK2, linewidth=0.7, linestyle=(0, (3, 2)), zorder=1)
    ax.scatter(val_long["proxy"], val_long["measured"], s=18, color=BLUE, alpha=0.6, edgecolors="none", zorder=2)
    ax.set_xlabel("proxy detection fraction (from ISP N_Detections)", fontsize=9, color=INK2)
    ax.set_ylabel("measured detection rate (T2, held-out atlas cells)", fontsize=9, color=INK2)
    ax.set_title(f"T cells — HK-gene proxy validation (r = {r:.3f}, n = {len(val_long)})",
                 fontsize=10, color=INK, loc="left", pad=8, fontweight="bold")
    ax.tick_params(labelsize=8, colors=INK2)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(INK2)
    fig.tight_layout()
    out = FIGURES / "t2_hk_proxy_validation.png"
    fig.savefig(out, dpi=DPI, bbox_inches="tight")
    print(f"Wrote {out.relative_to(ROOT)}  (Pearson r = {r:.3f})")


def main() -> None:
    pooled = pd.read_csv(RESULTS / "baseline_expression_pooled.csv")
    per_donor = pd.read_csv(RESULTS / "baseline_expression_per_donor.csv")
    FIGURES.mkdir(parents=True, exist_ok=True)

    summary = program_summary(pooled)
    summary.to_csv(RESULTS / "t2_program_summary.csv", index=False)
    print(summary.pivot(index="program", columns="state", values="detect_rate")[STATE_ORDER].round(4).to_string())
    print()
    print(check_t2_criterion(summary))
    print()

    plot_program_by_state(summary, per_donor)
    plot_proxy_validation(pooled)


if __name__ == "__main__":
    main()
