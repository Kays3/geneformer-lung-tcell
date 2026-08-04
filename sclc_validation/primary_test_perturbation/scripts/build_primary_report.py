#!/usr/bin/env python3
"""Build the locked, donor-held-out primary perturbation report."""

from __future__ import annotations

import html
import json
import os
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parents[1]
ROOT = Path(os.environ.get("SCLC_PERTURBATION_ROOT", Path.home() / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation"))
STATS = ROOT / "stats"
TABLES = HERE / "tables"
REPORTS = HERE / "reports"
COMPARISONS = ("sclc_to_luad", "sclc_to_normal", "luad_to_sclc", "luad_to_normal", "normal_to_sclc", "normal_to_luad")
ARMS = ("delete", "overexpress")
FDR = 0.05
MIN_DETECTIONS = 25


def path_for(arm: str, comparison: str) -> Path:
    return STATS / arm / f"heldout_allgene_{arm}_{comparison}.csv"


def read_table(arm: str, comparison: str) -> tuple[pd.DataFrame | None, dict]:
    path = path_for(arm, comparison)
    audit = {"arm": arm, "comparison": comparison, "path": str(path), "exists": path.exists(), "status": "missing", "rows": 0}
    if not path.exists():
        return None, audit
    frame = pd.read_csv(path).drop(columns=["Unnamed: 0"], errors="ignore")
    required = {"Gene_name", "Ensembl_ID", "Goal_end_FDR", "Shift_to_goal_end", "N_Detections"}
    missing = sorted(required.difference(frame.columns))
    audit.update({"rows": len(frame), "missing_columns": ";".join(missing), "status": "ok" if not missing else "invalid"})
    if missing:
        return None, audit
    frame["Goal_end_FDR"] = pd.to_numeric(frame["Goal_end_FDR"], errors="coerce")
    frame["Shift_to_goal_end"] = pd.to_numeric(frame["Shift_to_goal_end"], errors="coerce")
    frame["N_Detections"] = pd.to_numeric(frame["N_Detections"], errors="coerce")
    frame["comparison"] = comparison
    frame["arm"] = arm
    frame["qualified_primary"] = (frame["Goal_end_FDR"] < FDR) & (frame["Shift_to_goal_end"] > 0) & (frame["N_Detections"] >= MIN_DETECTIONS)
    audit["duplicate_ensembl_ids"] = int(frame["Ensembl_ID"].duplicated().sum())
    audit["n_qualified"] = int(frame["qualified_primary"].sum())
    audit["max_detections"] = int(frame["N_Detections"].max()) if len(frame) else 0
    return frame, audit


def build() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    TABLES.mkdir(parents=True, exist_ok=True)
    audits, summaries, frames = [], [], {}
    for arm in ARMS:
        for comparison in COMPARISONS:
            frame, audit = read_table(arm, comparison)
            audits.append(audit)
            if frame is not None:
                frames[(arm, comparison)] = frame
            qualified = frame[frame["qualified_primary"]] if frame is not None else pd.DataFrame()
            summaries.append({"arm": arm, "comparison": comparison, "status": audit["status"], "genes_tested": len(frame) if frame is not None else 0, "qualified_primary": len(qualified), "fdr_significant": int((frame["Goal_end_FDR"] < FDR).sum()) if frame is not None else 0, "positive_shift": int((frame["Shift_to_goal_end"] > 0).sum()) if frame is not None else 0})

    hits = []
    for comparison in COMPARISONS:
        delete = frames.get(("delete", comparison))
        over = frames.get(("overexpress", comparison))
        if delete is None or over is None:
            continue
        d = delete[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(columns={"Shift_to_goal_end": "delete_shift", "Goal_end_FDR": "delete_fdr", "N_Detections": "delete_n"})
        o = over[["Gene_name", "Ensembl_ID", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"]].rename(columns={"Shift_to_goal_end": "overexpress_shift", "Goal_end_FDR": "overexpress_fdr", "N_Detections": "overexpress_n"})
        merged = d.merge(o, on=["Gene_name", "Ensembl_ID"], how="outer")
        merged["comparison"] = comparison
        merged["concordant_primary"] = (merged["delete_fdr"] < FDR) & (merged["overexpress_fdr"] < FDR) & (merged["delete_shift"] * merged["overexpress_shift"] < 0) & (merged["delete_n"] >= MIN_DETECTIONS)
        hits.append(merged[merged["concordant_primary"]])
    hit_table = pd.concat(hits, ignore_index=True) if hits else pd.DataFrame(columns=["comparison", "Gene_name", "Ensembl_ID", "delete_shift", "overexpress_shift", "delete_fdr", "overexpress_fdr", "delete_n", "overexpress_n", "concordant_primary"])
    if not hit_table.empty:
        hit_table["abs_delete_shift"] = hit_table["delete_shift"].abs()
        hit_table = hit_table.sort_values(["comparison", "abs_delete_shift"], ascending=[True, False])
    summary = pd.DataFrame(summaries)
    audit = pd.DataFrame(audits)
    summary.to_csv(TABLES / "primary_arm_summary.csv", index=False)
    audit.to_csv(TABLES / "coverage_audit.csv", index=False)
    hit_table.to_csv(TABLES / "primary_concordant_hits.csv", index=False)
    return summary, audit, hit_table


def report(summary: pd.DataFrame, audit: pd.DataFrame, hits: pd.DataFrame) -> None:
    REPORTS.mkdir(parents=True, exist_ok=True)
    complete = bool(((audit["status"] == "ok") & (audit["rows"].fillna(0) > 0)).all())
    status = "COMPLETE" if complete else "PENDING: computation artifacts are incomplete"
    lines = ["# Primary donor-held-out SCLC perturbation report", "", f"**Status:** {status}", f"**Generated:** {datetime.now(timezone.utc).isoformat()}", "", "## Decision boundary", "", "Only test-cell perturbations with training-only reference centroids are primary evidence. Training, evaluation, and pooled whole-cohort perturbations are excluded from this report's primary estimates.", "", "## Completeness", "", f"- Required arm/comparison tables: {int((audit['status'] == 'ok').sum())}/{len(audit)}", f"- Primary concordant gene-comparison hits: {len(hits)}", f"- FDR threshold: {FDR}", f"- Minimum detections: {MIN_DETECTIONS}", "", "## Arm summary", "", "```text", summary.to_string(index=False), "```", "", "## Coverage audit", "", "```text", audit.to_string(index=False), "```", "", "## Interpretation", "", "A concordant hit requires significant, adequately detected, opposite-signed deletion and overexpression shifts. This is a model-level prioritization criterion, not evidence of experimental causality. Donor-level consistency, ambient-RNA/doublet sensitivity, and independent validation remain required before biological claims.", "", "" if complete else "The report is intentionally incomplete. Do not interpret missing comparisons as null effects."]
    markdown = "\n".join(lines)
    (REPORTS / "primary_test_perturbation_report.md").write_text(markdown)
    table_html = summary.to_html(index=False, classes="data") + audit.to_html(index=False, classes="data")
    document = f"<html><head><meta charset='utf-8'><title>Primary SCLC perturbation</title><style>body{{font-family:system-ui;max-width:1200px;margin:2rem auto;line-height:1.45}} table{{border-collapse:collapse;margin:1rem 0}} th,td{{border:1px solid #ccd;padding:.35rem}} th{{background:#eef2f5}} code{{background:#f4f4f4}}</style></head><body><pre style='white-space:pre-wrap'>{html.escape(markdown.split('## Arm summary')[0])}</pre><h2>Arm summary and audit</h2>{table_html}<h2>Interpretation</h2><p>A concordant hit is a model-level prioritization criterion, not experimental causal evidence. Donor-level consistency and contamination sensitivity remain required.</p></body></html>"
    (REPORTS / "primary_test_perturbation_report.html").write_text(document)
    (HERE / "analysis_manifest.json").write_text(json.dumps({"generated_utc": datetime.now(timezone.utc).isoformat(), "artifact_root": str(ROOT), "primary_dataset": "donor-held-out test cells", "reference_dataset": "training cells only", "complete": complete, "arms": list(ARMS), "comparisons": list(COMPARISONS), "fdr_threshold": FDR, "minimum_detections": MIN_DETECTIONS}, indent=2) + "\n")


if __name__ == "__main__":
    report(*build())
