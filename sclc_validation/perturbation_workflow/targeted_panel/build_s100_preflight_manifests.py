#!/usr/bin/env python3
"""CPU-only, gate-3(a) manifest build for s100-isp-execution-20260922.

Builds the per-gene, per-donor eligibility manifest for the 12-gene S100
core panel (sclc_validation/perturbation_workflow/s100_isp/s100_gene_panel_20260922.json)
against the two Module-A-eligible sources (sclc, luad) -- normal is skipped
deliberately: it has a single donor total in the held-out split, so
MIN_DONORS_ELIGIBLE=3 can never be met there for any gene (Pam's own
preflight already excludes normal from Module A), and this is CPU-only so
there's no cost saved by including it, only noise in the summary.

Uses run_targeted_panel.py's own paired_eligible_dataset() unmodified --
no GPU, no model load, same function the real run will use, so eligibility
here is exactly what the real run will see.

Run from sclc_validation/perturbation_workflow/targeted_panel/:
    TARGET_GENES_FILE_OVERRIDE=../s100_isp/s100_gene_panel_20260922.json \
        python3 build_s100_preflight_manifests.py
"""
from __future__ import annotations

import json
import sys

sys.argv = ["run_targeted_panel.py", "--run-tag", "isp_preflight_manifests"]
sys.path.insert(0, ".")
import run_targeted_panel as rtp  # noqa: E402

PANEL_FILE = "../s100_isp/s100_gene_panel_20260922.json"
SOURCES = ("sclc", "luad")


def main() -> None:
    panel = json.load(open(PANEL_FILE))
    genes = panel["genes"]
    rows = []
    for gene in genes:
        symbol, ensembl = gene["gene"], gene["ensembl_id"]
        for source in SOURCES:
            _, info = rtp.paired_eligible_dataset(source, {"gene": symbol, "ensembl_id": ensembl})
            rows.append({
                "gene": symbol,
                "ensembl_id": ensembl,
                "stratum": gene["stratum"],
                "role": gene["role"],
                "source": source,
                "eligible": info["eligible"],
                "n_token_positive_cells": info.get("n_cells"),
                "n_donors": info.get("n_donors"),
                "donor_cell_counts": info.get("donor_cell_counts"),
                "reason": info.get("reason"),
                "manifest": info.get("manifest"),
            })

    df = rtp.pd.DataFrame(rows)
    out_path = rtp.TABLE_ROOT / "s100_preflight_eligibility.csv"
    rtp.TABLE_ROOT.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    n_eligible = int(df["eligible"].sum())
    print(f"Wrote {out_path}: {len(df)} rows, {n_eligible} eligible / {len(df)} total")
    print()
    print("=== Per-row eligibility (gene / source / eligible / cells / donors) ===")
    for _, r in df.iterrows():
        print(f"{r['gene']:>10s} / {r['source']:>4s}: eligible={r['eligible']!s:5s} "
              f"cells={r['n_token_positive_cells']} donors={r['n_donors']} "
              f"reason={r['reason']}")

    print()
    ineligible = df[~df["eligible"]]
    if len(ineligible):
        print(f"=== {len(ineligible)} NOT-ESTIMABLE (gene, source) pairs ===")
        for _, r in ineligible.iterrows():
            print(f"  {r['gene']}/{r['source']}: {r['reason']}")
    else:
        print("All (gene, source) pairs eligible.")

    a2 = df[(df["gene"] == "S100A2")]
    print()
    print("=== S100A2 detail (the decisive non-circular high-risk test) ===")
    for _, r in a2.iterrows():
        print(f"  {r['source']}: eligible={r['eligible']} cells={r['n_token_positive_cells']} "
              f"donors={r['n_donors']} donor_cell_counts={r['donor_cell_counts']} reason={r['reason']}")


if __name__ == "__main__":
    main()
