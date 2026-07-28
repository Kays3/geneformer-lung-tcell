#!/usr/bin/env python3
"""Build the 50-gene targeted perturbation panel: the audit's pre-registered
21-gene signature panel, plus the top 29 driver genes (by |shift|, unique,
excluding contamination markers) from the prior LUAD/LUSC/normal all-gene
perturbation screen.
"""
import csv
import json
from pathlib import Path

PANEL_GENES = {
    "PDCD1": "ENSG00000188389", "CTLA4": "ENSG00000163599", "HAVCR2": "ENSG00000135077",
    "LAG3": "ENSG00000089692", "TIGIT": "ENSG00000181847", "TOX": "ENSG00000198846",
    "LAYN": "ENSG00000204381", "NKG7": "ENSG00000105374", "GNLY": "ENSG00000115523",
    "PRF1": "ENSG00000180644", "GZMB": "ENSG00000100453", "GZMH": "ENSG00000100450",
    "IFNG": "ENSG00000111537", "TCF7": "ENSG00000081059", "SLAMF6": "ENSG00000162739",
    "IL7R": "ENSG00000168685", "CCR7": "ENSG00000126353", "ASCL1": "ENSG00000139352",
    "NEUROD1": "ENSG00000162992", "POU2F3": "ENSG00000137709", "YAP1": "ENSG00000137693",
}

CONTAMINATION_MARKERS = {"SFTPC", "SFTPB", "NAPSA", "MUC1", "PIGR", "FBLN1", "DCN", "ACKR1"}

TOP_GENES_CSV = Path("top_goal_shift_genes.csv")
N_TOP_DRIVERS = 29


def main():
    rows = list(csv.DictReader(TOP_GENES_CSV.open()))
    contamination = set(CONTAMINATION_MARKERS)
    best = {}
    for r in rows:
        g = r["Gene_name"]
        if g.startswith("KRT"):
            contamination.add(g)
        shift = abs(float(r["Shift_to_goal_end"]))
        if g not in best or shift > best[g][0]:
            best[g] = (shift, r)

    filtered = {g: v for g, v in best.items() if g not in contamination and g not in PANEL_GENES}
    ranked = sorted(filtered.values(), key=lambda x: -x[0])[:N_TOP_DRIVERS]

    target = []
    for symbol, ensembl_id in PANEL_GENES.items():
        target.append({"gene": symbol, "ensembl_id": ensembl_id, "source": "panel"})
    for shift, r in ranked:
        target.append({
            "gene": r["Gene_name"],
            "ensembl_id": r["Ensembl_ID"],
            "source": "top_driver_luad_lusc_normal",
            "origin_comparison": r["comparison_label"],
            "origin_shift": shift,
        })

    assert len(target) == len(PANEL_GENES) + N_TOP_DRIVERS
    assert len({t["gene"] for t in target}) == len(target), "duplicate genes in target list"

    out = {
        "n_genes": len(target),
        "n_panel": len(PANEL_GENES),
        "n_top_drivers": N_TOP_DRIVERS,
        "contamination_markers_excluded": sorted(contamination),
        "genes": target,
    }
    Path("target_gene_panel.json").write_text(json.dumps(out, indent=2) + "\n")
    print(f"Wrote target_gene_panel.json: {len(target)} genes ({len(PANEL_GENES)} panel + {N_TOP_DRIVERS} top drivers)")
    for t in target:
        print(" ", t["gene"], t["ensembl_id"], t["source"])


if __name__ == "__main__":
    main()
