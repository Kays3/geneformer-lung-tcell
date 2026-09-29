"""Amendment 5 (P/R4): ambient flag for the drawn null genes, reusing ambient_loao.py's fitted
classifier and features UNCHANGED (non-circular, same training-pool cells, same LOAO-validated
threshold). Scores whichever ensembl_ids are passed in --genes-file's "draw_order" (or a prefix of
it, via --limit), instead of Panel A. CPU only.
"""
import argparse
import json
import os
import sys

import anndata as ad
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from ambient_loao import anchor_lists, classify, features  # noqa: E402


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--slim", required=True)
    p.add_argument("--symbols", required=True)
    p.add_argument("--diagnostic", required=True)
    p.add_argument("--genes-file", required=True, help="null_genes.json (or a prefix of it)")
    p.add_argument("--limit", type=int, default=None, help="score only the first N of draw_order")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    adata = ad.read_h5ad(a.slim)
    sym = pd.read_csv(a.symbols)
    s2e = sym.groupby("symbol").ensembl_id.agg(list)
    ka, kt = anchor_lists(a.diagnostic)
    amb = {s2e[g][0] for g in ka if g in s2e.index and len(s2e[g]) == 1}
    tc = {s2e[g][0] for g in kt if g in s2e.index and len(s2e[g]) == 1}

    draw = json.load(open(a.genes_file))["draw_order"]
    genes = draw[: a.limit] if a.limit else draw

    feat = features(adata)
    table, summary, _ = classify(feat, amb, tc, genes)
    summary["n_cells"] = int(adata.n_obs)
    summary["n_null_genes_scored_for"] = len(genes)
    json.dump({"summary": summary, "rows": table.to_dict(orient="records")}, open(a.out, "w"), indent=1)
    n_flagged = int(table.treated_as_flagged.sum())
    print(f"scored {len(genes)} null genes: {n_flagged} treated_as_flagged, {len(genes) - n_flagged} not")


if __name__ == "__main__":
    main()
