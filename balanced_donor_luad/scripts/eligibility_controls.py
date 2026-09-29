"""Step 0.4 (+ host assignment): eligibility, matched-control strata, gene->host map. CPU only.

Source cells: every donor's 100 held-out tumour analysis cells (k-fold: all 43
donors are held out once) = the cells ISP will perturb.

  * Estimable (g, d): gene token present in >= 10 of donor d's 100 cells.
    Eligible: estimable in >= d_min donors (Panel A 10, Panel B 11; s.4).
    Deletion and overexpression share this criterion (Geneformer perturbs
    only token-positive cells in both operations).
  * Gene stats for matching: detect_frac and median 0-based token rank over
    the same cells, via s100_isp.matched_controls.median_token_rank_and_detection.
  * Strata: the rule proposed to god BEFORE any gene statistic was computed
    (2026-09-24, pending acceptance, not changed after): per family, sort
    ELIGIBLE panel genes by (log2 detect, rank percentile); greedily add the
    next gene to the open stratum only if within 0.5 log2 detection and 5 rank
    percentile points of ALL members, else open a new one. 20 controls per
    stratum via build_matched_control_table (seed 20260924), excluding every
    panel gene and every anchor; < 20 candidates -> NOT_ESTIMABLE_CONTROLS.
  * Host map (Amendment 3.1): within each family (Panel A, Panel B, each
    stratum's controls) genes sorted by Ensembl ID alternate ts1, ts2, ...
    A control drawn by several strata is run once, on its first assignment.
"""
import argparse
import json
import pickle
import re
import sys

import numpy as np
import pandas as pd

D_MIN = {"A": 10, "B": 11}
MIN_CELLS = 10
DET_TOL, RANK_TOL = 0.5, 5.0
SEED = 20260924


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--token-dict", required=True)
    p.add_argument("--panel-a", required=True)
    p.add_argument("--panel-b", required=True)
    p.add_argument("--symbols", required=True)
    p.add_argument("--diagnostic", required=True)
    p.add_argument("--s100", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    sys.path.insert(0, a.s100)
    import matched_controls as mc
    from datasets import load_from_disk

    tok = pickle.load(open(a.token_dict, "rb"))
    ds = load_from_disk(a.dataset).filter(lambda x: int(x["in_analysis100"]) == 1 and x["origin"] == "tumor_primary")
    donors = np.array(ds["individual"]); cells = ds["input_ids"]
    sym = pd.read_csv(a.symbols).drop_duplicates("ensembl_id").set_index("ensembl_id").symbol

    fam = {}
    for g in json.load(open(a.panel_a))["genes"]:
        fam[g["ensembl_id"]] = "A"
    for g in json.load(open(a.panel_b))["genes"]:
        if g["perturbable"]:
            fam[g["ensembl_id"]] = "B"

    # estimable donors per panel gene
    pos_sets = [set(c) for c in cells]
    rows = []
    for e, f in fam.items():
        t = tok[e]
        has = np.array([t in s for s in pos_sets])
        per_donor = pd.Series(has).groupby(donors).sum()
        est = int((per_donor >= MIN_CELLS).sum())
        rows.append({"ensembl_id": e, "symbol": sym.get(e, e), "panel": f, "estimable_donors": est,
                     "d_min": D_MIN[f], "eligible": est >= D_MIN[f]})
    elig = pd.DataFrame(rows).sort_values(["panel", "ensembl_id"])

    # gene stats over the same cells (per-source-state, as the S100 ruling requires)
    stats = mc.median_token_rank_and_detection(cells, tok)
    est_stats = stats[stats.detect_frac > 0].dropna(subset=["median_token_rank"]).set_index("ensembl_id")
    log2d = np.log2(est_stats.detect_frac)
    pct = mc.rank_percentile_transform(est_stats.median_token_rank)

    # greedy strata over eligible genes, per family
    strata = {}
    for f in ("A", "B"):
        genes = [e for e in elig[(elig.panel == f) & elig.eligible].ensembl_id if e in est_stats.index]
        genes.sort(key=lambda e: (log2d[e], pct[e]))
        cur, k = [], 0
        for e in genes:
            if cur and all(abs(log2d[e] - log2d[m]) <= DET_TOL and abs(pct[e] - pct[m]) <= RANK_TOL for m in cur):
                cur.append(e)
            else:
                if cur:
                    strata[f"{f}{k:02d}"] = cur; k += 1
                cur = [e]
        if cur:
            strata[f"{f}{k:02d}"] = cur

    # anchors excluded from every candidate pool
    s = open(a.diagnostic).read()
    anchor_syms = set(re.findall(r'"([^"]+)"', re.search(r"KNOWN_AMBIENT = \[(.*?)\]", s, re.S).group(1)))
    sym_to_e = pd.read_csv(a.symbols).groupby("symbol").ensembl_id.first()
    panel_genes = [{"gene": e, "ensembl_id": e, "role": "panel"} for e in fam]
    panel_genes += [{"gene": sym_to_e[x], "ensembl_id": sym_to_e[x], "role": "ambient_anchor"}
                    for x in anchor_syms if x in sym_to_e.index]
    table = mc.build_matched_control_table(panel_genes, stats, strata={k: tuple(v) for k, v in strata.items()},
                                           detection_tolerance_log2=DET_TOL, rank_tolerance_percentile=RANK_TOL,
                                           min_common_controls=20, seed=SEED)

    # host map
    host, used = {}, set()
    def alternate(ids):
        i = 0
        for e in sorted(ids):
            if e in used:
                continue
            host[e] = "thinkstation1" if i % 2 == 0 else "thinkstation2"; used.add(e); i += 1
    alternate(elig[(elig.panel == "A") & elig.eligible].ensembl_id)
    alternate(elig[(elig.panel == "B") & elig.eligible].ensembl_id)
    for k in sorted(table):
        if table[k]["controls"]:
            alternate(table[k]["controls"])

    stratum_of = {e: k for k, v in strata.items() for e in v}
    elig["stratum"] = elig.ensembl_id.map(stratum_of)
    elig["controls_status"] = elig.stratum.map(lambda k: table[k]["status"] if isinstance(k, str) else None)
    elig["host"] = elig.ensembl_id.map(host)
    elig.to_csv(a.out + "_eligibility.csv", index=False)
    json.dump({k: {**v, "members": strata[k], "controls": list(v["controls"]) if v["controls"] else None}
               for k, v in table.items()}, open(a.out + "_strata.json", "w"), indent=1)
    pd.Series(host, name="host").rename_axis("ensembl_id").to_csv(a.out + "_host_map.csv")
    n_controls = len({c for v in table.values() if v["controls"] for c in v["controls"]})
    summary = {
        "n_cells": len(cells), "n_donors": int(len(set(donors))),
        "panel_A_eligible": int(elig[(elig.panel == "A")].eligible.sum()), "panel_A_total": int((elig.panel == "A").sum()),
        "panel_B_eligible": int(elig[(elig.panel == "B")].eligible.sum()), "panel_B_total": int((elig.panel == "B").sum()),
        "n_strata": len(strata),
        "n_strata_not_estimable": sum(v["status"] != "eligible" for v in table.values()),
        "n_unique_controls": n_controls,
        "genes_to_run": len(host),
        "host_counts": pd.Series(host).value_counts().to_dict(),
    }
    json.dump(summary, open(a.out + "_summary.json", "w"), indent=1)
    print(json.dumps(summary, indent=1))
    print(elig.to_string(index=False))


if __name__ == "__main__":
    main()
