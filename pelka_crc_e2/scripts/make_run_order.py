"""E2: the ISP run order, frozen before any GPU call (registration s.5.4), plus its cost prediction.

Order: the 100 null genes in draw order (H2b), then the eligible Panel B genes whose stratum has 20
controls, sorted by Ensembl ID (H2c), then every control gene once, sorted by Ensembl ID. Panel genes in a
control-less stratum are not run (they would be NOT_ESTIMABLE_CONTROLS rows). The no-op spot-check genes are
every 20th gene of the null list and every 20th gene of the panel+control list.
Cost: the balanced-donor Phase 6 per-call fit (cost_fit_thinkstation1.json), applied to each (gene, donor)
token-positive count, both operations. CPU only.
"""
import argparse
import json
import pickle


def order(null, strata, elig_rows):
    ok_strata = {k for k, v in strata.items() if v["status"] == "eligible" and len(v.get("controls") or []) == 20}
    panel = sorted(r["ensembl_id"] for r in elig_rows if r["eligible"] == "True" and r["stratum"] in ok_strata)
    controls = sorted({c for k in ok_strata for c in strata[k]["controls"]})
    rest = panel + [c for c in controls if c not in set(panel)]
    if set(null) & set(rest):
        raise SystemExit("null genes overlap panel/control genes")
    spot = [null[i] for i in range(len(null)) if (i + 1) % 20 == 0] + [rest[i] for i in range(len(rest)) if (i + 1) % 20 == 0]
    return null + rest, panel, controls, spot


def main():
    import csv
    ap = argparse.ArgumentParser()
    ap.add_argument("--null-genes", required=True)
    ap.add_argument("--strata", required=True)
    ap.add_argument("--eligibility", required=True)
    ap.add_argument("--dataset", required=True)
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--cost-fit", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    from datasets import load_from_disk
    null = json.load(open(a.null_genes))["draw_order"]
    strata = json.load(open(a.strata))
    elig = list(csv.DictReader(open(a.eligibility)))
    genes, panel, controls, spot = order(null, strata, elig)
    tok = pickle.load(open(a.token_dict, "rb"))
    fit = json.load(open(a.cost_fit))["fit"]
    ds = load_from_disk(a.dataset).filter(lambda x: int(x["in_analysis100"]) == 1 and x["origin"] == "tumor_primary")
    sets = {}
    for d, ids in zip(ds["individual"], ds["input_ids"]):
        sets.setdefault(d, []).append(set(ids))
    sec = 0.0
    for g in genes:
        t = tok[g]
        for d in sets:
            k = sum(t in s for s in sets[d])
            if k:
                sec += fit["delete"]["intercept_a"] + fit["delete"]["slope_b_per_cell"] * k
                sec += fit["overexpress"]["intercept_a"] + fit["overexpress"]["slope_b_per_cell"] * k
    out = {"about": __doc__.strip(), "n_genes": len(genes), "n_null": len(null), "n_panel": len(panel),
           "n_controls": len(controls), "n_donors": len(sets), "genes": genes, "panel": panel,
           "noop_spotcheck_genes": spot, "predicted_isp_gpu_h": sec / 3600}
    json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k not in ("genes", "panel", "about")}, indent=1))


if __name__ == "__main__":
    main()
