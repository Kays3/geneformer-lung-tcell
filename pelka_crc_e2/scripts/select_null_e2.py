"""E2 step 5: freeze the null-gene draw and its first 100 ESTIMABLE genes before any GPU output.

Combines balanced_donor_luad/scripts/select_null_genes.py (Amendment 4: the sampling frame and the seeded
shuffle) and select_estimable_prefix.py (Amendment 5: walk the draw order, keep the first N estimable genes).
Changes: a new seed (20261001), the E2 exclusions (Panel B and every E2 control / stratum member), and the
token sets come from the tokenised cohort's tumour analysis cells per donor, which are exactly the cells of
each donor's ISP input (goal_embeddings.py writes the same filter), so n_token_cells is the quantity
run_isp.py will record in its markers.

Estimable gene: token present in >= 10 of a donor's 100 cells in >= 10 donors (unchanged numbers; with 19
donors this asks for detection in more than half the cohort, where the lung study asked 10 of 43).
Outputs (both committed): null_genes_100_estimable.json ({"n": 100, "draw_order": [...]}, the format
null_analysis.py asserts) and frozen_100_estimable.json ({"rows": [...]} with npos_by_donor per walked gene,
the gate-identity source null_analysis.py checks every marker against).
"""
import argparse
import hashlib
import json
import pickle
import random

CONTROL_MIN_CELLS = 10
D_MIN = 10
SPECIAL = {"<pad>", "<mask>", "<cls>", "<eos>"}


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def frame(tok, excluded):
    return [g for g in sorted(tok) if g not in SPECIAL and g not in excluded]


def walk(draw, tok, token_sets, target):
    """token_sets: {donor: [set(token ids) per cell]}. Returns (rows, chosen)."""
    donors = sorted(token_sets)
    rows, chosen = [], []
    for i, g in enumerate(draw, start=1):
        t = tok[g]
        npos = {d: sum(t in s for s in token_sets[d]) for d in donors}
        n_est = sum(v >= CONTROL_MIN_CELLS for v in npos.values())
        est = n_est >= D_MIN
        rows.append({"gene": g, "position": i, "estimable": est, "n_estimable_donors": n_est, "npos_by_donor": npos})
        if est:
            chosen.append(g)
            if len(chosen) == target:
                break
    if len(chosen) < target:
        raise SystemExit(f"draw exhausted with {len(chosen)} estimable genes < {target}")
    return rows, chosen


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--dataset", required=True, help="tokenised pelka_crc_e2_pool300.dataset")
    ap.add_argument("--panel-b", required=True)
    ap.add_argument("--strata", required=True, help="controls/pre_gpu_strata.json")
    ap.add_argument("--seed", type=int, default=20261001)
    ap.add_argument("--target", type=int, default=100)
    ap.add_argument("--out-dir", required=True)
    a = ap.parse_args()
    from datasets import load_from_disk

    tok = pickle.load(open(a.token_dict, "rb"))
    panel_b = {g["ensembl_id"] for g in json.load(open(a.panel_b))["genes"]}
    strata = json.load(open(a.strata))
    controls = {c for v in strata.values() for c in (v.get("controls") or [])}
    members = {m for v in strata.values() for m in v["members"]}
    excluded = panel_b | controls | members
    fr = frame(tok, excluded)
    draw = fr.copy()
    random.Random(a.seed).shuffle(draw)

    ds = load_from_disk(a.dataset).filter(lambda x: int(x["in_analysis100"]) == 1 and x["origin"] == "tumor_primary")
    token_sets = {}
    for d, ids in zip(ds["individual"], ds["input_ids"]):
        token_sets.setdefault(d, []).append(set(ids))
    assert all(len(v) == 100 for v in token_sets.values()), {d: len(v) for d, v in token_sets.items()}
    rows, chosen = walk(draw, tok, token_sets, a.target)

    common = {"seed": a.seed, "token_dict_sha256": sha256_file(a.token_dict), "frame_n": len(fr),
              "excluded": {"panel_b": len(panel_b), "controls": len(controls), "members": len(members),
                           "union": len(excluded)},
              "rule": f"estimable = token in >= {CONTROL_MIN_CELLS} of a donor's 100 tumour analysis cells in >= {D_MIN} donors"}
    json.dump({"about": "E2 null genes: first 100 estimable genes of the seed-20261001 draw", "n": len(chosen),
               "draw_order": chosen, **common}, open(f"{a.out_dir}/null_genes_100_estimable.json", "w"), indent=1)
    json.dump({"about": __doc__.strip(), "position_of_nth_estimable": rows[-1]["position"],
               "n_skipped_ineligible": rows[-1]["position"] - len(chosen), "rows": rows, **common},
              open(f"{a.out_dir}/frozen_100_estimable.json", "w"), indent=1)
    print(json.dumps({"n": len(chosen), "walked": len(rows), **common}, indent=1))


if __name__ == "__main__":
    main()
