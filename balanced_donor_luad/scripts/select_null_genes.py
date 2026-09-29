"""Amendment 4: freeze the genome-wide null-gene draw before any GPU output.

Sampling frame = every gene in the pinned V2 token dictionary MINUS the four
special tokens, MINUS Panel A (15), MINUS Panel B (39 listed, incl. the 5 not
in the dictionary -- excluded harmlessly since they can't be drawn anyway),
MINUS every gene ever used as a matched control in this analysis (design.json
strata, de-duplicated). Simple random sample without replacement, seed fixed
before any output exists, order is the draw order (first-drawn = first-run).
"""
import argparse
import hashlib
import json
import pickle
import random


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--panel-a", required=True)
    ap.add_argument("--panel-b", required=True)
    ap.add_argument("--design", required=True)
    ap.add_argument("--n-draw", type=int, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tok = pickle.load(open(a.token_dict, "rb"))
    special = {"<pad>", "<mask>", "<cls>", "<eos>"}
    universe = sorted(g for g in tok if g not in special)

    panel_a = json.load(open(a.panel_a))
    panel_a_ids = {g["ensembl_id"] for g in panel_a["genes"]}

    panel_b = json.load(open(a.panel_b))
    panel_b_ids = set(panel_b["genes"]) if isinstance(panel_b["genes"], list) and panel_b["genes"] and isinstance(panel_b["genes"][0], str) else {g["ensembl_id"] for g in panel_b["genes"]}

    design = json.load(open(a.design))
    control_ids = set()
    member_ids = set()
    for st in design["strata"].values():
        control_ids.update(st.get("controls", []))
        member_ids.update(st.get("members", []))

    excluded = panel_a_ids | panel_b_ids | control_ids | member_ids
    frame = [g for g in universe if g not in excluded]

    rng = random.Random(a.seed)
    draw = frame.copy()
    rng.shuffle(draw)
    draw = draw[: a.n_draw]

    out = {
        "about": __doc__.strip(),
        "token_dict_path": a.token_dict,
        "token_dict_sha256": sha256_file(a.token_dict),
        "token_dict_n_genes": len(tok),
        "universe_n": len(universe),
        "excluded_n": len(excluded),
        "excluded_breakdown": {
            "panel_a": len(panel_a_ids),
            "panel_b": len(panel_b_ids),
            "design_controls": len(control_ids),
            "design_members": len(member_ids),
            "union": len(excluded),
        },
        "frame_n": len(frame),
        "seed": a.seed,
        "n_draw": a.n_draw,
        "draw_order": draw,
    }
    json.dump(out, open(a.out, "w"), indent=1)
    print(f"universe={len(universe)} excluded={len(excluded)} frame={len(frame)} drew={len(draw)}")


if __name__ == "__main__":
    main()
