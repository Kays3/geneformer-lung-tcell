"""Amendment 5: walk the Amendment 4 seed-20260929 draw order and freeze the first N ESTIMABLE genes,
per the human's 2026-09-29 19:31 JST ruling as read by god and confirmed by Stanley -- N counts
estimable genes, not genes drawn; ineligible genes are recorded and skipped, never sent to the GPU.

Eligibility uses the SAME n_token_cells run_isp.py itself computes (tok[g] in s for s in
token_sets[d]), read from the tokenised isp_input datasets directly -- CPU only, no GPU, no model
load, and identical to what the post-run marker will record (the gate-identity assertion
null_analysis.py makes uses this same quantity, per Stanley's item 2).
"""
import argparse
import json
import pickle

from datasets import load_from_disk

CONTROL_MIN_CELLS = 10
D_MIN = 10


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--token-dict", required=True)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--null-genes", required=True, help="Amendment 4's frozen draw (null_genes.json), the "
                                                          "same seed-20260929 order; this script only reads "
                                                          "further down it, never re-draws")
    ap.add_argument("--target-estimable", type=int, required=True)
    ap.add_argument("--cost-fit", required=True, help="fit_cost_model.py's output")
    ap.add_argument("--margin", type=float, default=1.5, help="flat multiplicative margin on the predicted GPU-h")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    tok = pickle.load(open(a.token_dict, "rb"))
    manifest = json.load(open(a.manifest))
    donors = sorted(manifest)
    draw = json.load(open(a.null_genes))["draw_order"]
    token_sets = {d: [set(x) for x in load_from_disk(manifest[d]["isp_input"])["input_ids"]] for d in donors}
    fit = json.load(open(a.cost_fit))["fit"]

    rows = []
    n_estimable = 0
    predicted_seconds = 0.0
    for i, g in enumerate(draw, start=1):
        t = tok[g]
        npos_by_donor = {}
        n_est_donors = 0
        for d in donors:
            n_pos = sum(t in s for s in token_sets[d])
            npos_by_donor[d] = n_pos
            if n_pos >= CONTROL_MIN_CELLS:
                n_est_donors += 1
        is_estimable = n_est_donors >= D_MIN
        rows.append({"gene": g, "position": i, "estimable": is_estimable,
                     "n_estimable_donors": n_est_donors, "npos_by_donor": npos_by_donor})
        if is_estimable:
            n_estimable += 1
            for d in donors:
                npos = npos_by_donor[d]
                if npos > 0:
                    predicted_seconds += (fit["delete"]["intercept_a"] + fit["delete"]["slope_b_per_cell"] * npos)
                    predicted_seconds += (fit["overexpress"]["intercept_a"] + fit["overexpress"]["slope_b_per_cell"] * npos)
        if n_estimable == a.target_estimable:
            break

    out = {
        "about": __doc__.strip(),
        "target_estimable": a.target_estimable,
        "position_of_nth_estimable": i,
        "n_skipped_ineligible": i - a.target_estimable,
        "predicted_gpu_h": predicted_seconds / 3600,
        "margin_factor": a.margin,
        "predicted_gpu_h_with_margin": predicted_seconds / 3600 * a.margin,
        "cost_fit_source": a.cost_fit,
        "rows": rows,
    }
    json.dump(out, open(a.out, "w"), indent=1)
    print(json.dumps({k: v for k, v in out.items() if k != "rows"}, indent=1))


if __name__ == "__main__":
    main()
