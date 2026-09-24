"""Phase 6: in-silico perturbation, BOTH operations, per (gene, donor), on one host.

For each gene assigned to this host (controls/pre_gpu_host_map.csv, Amendment
3.1: split by gene, so each host runs delete AND overexpress for its genes),
for each operation, for each of the 43 donors:
  InSilicoPerturber(perturb_type=op, genes_to_perturb=[gene],
                    cell_states_to_model = start tumor_primary -> goal normal_adjacent
                                           (the DONOR's own centroid), alt global_normal_adjacent (S1),
                    state_embs_dict = that donor's goals.pkl, model = the donor's fold model,
                    emb_mode cls, emb_layer 0, nproc 1, bf16)
on the donor's 100 tumour analysis cells. Donor values are computed later from
these raw per-cell pickles (never Shift_to_goal_end).

A (gene, donor) with no token-positive cell is recorded as `no_token_cells`
without calling the perturber (Geneformer raises on an empty filter).
Resumable: one completion marker per (op, gene, donor). The model is loaded
once per fold (model_cache) and its checksum re-verified after every gene.
"""
import argparse
import json
import os
import pickle
import sys
import time

import pandas as pd

STATES = {"state_key": "origin", "start_state": "tumor_primary", "goal_state": "normal_adjacent",
          "alt_states": ["global_normal_adjacent"]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="goal_embeddings donor_manifest.json")
    p.add_argument("--host-map", required=True)
    p.add_argument("--host", required=True)
    p.add_argument("--exclude", default="ENSG00000133639", help="genes not run (BTG1: NOT_ESTIMABLE_CONTROLS)")
    p.add_argument("--genes", default=None, help="comma list overriding the host map (probe)")
    p.add_argument("--donors", default=None, help="comma list overriding all donors (probe)")
    p.add_argument("--token-dict", required=True)
    p.add_argument("--bf16-bench", required=True)
    p.add_argument("--scripts", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    sys.path[:0] = [a.bf16_bench, a.scripts]
    from dtype_cast import install_dtype_cast
    from model_cache import install_model_cache, verify
    install_dtype_cast("bf16"); install_model_cache()
    from datasets import load_from_disk
    from geneformer import InSilicoPerturber

    manifest = json.load(open(a.manifest))
    tok = pickle.load(open(a.token_dict, "rb"))
    if a.genes:
        genes = a.genes.split(",")
    else:
        hm = pd.read_csv(a.host_map)
        genes = sorted(set(hm[hm.host == a.host].ensembl_id) - set(a.exclude.split(",")))
    donors = a.donors.split(",") if a.donors else sorted(manifest)
    token_sets = {d: [set(x) for x in load_from_disk(manifest[d]["isp_input"])["input_ids"]] for d in donors}
    goals = {d: pickle.load(open(manifest[d]["goals"], "rb")) for d in donors}
    log = open(os.path.join(a.out, f"run_log_{a.host}.jsonl"), "a") if os.path.isdir(a.out) else None
    os.makedirs(a.out, exist_ok=True)
    log = log or open(os.path.join(a.out, f"run_log_{a.host}.jsonl"), "a")
    for g in genes:
        t_gene = time.time()
        for op in ("delete", "overexpress"):
            for d in donors:
                safe = d.replace("/", "_")
                od = os.path.join(a.out, op, g)
                marker = os.path.join(od, f"{safe}.complete.json")
                if os.path.exists(marker):
                    continue
                os.makedirs(od, exist_ok=True)
                n_pos = sum(tok[g] in s for s in token_sets[d])
                t0 = time.time()
                if n_pos == 0:
                    rec = {"gene": g, "op": op, "donor": d, "status": "no_token_cells", "n_token_cells": 0}
                else:
                    isp = InSilicoPerturber(perturb_type=op, genes_to_perturb=[g], combos=0, anchor_gene=None,
                                            model_type="CellClassifier", num_classes=2, emb_mode="cls",
                                            filter_data=None, cell_states_to_model=STATES, state_embs_dict=goals[d],
                                            max_ncells=None, emb_layer=0, forward_batch_size=128, nproc=1,
                                            model_version="V2", clear_mem_ncells=1000)
                    isp.perturb_data(model_directory=manifest[d]["model_dir"], input_data_file=manifest[d]["isp_input"],
                                     output_directory=od, output_prefix=safe)
                    rec = {"gene": g, "op": op, "donor": d, "status": "done", "n_token_cells": int(n_pos)}
                rec.update({"seconds": time.time() - t0, "host": a.host, "fold": manifest[d]["fold"]})
                json.dump(rec, open(marker, "w"))
                log.write(json.dumps(rec) + "\n"); log.flush()
        v = verify()
        log.write(json.dumps({"gene_done": g, "seconds": time.time() - t_gene, **v}) + "\n"); log.flush()


if __name__ == "__main__":
    main()
