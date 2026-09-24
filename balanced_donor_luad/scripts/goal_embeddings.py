"""Phase 6 prep: per-donor goal centroids and per-donor ISP inputs (GPU, ~0.5 h).

For every donor d (held out in fold k, scored by fold k's model):
  tumor_primary          = exact-mean CLS of d's tumour training-pool cells
  normal_adjacent        = exact-mean CLS of d's adjacent-normal pool cells (the registered GOAL, s.3)
  global_normal_adjacent = exact-mean CLS of fold k's TRAINING donors' normal cells (S1 sensitivity)
All via Geneformer EmbExtractor.get_state_embs with the settings the 316M
state embeddings used (CLS, emb_layer 0, exact_mean), in bf16 (dtype_cast)
and nproc=1 (a model already on CUDA cannot survive a forked Dataset.map).
Also writes each donor's ISP input: its 100 tumour ANALYSIS cells.
"""
import argparse
import json
import os
import pickle
import sys
import time


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--folds", required=True)
    p.add_argument("--phase4", required=True, help="phase4 work dir with fold*/fold_record.json")
    p.add_argument("--bf16-bench", required=True)
    p.add_argument("--scripts", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    sys.path[:0] = [a.bf16_bench, a.scripts]
    from dtype_cast import install_dtype_cast
    from model_cache import install_model_cache, verify
    install_dtype_cast("bf16"); install_model_cache()
    from datasets import load_from_disk
    from geneformer import EmbExtractor

    folds = json.load(open(a.folds))
    ds = load_from_disk(a.dataset)
    os.makedirs(a.out, exist_ok=True)
    states = {"state_key": "origin", "start_state": "tumor_primary", "goal_state": "normal_adjacent", "alt_states": []}
    ex = EmbExtractor(model_type="CellClassifier", num_classes=2, emb_mode="cls", filter_data=None,
                      max_ncells=None, emb_layer=0, summary_stat="exact_mean", forward_batch_size=16,
                      nproc=1, model_version="V2")
    t0 = time.time(); manifest = {}
    for f in folds["folds"]:
        k = f["fold"]
        model_dir = json.load(open(os.path.join(a.phase4, f"fold{k}", "fold_record.json")))["model_dir"]
        train = set(f["train"])
        tr_path = os.path.join(a.out, f"fold{k}_train_pool.dataset")
        ds.filter(lambda x: x["individual"] in train, num_proc=1).save_to_disk(tr_path)
        glob_states = ex.get_state_embs(states, model_directory=model_dir, input_data_file=tr_path,
                                        output_directory=a.out, output_prefix=f"fold{k}_global")
        for d in f["test"]:
            dd = os.path.join(a.out, "donors", d.replace("/", "_"))
            os.makedirs(dd, exist_ok=True)
            pool_path = os.path.join(dd, "pool.dataset")
            ds.filter(lambda x: x["individual"] == d, num_proc=1).save_to_disk(pool_path)
            st = ex.get_state_embs(states, model_directory=model_dir, input_data_file=pool_path,
                                   output_directory=dd, output_prefix="donor")
            goals = {"tumor_primary": st["tumor_primary"], "normal_adjacent": st["normal_adjacent"],
                     "global_normal_adjacent": glob_states["normal_adjacent"]}
            pickle.dump(goals, open(os.path.join(dd, "goals.pkl"), "wb"))
            isp_path = os.path.join(dd, "isp_input.dataset")
            ds.filter(lambda x: x["individual"] == d and x["origin"] == "tumor_primary"
                      and int(x["in_analysis100"]) == 1, num_proc=1).save_to_disk(isp_path)
            manifest[d] = {"fold": k, "model_dir": model_dir, "goals": os.path.join(dd, "goals.pkl"),
                           "isp_input": isp_path}
    manifest_path = os.path.join(a.out, "donor_manifest.json")
    json.dump(manifest, open(manifest_path, "w"), indent=1)
    rec = {"n_donors": len(manifest), "seconds": time.time() - t0, **verify()}
    json.dump(rec, open(os.path.join(a.out, "goal_embeddings_record.json"), "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
