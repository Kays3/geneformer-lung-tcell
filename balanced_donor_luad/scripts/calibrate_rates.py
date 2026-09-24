"""Calibration run (human-approved 2026-09-24, <= 0.15 GPU-h hard stop).

Measures, and records ONLY, two rates on the thinkstation1 GB10 with the
model the human directed on 2026-09-24 ("all future work use bf16 model
geneformer 316M"): Geneformer-V2-316M, bf16. Training uses the Trainer's bf16
flag. ISP uses the lab's dtype_cast wrapper (model cast to bfloat16 at load)
on the vendored Geneformer checkout f45a6c7 + geneformer-bf16-export.patch,
which is the same stack the measured 316M benchmark used:
  1. fine-tune throughput: training cells/s for 1 epoch on 1,000 cells;
  2. ISP throughput: seconds per 300-cell unit, deletion and overexpression.

This is not an experiment. The model is a throwaway. Its evaluation metrics
are never produced: train_all_data trains with no evaluation split. No analysis cell (in_analysis100 == 1) is
used anywhere. The ISP gene (B2M) is in neither panel. The ISP outputs are
deleted unread; only the elapsed time is kept. Run the script under
`timeout 540` so the budget is enforced by the OS, not by this code.
"""
import argparse
import glob
import json
import os
import shutil
import time

import numpy as np

SEED = 20260924
N_TRAIN, N_ISP = 1000, 300
ISP_GENE = "ENSG00000166710"  # B2M; ubiquitous in T cells; not in panel A or B
PREFIX = "calib"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--base-model", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    os.environ["WANDB_DISABLED"] = "true"
    t_start = time.time()
    rec = {"purpose": "rate calibration only; model and ISP outputs discarded unread",
           "seed": SEED, "precision": "bf16", "model": os.path.basename(a.base_model.rstrip("/"))}

    from datasets import load_from_disk
    ds = load_from_disk(a.dataset)
    pool = ds.filter(lambda x: int(x["in_analysis100"]) == 0, num_proc=4)
    donors = sorted(set(pool["individual"]))
    rng = np.random.default_rng(SEED)
    rng.shuffle(donors)
    eval_donors, train_donors = donors[:4], donors[4:]
    tr = pool.filter(lambda x: x["individual"] in set(train_donors), num_proc=4).shuffle(seed=SEED).select(range(N_TRAIN))
    os.makedirs(a.work, exist_ok=True)
    calib_path = os.path.join(a.work, "calib_input.dataset")
    tr.save_to_disk(calib_path)
    rec["train_cells"] = len(tr)
    rec["median_tokens_train"] = float(np.median(tr["length"]))

    from geneformer import Classifier
    cc = Classifier(classifier="cell",
                    cell_state_dict={"state_key": "origin", "states": ["tumor_primary", "normal_adjacent"]},
                    filter_data=None,
                    training_args={"num_train_epochs": 1, "learning_rate": 5e-5,
                                   "per_device_train_batch_size": 8, "seed": 43,
                                   "save_strategy": "no", "logging_steps": 20, "report_to": "none",
                                   "bf16": True},
                    max_ncells=None, freeze_layers=6, num_crossval_splits=1,
                    forward_batch_size=16, nproc=4, model_version="V2")
    cc.prepare_data(input_data_file=calib_path, output_directory=a.work, output_prefix=PREFIX)
    # prepare_data splits 90/10 by default here; only the train part is used and
    # the 10% part is never evaluated. The rate uses the actual trained count.
    from datasets import load_from_disk as _lfd
    train_file = os.path.join(a.work, f"{PREFIX}_labeled_train.dataset")
    n_trained = len(_lfd(train_file))
    rec["train_cells"] = n_trained
    t0 = time.time()
    # train_all_data: no evaluation split at all, so no accuracy is ever produced
    cc.train_all_data(model_directory=a.base_model,
                      prepared_input_data_file=train_file,
                      id_class_dict_file=os.path.join(a.work, f"{PREFIX}_id_class_dict.pkl"),
                      output_directory=a.work, output_prefix=PREFIX, save_eval_output=False)
    rec["finetune_seconds"] = time.time() - t0
    rec["train_cells_per_second"] = n_trained / rec["finetune_seconds"]

    configs = sorted(glob.glob(os.path.join(a.work, "**", "config.json"), recursive=True))
    configs = [c for c in configs if os.path.exists(os.path.join(os.path.dirname(c), "model.safetensors"))
               or os.path.exists(os.path.join(os.path.dirname(c), "pytorch_model.bin"))]
    model_dir = os.path.dirname(configs[-1])
    isp_src = pool.filter(lambda x: x["origin"] == "tumor_primary" and x["individual"] in set(eval_donors),
                          num_proc=4)
    isp_src = isp_src.select(range(min(N_ISP, len(isp_src))))
    isp_path = os.path.join(a.work, "isp_input.dataset")
    isp_src.save_to_disk(isp_path)
    rec["isp_cells"] = len(isp_src)

    from geneformer import InSilicoPerturber
    from dtype_cast import install_dtype_cast
    install_dtype_cast("bf16")
    for op in ("delete", "overexpress"):
        outdir = os.path.join(a.work, f"isp_{op}")
        os.makedirs(outdir, exist_ok=True)
        isp = InSilicoPerturber(perturb_type=op, genes_to_perturb=[ISP_GENE], combos=0, anchor_gene=None,
                                model_type="CellClassifier", num_classes=2, emb_mode="cls",
                                filter_data=None, cell_states_to_model=None, max_ncells=None,
                                emb_layer=0, forward_batch_size=128, nproc=1, model_version="V2",
                                clear_mem_ncells=1000)
        t0 = time.time()
        isp.perturb_data(model_directory=model_dir, input_data_file=isp_path,
                         output_directory=outdir, output_prefix=f"{PREFIX}_{op}")
        rec[f"isp_{op}_seconds"] = time.time() - t0
        rec[f"isp_{op}_seconds_per_300_cells"] = rec[f"isp_{op}_seconds"] * 300 / len(isp_src)
        shutil.rmtree(outdir)  # ISP outputs deleted unread

    rec["total_wall_seconds"] = time.time() - t_start
    rec["gpu_hours_used_upper_bound"] = rec["total_wall_seconds"] / 3600
    json.dump(rec, open(a.out, "w"), indent=1)
    print(json.dumps(rec, indent=1))


if __name__ == "__main__":
    main()
