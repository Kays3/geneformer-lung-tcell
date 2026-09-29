"""Phase 4: fine-tune one cross-fitting fold (Geneformer-V2-316M, bf16) and score its held-out donors.

Recipe (registration s.2 + Amendment 1, fixed): 1 epoch, lr 5e-5, batch 8,
freeze_layers 6, seed 43, label `origin` (tumor_primary vs normal_adjacent),
no hyperparameter search, eval donors reported but never used for selection.
Follows sclc_validation/bf16_bench/run_finetune_316m.py: Classifier ->
prepare_data(split by donor) -> validate -> evaluate_saved_model, with
eval_tie_fix installed (bf16 logits can tie exactly) and tie counts recorded.

Cells:
  train + eval donors : their training-pool cells (<= 300 per tissue);
  held-out (test) donors : ONLY their analysis cells (exactly 100 per tissue),
                           so every held-out donor is scored on a balanced set.
Per-cell predictions for the held-out donors are saved with donor and cell ids
for the classifier gate (scripts/classifier_gate.py).
"""
import argparse
import glob
import json
import os
import pickle
import sys
import time

SEED = 43
STATES = ["tumor_primary", "normal_adjacent"]



def single_split_value(metrics, key):
    """Classifier.validate returns one value PER k-split as a list; this run has
    exactly one split (num_crossval_splits=1), so the list must have length 1."""
    if key not in metrics:
        return None
    v = metrics[key]
    if isinstance(v, (list, tuple)):
        if len(v) != 1:
            raise ValueError(f"expected one k-split value for {key}, got {len(v)}")
        v = v[0]
    return float(v)

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dataset", required=True)
    p.add_argument("--folds", required=True)
    p.add_argument("--fold", type=int, required=True)
    p.add_argument("--base-model", required=True)
    p.add_argument("--work", required=True)
    p.add_argument("--bf16-bench", required=True, help="path to sclc_validation/bf16_bench (eval_tie_fix)")
    a = p.parse_args()
    os.environ["WANDB_DISABLED"] = "true"
    sys.path.insert(0, a.bf16_bench)
    from eval_tie_fix import install_eval_tie_fix, reset_tie_counter, tie_stats
    from datasets import load_from_disk
    from geneformer import Classifier

    folds = json.load(open(a.folds))
    f = folds["folds"][a.fold]
    train, ev, test = set(f["train"]), set(f["eval"]), set(f["test"])
    work = os.path.join(a.work, f"fold{a.fold}")
    os.makedirs(work, exist_ok=True)
    t0 = time.time()
    ds = load_from_disk(a.dataset)
    fold_ds = ds.filter(lambda x: x["individual"] in train or x["individual"] in ev
                        or (x["individual"] in test and int(x["in_analysis100"]) == 1), num_proc=4)
    inp = os.path.join(work, "fold_input.dataset")
    fold_ds.save_to_disk(inp)

    install_eval_tie_fix()
    cc = Classifier(classifier="cell", cell_state_dict={"state_key": "origin", "states": STATES},
                    filter_data=None,
                    training_args={"num_train_epochs": 1, "learning_rate": 5e-5,
                                   "per_device_train_batch_size": 8, "seed": SEED,
                                   "save_strategy": "epoch", "logging_steps": 20,
                                   "report_to": "none", "bf16": True},
                    max_ncells=None, freeze_layers=6, num_crossval_splits=1,
                    forward_batch_size=16, nproc=4, model_version="V2")
    prefix = f"bdl_fold{a.fold}"
    cc.prepare_data(input_data_file=inp, output_directory=work, output_prefix=prefix,
                    split_id_dict={"attr_key": "individual", "train": sorted(train | ev), "test": sorted(test)})
    t_train = time.time()
    eval_metrics = cc.validate(model_directory=a.base_model,
                               prepared_input_data_file=os.path.join(work, f"{prefix}_labeled_train.dataset"),
                               id_class_dict_file=os.path.join(work, f"{prefix}_id_class_dict.pkl"),
                               output_directory=work, output_prefix=prefix,
                               split_id_dict={"attr_key": "individual", "train": sorted(train), "eval": sorted(ev)},
                               n_hyperopt_trials=0)
    train_seconds = time.time() - t_train
    eval_ties = tie_stats(); reset_tie_counter()
    model_dir = sorted(glob.glob(os.path.join(work, f"*geneformer_cellClassifier_{prefix}", "ksplit1")))[-1]
    t_test = time.time()
    test_metrics = cc.evaluate_saved_model(model_directory=model_dir,
                                           id_class_dict_file=os.path.join(work, f"{prefix}_id_class_dict.pkl"),
                                           test_data_file=os.path.join(work, f"{prefix}_labeled_test.dataset"),
                                           output_directory=work, output_prefix=f"{prefix}_test",
                                           predict=True, predict_metadata=["individual", "cell_id"])
    test_seconds = time.time() - t_test
    rec = {"fold": a.fold, "model_dir": model_dir, "n_train_donors": len(train), "n_eval_donors": len(ev),
           "n_test_donors": len(test), "n_fold_cells": len(fold_ds),
           "eval_macro_f1": single_split_value(eval_metrics, "macro_f1"),
           "eval_acc": single_split_value(eval_metrics, "acc"),
           "eval_tie_stats": eval_ties, "test_tie_stats": tie_stats(),
           "train_plus_eval_seconds": train_seconds, "test_seconds": test_seconds,
           "total_seconds": time.time() - t0,
           "pred_dict": os.path.join(work, f"{prefix}_test_pred_dict.pkl"),
           "id_class_dict": os.path.join(work, f"{prefix}_id_class_dict.pkl")}
    json.dump(rec, open(os.path.join(work, "fold_record.json"), "w"), indent=1, default=str)
    print(json.dumps(rec, indent=1, default=str))


if __name__ == "__main__":
    main()
