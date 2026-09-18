"""
Real Geneformer-V2-316M fine-tune for the second-stage bf16 replication
(2026-09-18). Same recipe, seed, tokenized dataset, and donor-based split
as the committed 104M fine-tune
(/srv/lab/KD/sclc_luad_normal_htan_finetune/scripts/run_finetune.py,
READ-ONLY reference, never modified), except: base model is
Geneformer-V2-316M, training runs in bf16, and this is NOT capped by
max_ncells (this is the real classifier, not a calibration throwaway).
Output goes to a fresh directory under bf16_bench/runs/, never touching
/srv/lab/KD (the 104M fine-tune's committed outputs).

Held-out test-set accuracy/macro-F1 (from classifier.evaluate_saved_model)
is written to tables/test_metrics.json -- this is the sanity-abort input
for the stage chain (god's ruling: abort the chain if macro_f1 < 0.60,
log the metrics either way, and report them in RESULTS_BF16_316M.md next
to the ISP comparison per the plan's hard boundary on 316M classifier
quality).

Applies eval_tie_fix.py (2026-09-18, added after a live crash on the first
full run of this script): geneformer.evaluation_utils.vote() returns the
string "tie" on an exact logit tie, which crashes sklearn's
confusion_matrix once any tie occurs -- and bf16 TRAINING (unlike the
104M classifier, trained fp32) makes exact ties plausible on real data.
See eval_tie_fix.py's docstring for the full root cause. Tie counts are
recorded in test_metrics.json's eval_tie_stats.
"""
from __future__ import annotations

import json
import os
import pickle
import random
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eval_tie_fix import install_eval_tie_fix, reset_tie_counter, tie_stats  # noqa: E402

SEED = 43

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
DISEASE_STATES = [SCLC, LUAD, NORMAL]

HOME = Path.home()
SOURCE_WORK_DIR = HOME / "workspace/KD/sclc_luad_normal_htan_finetune"
SOURCE_TOKENIZED = SOURCE_WORK_DIR / "data/sclc_luad_normal_htan_tcells.dataset"

WORK_DIR = HOME / "workspace/geneformer-lung-tcell/sclc_validation/bf16_bench/runs/316m_bf16_finetune"
RUN_DIR = WORK_DIR / "runs"
TABLE_DIR = WORK_DIR / "tables"

BASE_MODEL_DIR = HOME / "workspace/geneformer-uv-starter/Geneformer/Geneformer-V2-316M"
OUTPUT_PREFIX = "sclc_luad_normal_htan_316m"


def normalize_metric(value):
    if isinstance(value, (np.floating, np.integer)):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, list):
        return [normalize_metric(v) for v in value]
    if isinstance(value, dict):
        return {str(k): normalize_metric(v) for k, v in value.items()}
    return value


def main() -> None:
    random.seed(SEED)
    np.random.seed(SEED)
    os.environ["WANDB_DISABLED"] = "true"
    os.environ.setdefault("MPLBACKEND", "Agg")
    install_eval_tie_fix()
    reset_tie_counter()

    for d in [RUN_DIR, TABLE_DIR]:
        d.mkdir(parents=True, exist_ok=True)

    print("Python:", sys.executable)
    print("Work dir:", WORK_DIR)
    print("Source tokenized dataset (read-only):", SOURCE_TOKENIZED)
    print("Base model:", BASE_MODEL_DIR)

    assert SOURCE_TOKENIZED.exists(), f"Missing tokenized dataset: {SOURCE_TOKENIZED}"
    assert BASE_MODEL_DIR.exists(), f"Missing base model: {BASE_MODEL_DIR}"

    from datasets import load_from_disk

    ds = load_from_disk(str(SOURCE_TOKENIZED))
    meta = pd.DataFrame(
        {
            "row_index": list(range(len(ds))),
            "cell_id": ds["cell_id"],
            "individual": ds["individual"],
            "celltype": ds["celltype"],
            "disease": ds["disease"],
            "split": ds["split"],
        }
    )
    meta.to_csv(TABLE_DIR / "tokenized_cell_metadata.csv", index=False)

    train_ids = sorted(meta.loc[meta["split"] == "train", "individual"].unique().tolist())
    eval_ids = sorted(meta.loc[meta["split"] == "eval", "individual"].unique().tolist())
    test_ids = sorted(meta.loc[meta["split"] == "test", "individual"].unique().tolist())

    from geneformer import Classifier

    classifier = Classifier(
        classifier="cell",
        cell_state_dict={"state_key": "disease", "states": DISEASE_STATES},
        filter_data=None,
        training_args={
            "num_train_epochs": 1,
            "learning_rate": 5e-5,
            "per_device_train_batch_size": 8,
            "seed": SEED,
            "save_strategy": "epoch",
            "logging_steps": 20,
            "report_to": "none",
            "bf16": True,
        },
        max_ncells=None,
        freeze_layers=6,
        num_crossval_splits=1,
        forward_batch_size=16,
        nproc=4,
        model_version="V2",
    )

    classifier.prepare_data(
        input_data_file=str(SOURCE_TOKENIZED),
        output_directory=str(RUN_DIR),
        output_prefix=OUTPUT_PREFIX,
        split_id_dict={
            "attr_key": "individual",
            "train": train_ids + eval_ids,
            "test": test_ids,
        },
    )

    eval_metrics = classifier.validate(
        model_directory=str(BASE_MODEL_DIR),
        prepared_input_data_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_labeled_train.dataset"),
        id_class_dict_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_id_class_dict.pkl"),
        output_directory=str(RUN_DIR),
        output_prefix=OUTPUT_PREFIX,
        split_id_dict={
            "attr_key": "individual",
            "train": train_ids,
            "eval": eval_ids,
        },
        n_hyperopt_trials=0,
    )
    eval_split_tie_stats = tie_stats()
    reset_tie_counter()

    saved_models = sorted(RUN_DIR.glob(f"*geneformer_cellClassifier_{OUTPUT_PREFIX}/ksplit1"))
    assert saved_models, "No saved classifier model found"
    model_path = saved_models[-1]
    (RUN_DIR / "MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt").write_text(str(model_path))
    print("\nSelected model:", model_path)

    test_metrics = classifier.evaluate_saved_model(
        model_directory=str(model_path),
        id_class_dict_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_id_class_dict.pkl"),
        test_data_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_labeled_test.dataset"),
        output_directory=str(RUN_DIR),
        output_prefix=f"{OUTPUT_PREFIX}_test",
    )

    conf_mat = test_metrics.get("conf_matrix")
    if conf_mat is not None:
        conf_mat.to_csv(TABLE_DIR / "test_confusion_matrix.csv")

    with open(RUN_DIR / f"{OUTPUT_PREFIX}_id_class_dict.pkl", "rb") as f:
        id_class_dict = pickle.load(f)
    with open(TABLE_DIR / "id_class_dict.json", "w") as f:
        json.dump({str(k): v for k, v in id_class_dict.items()}, f, indent=2)

    eval_metrics_clean = {
        k: normalize_metric(v) for k, v in eval_metrics.items() if k != "conf_matrix"
    }
    test_metrics_clean = {
        k: normalize_metric(v) for k, v in test_metrics.items() if k != "conf_matrix"
    }
    payload = {
        "model_path": str(model_path),
        "base_model_dir": str(BASE_MODEL_DIR),
        "seed": SEED,
        "disease_states": DISEASE_STATES,
        "source_tokenized_dataset": str(SOURCE_TOKENIZED),
        "bf16": True,
        "eval_split_tie_stats": eval_split_tie_stats,
        "test_split_tie_stats": tie_stats(),
        "eval_metrics": eval_metrics_clean,
        "test_metrics": test_metrics_clean,
    }
    with open(TABLE_DIR / "test_metrics.json", "w") as f:
        json.dump(payload, f, indent=2)

    print("\nDone.")
    print(json.dumps(payload, indent=2))
    print(f"\nEVAL_SPLIT_TIE_STATS {eval_split_tie_stats}")
    print(f"TEST_SPLIT_TIE_STATS {tie_stats()}")


if __name__ == "__main__":
    main()
