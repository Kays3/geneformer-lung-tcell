"""
Recovery script (2026-09-18): the real 316M fine-tune (run_finetune_316m.py)
completed training successfully -- checkpoint saved -- but crashed during
post-training eval-metrics computation on a Geneformer library bug (see
eval_tie_fix.py). Rather than re-running the ~2h training, this reuses the
already-saved checkpoint and already-prepared train/test datasets (both
produced by prepare_data() before the crash) to compute held-out test
metrics via evaluate_saved_model(), with the tie-fix applied.
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from eval_tie_fix import install_eval_tie_fix, reset_tie_counter, tie_stats  # noqa: E402

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
DISEASE_STATES = [SCLC, LUAD, NORMAL]

HOME = Path.home()
WORK_DIR = HOME / "workspace/geneformer-lung-tcell/sclc_validation/bf16_bench/runs/316m_bf16_finetune"
RUN_DIR = WORK_DIR / "runs"
TABLE_DIR = WORK_DIR / "tables"
OUTPUT_PREFIX = "sclc_luad_normal_htan_316m"

MODEL_PATH = RUN_DIR / "260918_geneformer_cellClassifier_sclc_luad_normal_htan_316m/ksplit1"


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
    os.environ["WANDB_DISABLED"] = "true"
    os.environ.setdefault("MPLBACKEND", "Agg")

    assert MODEL_PATH.exists(), f"Missing trained checkpoint: {MODEL_PATH}"
    install_eval_tie_fix()
    reset_tie_counter()

    from geneformer import Classifier

    classifier = Classifier(
        classifier="cell",
        cell_state_dict={"state_key": "disease", "states": DISEASE_STATES},
        filter_data=None,
        training_args={
            "num_train_epochs": 1,
            "learning_rate": 5e-5,
            "per_device_train_batch_size": 8,
            "seed": 43,
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

    (RUN_DIR / "MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt").write_text(str(MODEL_PATH))

    test_metrics = classifier.evaluate_saved_model(
        model_directory=str(MODEL_PATH),
        id_class_dict_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_id_class_dict.pkl"),
        test_data_file=str(RUN_DIR / f"{OUTPUT_PREFIX}_labeled_test.dataset"),
        output_directory=str(RUN_DIR),
        output_prefix=f"{OUTPUT_PREFIX}_test",
    )

    conf_mat = test_metrics.get("conf_matrix")
    if conf_mat is not None:
        conf_mat.to_csv(TABLE_DIR / "test_confusion_matrix.csv")

    test_metrics_clean = {
        k: normalize_metric(v) for k, v in test_metrics.items() if k != "conf_matrix"
    }
    payload = {
        "model_path": str(MODEL_PATH),
        "seed": 43,
        "disease_states": DISEASE_STATES,
        "bf16": True,
        "recovered_from_eval_crash": True,
        "eval_tie_fix_applied": True,
        "eval_tie_stats": tie_stats(),
        "test_metrics": test_metrics_clean,
    }
    with open(TABLE_DIR / "test_metrics.json", "w") as f:
        json.dump(payload, f, indent=2)

    print("\nRecovered test metrics (no retraining):")
    print(json.dumps(payload, indent=2))
    print(f"\nEVAL_TIE_STATS {tie_stats()}")


if __name__ == "__main__":
    main()
