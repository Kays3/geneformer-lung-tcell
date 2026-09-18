"""
THROWAWAY calibration fine-tune, 316M second-stage step-0 (2026-09-18).

Produces a structurally-valid Geneformer-V2-316M classifier checkpoint
FAST (small max_ncells cap) purely to time one real targeted-panel ISP
unit before committing GPU time to the full ~2h bf16 fine-tune. This is
NOT the real 316M classifier -- its accuracy is not evaluated or
reported anywhere. Mirrors /srv/lab/KD/sclc_luad_normal_htan_finetune/
scripts/run_finetune.py exactly except: (a) base model swapped to
Geneformer-V2-316M, (b) bf16 training, (c) max_ncells capped for speed,
(d) output directory is a throwaway path under bf16_bench/runs/, never
touching /srv/lab/KD (Pam's committed 104M fine-tune outputs).

Committed after the fact (2026-09-19), recovered from a stale ts1 working
checkout during a checkout cleanup -- it is a step-0 SIZING tool, not part
of the science: its only purpose is producing one real targeted-panel ISP
unit's wall-clock/GPU numbers fast enough to size the real ~2h run before
committing GPU time to it, which is what the 316M sizing report (that
approved the full run) was built on. The checkpoint it trains is never
evaluated for accuracy anywhere and must never be cited as a 316M result --
that is run_finetune_316m.py's job, the real (committed, uncapped) fine-tune
this task's actual RESULTS_BF16_316M.md numbers come from. Kept committed,
not deleted, so the next model-size variant has this sizing step ready to
reuse; its outputs stay under the gitignored runs/, so nothing bulky lands.
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

SEED = 43

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
DISEASE_STATES = [SCLC, LUAD, NORMAL]

HOME = Path.home()
SOURCE_WORK_DIR = HOME / "workspace/KD/sclc_luad_normal_htan_finetune"
SOURCE_TOKENIZED = SOURCE_WORK_DIR / "data/sclc_luad_normal_htan_tcells.dataset"

WORK_DIR = (
    HOME
    / "workspace/geneformer-lung-tcell/sclc_validation/bf16_bench/runs/_calibration_316m_finetune"
)
RUN_DIR = WORK_DIR / "runs"
TABLE_DIR = WORK_DIR / "tables"

BASE_MODEL_DIR = (
    HOME / "workspace/geneformer-uv-starter/Geneformer/Geneformer-V2-316M"
)
OUTPUT_PREFIX = "sclc_luad_normal_htan_316m_calib"

CALIB_MAX_NCELLS = 400  # throwaway cap: speed only, not the real fine-tune


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
        max_ncells=CALIB_MAX_NCELLS,
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

    saved_models = sorted(RUN_DIR.glob(f"*geneformer_cellClassifier_{OUTPUT_PREFIX}/ksplit1"))
    assert saved_models, "No saved classifier model found"
    model_path = saved_models[-1]
    (RUN_DIR / "MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt").write_text(str(model_path))
    print("\nCalibration-only model (NOT the real 316M classifier):", model_path)

    payload = {
        "model_path": str(model_path),
        "base_model_dir": str(BASE_MODEL_DIR),
        "seed": SEED,
        "disease_states": DISEASE_STATES,
        "source_tokenized_dataset": str(SOURCE_TOKENIZED),
        "calibration_only": True,
        "calib_max_ncells": CALIB_MAX_NCELLS,
        "note": "THROWAWAY checkpoint for ISP timing calibration only. Not evaluated for accuracy. Do not cite as the 316M classifier result.",
        "eval_metrics": {
            k: normalize_metric(v) for k, v in eval_metrics.items() if k != "conf_matrix"
        },
    }
    with open(TABLE_DIR / "calibration_finetune_metadata.json", "w") as f:
        json.dump(payload, f, indent=2)

    print("\nDone.")
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
