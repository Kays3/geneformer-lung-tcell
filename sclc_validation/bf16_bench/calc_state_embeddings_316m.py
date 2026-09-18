"""
Compute 316M-specific disease-state centroid embeddings for the targeted-
panel ISP calibration (2026-09-18). Recipe copied exactly from
KD/sclc_luad_normal_htan_heldout_allgene_perturbation/scripts/
run_heldout_allgene.py's calculate_state_embeddings(), pointed at whatever
MODEL_DIRECTORY env var gives (the throwaway calibration classifier for
this step-0 sizing pass; the real fine-tuned 316M classifier later).
Reads ALLGENE_ROOT's train_reference.dataset READ-ONLY; writes only to
OUTPUT_DIR (never touches ALLGENE_ROOT).
"""
from __future__ import annotations

import os
import pickle
import time
from pathlib import Path

HOME = Path.home()
ALLGENE_ROOT = Path(
    os.environ.get(
        "SCLC_PERTURBATION_ROOT",
        HOME / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation",
    )
)
TRAIN_DATASET = ALLGENE_ROOT / "data/train_reference.dataset"

MODEL_DIRECTORY = Path(os.environ["MODEL_DIRECTORY"])
OUTPUT_DIR = Path(os.environ["OUTPUT_DIR"])
OUTPUT_PREFIX = os.environ.get("OUTPUT_PREFIX", "state_centroids")

SCLC = "small cell lung carcinoma"
LUAD = "lung adenocarcinoma"
NORMAL = "normal"
STATES = [SCLC, LUAD, NORMAL]


def canonical_states(start_state: str) -> dict:
    others = [state for state in STATES if state != start_state]
    return {
        "state_key": "disease",
        "start_state": start_state,
        "goal_state": others[0],
        "alt_states": [others[1]],
    }


def main() -> None:
    assert TRAIN_DATASET.exists(), f"Missing {TRAIN_DATASET}"
    assert MODEL_DIRECTORY.exists(), f"Missing {MODEL_DIRECTORY}"
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    from geneformer import EmbExtractor

    t0 = time.time()
    extractor = EmbExtractor(
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",
        filter_data=None,
        max_ncells=None,
        emb_layer=0,
        summary_stat="exact_mean",
        forward_batch_size=16,
        nproc=4,
        model_version="V2",
    )
    state_embs = extractor.get_state_embs(
        canonical_states(SCLC),
        model_directory=str(MODEL_DIRECTORY),
        input_data_file=str(TRAIN_DATASET),
        output_directory=str(OUTPUT_DIR),
        output_prefix=OUTPUT_PREFIX,
    )
    elapsed = time.time() - t0

    out_pkl = OUTPUT_DIR / f"{OUTPUT_PREFIX}.pkl"
    with out_pkl.open("wb") as f:
        pickle.dump(state_embs, f)

    print("STATE_EMB_ELAPSED_SECONDS", elapsed)
    print("STATE_EMB_STATES", sorted(state_embs))
    print("STATE_EMB_FILE", out_pkl)


if __name__ == "__main__":
    main()
