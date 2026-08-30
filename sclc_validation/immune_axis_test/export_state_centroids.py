#!/usr/bin/env python3
"""Export the trusted Geneformer state-centroid pickle as portable NumPy.

Run this on a compute host where the original pickle and PyTorch environment
exist. The output contains only three aggregate vectors; it contains no cells,
donor identifiers, or model weights. Never run this on an untrusted pickle.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import pickle
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

STATE_ALIASES = {
    "normal": "normal",
    "small cell lung carcinoma": "sclc",
    "sclc": "sclc",
    "lung adenocarcinoma": "luad",
    "luad": "luad",
}
REQUIRED_STATES = ("normal", "sclc", "luad")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def as_vector(value: object) -> np.ndarray:
    if hasattr(value, "detach"):
        value = value.detach()
    if hasattr(value, "cpu"):
        value = value.cpu()
    if hasattr(value, "numpy"):
        value = value.numpy()
    vector = np.asarray(value, dtype=np.float64).squeeze()
    if vector.ndim != 1 or vector.size < 2:
        raise ValueError(f"Expected one centroid vector, got shape {vector.shape}")
    if not np.isfinite(vector).all():
        raise ValueError("Centroid contains NaN or infinite values")
    return vector


def export_centroids(input_path: Path, output_path: Path) -> dict:
    # This is an internal, trusted experiment artifact. Pickle can execute code;
    # do not relax the provenance requirement or accept arbitrary downloads.
    with input_path.open("rb") as handle:
        payload = pickle.load(handle)
    if not isinstance(payload, dict):
        raise TypeError(f"Expected dict in {input_path}, got {type(payload).__name__}")

    vectors: dict[str, np.ndarray] = {}
    original_keys: dict[str, str] = {}
    for raw_key, raw_value in payload.items():
        slug = STATE_ALIASES.get(str(raw_key).strip().lower())
        if slug is None:
            continue
        if slug in vectors:
            raise ValueError(f"Duplicate centroid for {slug!r}")
        vectors[slug] = as_vector(raw_value)
        original_keys[slug] = str(raw_key)

    missing = set(REQUIRED_STATES) - set(vectors)
    if missing:
        raise ValueError(f"Missing required states: {sorted(missing)}; found {sorted(payload)}")
    dimensions = {vector.shape for vector in vectors.values()}
    if len(dimensions) != 1:
        raise ValueError(f"Centroid shapes differ: {sorted(dimensions)}")
    if any(np.linalg.norm(vector) == 0 for vector in vectors.values()):
        raise ValueError("Centroid vectors must be non-zero")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(output_path, **{state: vectors[state] for state in REQUIRED_STATES})
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_file": input_path.name,
        "source_sha256": sha256(input_path),
        "output_file": output_path.name,
        "output_sha256": sha256(output_path),
        "states": list(REQUIRED_STATES),
        "original_state_keys": original_keys,
        "embedding_dimension": int(next(iter(vectors.values())).size),
        "vector_norms": {state: float(np.linalg.norm(vectors[state])) for state in REQUIRED_STATES},
        "privacy": "three aggregate training-donor centroids; no cell- or donor-level records",
    }
    output_path.with_suffix(".json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path, help="Trusted Geneformer centroid pickle")
    parser.add_argument("--output", required=True, type=Path, help="Destination .npz path")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    print(json.dumps(export_centroids(args.input, args.output), indent=2))


if __name__ == "__main__":
    main()
