#!/usr/bin/env python3
"""Verify migrated Geneformer assets, statistical outputs, and runtime."""

from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
from pathlib import Path


EXPECTED_GENEFORMER_COMMIT = "f45a6c7"
CRITICAL_HASHES = {
    "Geneformer-V2-104M/model.safetensors": "fff5cba29ddd8792991fa77b4872246fbe548a178cebda3775cdc72b67780e7f",
    "KD/tcell_luad_lusc_normal_luscmax7000_finetune/runs/260717_geneformer_cellClassifier_tcell_luad_lusc_normal_luscmax7000/ksplit1/model.safetensors": "039509ddf56121e5320ca24079aeaf05d6aeb57505e3d1e4d368b7744a190c0e",
}
OPTIONAL_ATLAS = (
    "KD/data/nsclc/nsclc_integrated.h5ad",
    "141db65b76b1e34f895131e36c74cd829db05fc037f8cd2f422c2960a5a266cd",
)
EXPECTED_STATS_ROWS = {
    "heldout_allgene_lusc_to_luad.csv": 11242,
    "heldout_allgene_lusc_to_normal.csv": 11242,
    "heldout_allgene_luad_to_lusc.csv": 13458,
    "heldout_allgene_luad_to_normal.csv": 13458,
    "heldout_allgene_normal_to_luad.csv": 14923,
    "heldout_allgene_normal_to_lusc.csv": 14923,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def check_hash(root: Path, relative: str, expected: str, *, announce: bool = True) -> None:
    path = root / relative
    if not path.is_file():
        raise AssertionError(f"Missing file: {path}")
    actual = sha256(path)
    if actual != expected:
        raise AssertionError(f"Checksum mismatch for {relative}: {actual}")
    if announce:
        print(f"OK hash: {relative}")


def check_manifest(root: Path, manifest: Path) -> None:
    checked = 0
    for line_number, line in enumerate(manifest.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        try:
            expected, relative = line.split(maxsplit=1)
        except ValueError as error:
            raise AssertionError(f"Invalid manifest line {line_number}") from error
        relative = relative.lstrip("* ")
        check_hash(root, relative, expected, announce=False)
        checked += 1
    print(f"OK manifest: {checked} files")


def check_stats(root: Path) -> None:
    stats = root / (
        "KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/stats"
    )
    for filename, expected_rows in EXPECTED_STATS_ROWS.items():
        path = stats / filename
        if not path.is_file():
            raise AssertionError(f"Missing statistics table: {path}")
        with path.open(newline="") as handle:
            reader = csv.DictReader(handle)
            rows = sum(1 for _ in reader)
            columns = set(reader.fieldnames or [])
        required = {"Gene_name", "Shift_to_goal_end", "Goal_end_FDR", "N_Detections"}
        if not required.issubset(columns):
            raise AssertionError(f"Missing columns in {filename}: {required - columns}")
        if rows != expected_rows:
            raise AssertionError(
                f"Row-count mismatch for {filename}: expected {expected_rows}, found {rows}"
            )
        print(f"OK stats: {filename} ({rows:,} rows)")


def check_git(root: Path) -> None:
    commit = subprocess.run(
        ["git", "-C", str(root), "rev-parse", "--short=7", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    ).stdout.strip()
    if commit != EXPECTED_GENEFORMER_COMMIT:
        raise AssertionError(
            f"Geneformer commit mismatch: expected {EXPECTED_GENEFORMER_COMMIT}, found {commit}"
        )
    print(f"OK Geneformer commit: {commit}")


def check_monitor(root: Path) -> None:
    required = [
        "current_workflow/monitoring/GPU_PROGRESS_REPORT.md",
        "current_workflow/perturbation_statistics/perturbation_statistics_report.html",
        "current_workflow/perturbation_statistics/perturbation_statistics.ipynb",
        "current_workflow/perturbation_statistics/evaluation/README.md",
    ]
    for relative in required:
        path = root / relative
        if not path.is_file() or path.stat().st_size == 0:
            raise AssertionError(f"Missing monitor artifact: {path}")
        print(f"OK monitor artifact: {relative}")


def check_runtime(root: Path) -> None:
    python = root / ".venv/bin/python"
    if not python.is_file():
        raise AssertionError(f"Missing target environment: {python}")
    code = """
import sys
sys.path.insert(0, '.')
import datasets, geneformer, torch, transformers
print('torch', torch.__version__)
print('torch_cuda', torch.version.cuda)
print('cuda_available', torch.cuda.is_available())
print('gpu', torch.cuda.get_device_name(0) if torch.cuda.is_available() else None)
print('transformers', transformers.__version__)
print('datasets', datasets.__version__)
print('geneformer_import', 'OK')
"""
    subprocess.run([str(python), "-c", code], cwd=root, check=True)
    print("OK runtime")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--geneformer-root", type=Path, required=True)
    parser.add_argument("--monitor-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--runtime", action="store_true")
    args = parser.parse_args()

    geneformer_root = args.geneformer_root.resolve()
    monitor_root = args.monitor_root.resolve()
    check_git(geneformer_root)
    for relative, expected in CRITICAL_HASHES.items():
        check_hash(geneformer_root, relative, expected)
    atlas_path = geneformer_root / OPTIONAL_ATLAS[0]
    if atlas_path.exists():
        check_hash(geneformer_root, *OPTIONAL_ATLAS)
    else:
        print("SKIP optional atlas hash: atlas not transferred")
    check_stats(geneformer_root)
    check_monitor(monitor_root)
    if args.manifest:
        check_manifest(geneformer_root, args.manifest.resolve())
    if args.runtime:
        check_runtime(geneformer_root)
    print("Migration verification passed")


if __name__ == "__main__":
    main()
