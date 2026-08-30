#!/usr/bin/env python3
"""Run resumable T4 program-level overexpression on a Geneformer GPU host.

The manifest separates primary programs, the exhaustion nested-set titration,
and matched nulls. Use --dry-run on any machine; Geneformer/PyTorch/datasets are
imported only after the plan and external assets pass preflight.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import multiprocessing
import os
import pickle
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

multiprocessing.set_start_method("spawn", force=True)

HERE = Path(__file__).resolve().parent
HOME = Path.home()
MANIFEST = HERE / "t4_program_manifest.json"
ALLGENE_ROOT = Path(
    os.environ.get(
        "SCLC_PERTURBATION_ROOT",
        HOME / "workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation",
    )
)
FINETUNE_ROOT = Path(
    os.environ.get(
        "HTAN_FINETUNE_ROOT",
        HOME / "workspace/KD/sclc_luad_normal_htan_finetune",
    )
)
RUN_DIR = Path(
    os.environ.get(
        "T4_RUN_DIR",
        HOME / "workspace/KD/sclc_luad_normal_htan_program_overexpression",
    )
)
GENEFORMER_ROOT = Path(
    os.environ.get(
        "GENEFORMER_ROOT",
        HOME / "workspace/geneformer-uv-starter/geneformer-workspace/Geneformer",
    )
)
TOKEN_DICT = Path(
    os.environ.get(
        "GENEFORMER_TOKEN_DICT",
        GENEFORMER_ROOT / "geneformer/token_dictionary_gc104M.pkl",
    )
)
MODEL_PATH_FILE = FINETUNE_ROOT / "runs/MODEL_SCLC_LUAD_NORMAL_HTAN_PATH.txt"
TEST_DATASET = ALLGENE_ROOT / "data/heldout_test.dataset"
STATE_EMBEDDINGS = ALLGENE_ROOT / "state_embeddings/training_donor_disease_centroids.pkl"

STATE_NAMES = {
    "normal": "normal",
    "sclc": "small cell lung carcinoma",
    "luad": "lung adenocarcinoma",
}
STATE_SLUGS = {value: key for key, value in STATE_NAMES.items()}
PHASES = ("program", "nested", "null")
FORWARD_BATCH_SIZE = 128
NPROC = 4


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest(path: Path) -> dict:
    payload = json.loads(path.read_text())
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported T4 manifest schema: {payload.get('schema_version')}")
    runs = payload.get("runs", [])
    ids = [run.get("id") for run in runs]
    if len(ids) != len(set(ids)) or not all(ids):
        raise ValueError("T4 run IDs must be present and unique")
    for run in runs:
        if run.get("phase") not in PHASES:
            raise ValueError(f"Invalid phase in {run.get('id')}: {run.get('phase')}")
        genes = run.get("genes", [])
        if not genes or len({gene["ensembl_id"] for gene in genes}) != len(genes):
            raise ValueError(f"Missing or duplicate genes in {run.get('id')}")
    return payload


def completion_marker(run_dir: Path, run: dict, source: str) -> Path:
    return run_dir / "markers" / run["phase"] / run["id"] / f"{source}.complete.json"


def select_work(
    manifest: dict,
    run_dir: Path,
    phase: str,
    sources: list[str],
    items: list[str] | None = None,
    include_completed: bool = False,
    max_runs: int | None = None,
) -> list[tuple[dict, str]]:
    phases = set(PHASES if phase == "all" else [phase])
    item_filter = set(items or [])
    if item_filter:
        unknown = item_filter - {run["id"] for run in manifest["runs"]}
        if unknown:
            raise ValueError(f"Unknown T4 item IDs: {sorted(unknown)}")
    work = []
    for run in manifest["runs"]:
        if run["phase"] not in phases or (item_filter and run["id"] not in item_filter):
            continue
        for source in sources:
            if include_completed or not completion_marker(run_dir, run, source).exists():
                work.append((run, source))
    return work[:max_runs] if max_runs is not None else work


def summarize_work(work: list[tuple[dict, str]]) -> dict:
    by_phase = {phase: 0 for phase in PHASES}
    by_source = {source: 0 for source in STATE_NAMES}
    for run, source in work:
        by_phase[run["phase"]] += 1
        by_source[source] += 1
    return {
        "n_gpu_runs": len(work),
        "by_phase": by_phase,
        "by_source": by_source,
        "first_runs": [f"{run['id']}/{source}" for run, source in work[:12]],
    }


def preflight(manifest_path: Path, manifest: dict) -> tuple[Path, dict]:
    required = {
        "manifest": manifest_path,
        "Geneformer checkout": GENEFORMER_ROOT,
        "token dictionary": TOKEN_DICT,
        "model path file": MODEL_PATH_FILE,
        "held-out dataset": TEST_DATASET,
        "state embeddings": STATE_EMBEDDINGS,
    }
    missing = {name: str(path) for name, path in required.items() if not path.exists()}
    if missing:
        raise FileNotFoundError(f"T4 preflight missing assets: {json.dumps(missing, indent=2)}")
    model_dir = Path(MODEL_PATH_FILE.read_text().strip())
    if not model_dir.exists():
        raise FileNotFoundError(f"Model directory from {MODEL_PATH_FILE} does not exist: {model_dir}")
    with TOKEN_DICT.open("rb") as handle:
        token_dict = pickle.load(handle)
    requested = {
        gene["ensembl_id"]
        for run in manifest["runs"]
        for gene in run["genes"]
    }
    absent = sorted(requested - set(token_dict))
    if absent:
        raise ValueError(f"{len(absent)} manifest genes are absent from the token dictionary: {absent[:20]}")
    return model_dir, token_dict


def source_dataset_path(run_dir: Path, source: str, load_from_disk) -> Path:
    path = run_dir / "data/sources" / f"{source}.dataset"
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    dataset = load_from_disk(str(TEST_DATASET))
    disease = STATE_NAMES[source]
    subset = dataset.filter(lambda row, value=disease: row["disease"] == value, num_proc=1)
    if not len(subset):
        raise ValueError(f"No held-out cells for {source} ({disease})")
    subset.save_to_disk(str(path))
    return path


def canonical_states(source: str) -> dict:
    start = STATE_NAMES[source]
    others = [state for state in STATE_NAMES.values() if state != start]
    return {"state_key": "disease", "start_state": start, "goal_state": others[0], "alt_states": [others[1]]}


def run_one(
    run: dict,
    source: str,
    run_dir: Path,
    model_dir: Path,
    state_embeddings: dict,
    load_from_disk,
    InSilicoPerturber,
    forward_batch_size: int,
    nproc: int,
    force: bool,
) -> dict:
    marker = completion_marker(run_dir, run, source)
    if marker.exists() and not force:
        return json.loads(marker.read_text())
    raw_dir = run_dir / "raw" / run["phase"] / run["id"] / source
    raw_dir.mkdir(parents=True, exist_ok=True)
    prefix = f"t4_{run['id']}_{source}"
    if force:
        for path in raw_dir.glob(f"in_silico_overexpress_{prefix}_*_raw.pickle"):
            path.unlink()
    data_path = source_dataset_path(run_dir, source, load_from_disk)
    genes = [gene["ensembl_id"] for gene in run["genes"]]
    started = time.time()
    perturber = InSilicoPerturber(
        perturb_type="overexpress",
        genes_to_perturb=genes,
        combos=0,
        anchor_gene=None,
        model_type="CellClassifier",
        num_classes=3,
        emb_mode="cls",
        filter_data=None,
        cell_states_to_model=canonical_states(source),
        state_embs_dict=state_embeddings,
        max_ncells=None,
        emb_layer=0,
        forward_batch_size=forward_batch_size,
        nproc=nproc,
        model_version="V2",
        clear_mem_ncells=1000,
    )
    perturber.perturb_data(
        model_directory=str(model_dir),
        input_data_file=str(data_path),
        output_directory=str(raw_dir),
        output_prefix=prefix,
    )
    outputs = sorted(raw_dir.glob(f"in_silico_overexpress_{prefix}_*_raw.pickle"))
    if len(outputs) != 1:
        raise RuntimeError(f"Expected one raw pickle for {run['id']}/{source}, found {len(outputs)}")
    payload = {
        "completed_utc": utc_now(),
        "item_id": run["id"],
        "phase": run["phase"],
        "program": run["program"],
        "source": source,
        "genes": run["genes"],
        "elapsed_seconds": time.time() - started,
        "raw_file": str(outputs[0].relative_to(run_dir)),
        "raw_sha256": sha256(outputs[0]),
    }
    marker.parent.mkdir(parents=True, exist_ok=True)
    marker.write_text(json.dumps(payload, indent=2) + "\n")
    return payload


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    parser.add_argument("--run-dir", type=Path, default=RUN_DIR)
    parser.add_argument("--phase", choices=("all",) + PHASES, default="program")
    parser.add_argument("--source", choices=tuple(STATE_NAMES), action="append")
    parser.add_argument("--item", action="append", help="Restrict to a manifest item ID; repeatable")
    parser.add_argument("--max-runs", type=int)
    parser.add_argument("--forward-batch-size", type=int, default=FORWARD_BATCH_SIZE)
    parser.add_argument("--nproc", type=int, default=NPROC)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    manifest = load_manifest(args.manifest)
    sources = args.source or list(STATE_NAMES)
    work = select_work(
        manifest,
        args.run_dir,
        args.phase,
        sources,
        args.item,
        include_completed=args.force,
        max_runs=args.max_runs,
    )
    print(json.dumps(summarize_work(work), indent=2))
    if args.dry_run or not work:
        return

    model_dir, _ = preflight(args.manifest, manifest)
    sys.path.insert(0, str(GENEFORMER_ROOT))
    import torch
    from datasets import load_from_disk
    from geneformer import InSilicoPerturber

    if not torch.cuda.is_available():
        raise RuntimeError("T4 requires CUDA; use --dry-run to inspect the plan")
    with STATE_EMBEDDINGS.open("rb") as handle:
        state_embeddings = pickle.load(handle)
    args.run_dir.mkdir(parents=True, exist_ok=True)
    (args.run_dir / "run_config.json").write_text(
        json.dumps(
            {
                "started_utc": utc_now(),
                "manifest": str(args.manifest),
                "manifest_sha256": sha256(args.manifest),
                "model_directory": str(model_dir),
                "forward_batch_size": args.forward_batch_size,
                "nproc": args.nproc,
                "selection": summarize_work(work),
            },
            indent=2,
        )
        + "\n"
    )
    for index, (run, source) in enumerate(work, start=1):
        print(f"[{index}/{len(work)}] {run['phase']}/{run['id']}/{source}", flush=True)
        result = run_one(
            run,
            source,
            args.run_dir,
            model_dir,
            state_embeddings,
            load_from_disk,
            InSilicoPerturber,
            args.forward_batch_size,
            args.nproc,
            args.force,
        )
        print(f"  complete in {result['elapsed_seconds']:.1f}s", flush=True)


if __name__ == "__main__":
    main()
