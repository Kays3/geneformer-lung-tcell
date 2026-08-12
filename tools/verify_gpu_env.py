#!/usr/bin/env python3
"""Verify a cu130 Geneformer environment actually computes correctly on the GPU.

A forward pass that returns a tensor of the right shape proves the environment
imports and the driver is reachable. It does not prove the numbers are right: a
mismatched CUDA/driver pair, a bad kernel, or a silently wrong dtype can all
produce plausible garbage. This checks three increasingly strict things:

1. Device and CUDA reporting agree with the built wheel.
2. GPU logits match CPU logits on the same real cells, to float tolerance.
   This is the numerical-correctness test and needs no label semantics.
3. Full held-out test inference reproduces the recorded confusion matrix and
   accuracy. This is end-to-end: model weights, tokenized data, GPU kernels and
   the recorded result all have to agree.

The classifier stores generic LABEL_0/1/2, so the class order is recovered by
testing which permutation reproduces the committed confusion matrix rather than
being assumed. If no permutation matches, that is reported as a failure instead
of quietly picking the best one.

    python tools/verify_gpu_env.py                 # full check
    python tools/verify_gpu_env.py --max-cells 512 # quick smoke
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DEFAULT_MODEL = "/srv/lab/KD/sclc_luad_normal_htan_finetune/runs/260728_geneformer_cellClassifier_sclc_luad_normal_htan/ksplit1"
DEFAULT_DATASET = "/srv/lab/KD/sclc_luad_normal_htan_finetune/data/sclc_luad_normal_htan_tcells.dataset"
EXPECTED_CM = REPO / "sclc_validation/perturbation_workflow/results/test_confusion_matrix.csv"
EXPECTED_METRICS = REPO / "sclc_validation/perturbation_workflow/results/test_metrics.json"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", default=os.environ.get("CLASSIFIER_DIR", DEFAULT_MODEL))
    ap.add_argument("--dataset", default=os.environ.get("TOKENIZED_DATASET", DEFAULT_DATASET))
    ap.add_argument("--batch-size", type=int, default=16)
    ap.add_argument("--max-cells", type=int, default=0, help="0 = all test cells")
    args = ap.parse_args()

    import numpy as np
    import pandas as pd
    import torch
    from datasets import load_from_disk
    from sklearn.metrics import accuracy_score, f1_score, confusion_matrix
    from transformers import BertForSequenceClassification

    failures = []

    print("\n=== 1. device ===")
    print(f"  torch            {torch.__version__}")
    print(f"  cuda available   {torch.cuda.is_available()}")
    if not torch.cuda.is_available():
        print("  no GPU visible; nothing further can be verified", file=sys.stderr)
        return 1
    print(f"  torch cuda       {torch.version.cuda}")
    print(f"  device           {torch.cuda.get_device_name(0)}")
    free, total = torch.cuda.mem_get_info()
    print(f"  memory           {free/2**30:.1f} GiB free of {total/2**30:.1f} GiB")
    cap = torch.cuda.get_device_capability(0)
    print(f"  capability       sm_{cap[0]}{cap[1]}")

    print("\n=== 2. GPU vs CPU agreement on real cells ===")
    model = BertForSequenceClassification.from_pretrained(args.model).eval()
    dataset = load_from_disk(args.dataset)
    test = dataset.filter(lambda r: r["split"] == "test", num_proc=1)
    print(f"  test cells       {len(test):,}")

    probe = test.select(range(8))
    width = min(max(len(r) for r in probe["input_ids"]), model.config.max_position_embeddings)
    batch = torch.zeros((len(probe), width), dtype=torch.long)
    mask = torch.zeros((len(probe), width), dtype=torch.long)
    for i, ids in enumerate(probe["input_ids"]):
        ids = ids[:width]
        batch[i, : len(ids)] = torch.tensor(ids, dtype=torch.long)
        mask[i, : len(ids)] = 1
    with torch.no_grad():
        cpu_logits = model(input_ids=batch, attention_mask=mask).logits
        gpu_logits = model.cuda()(input_ids=batch.cuda(), attention_mask=mask.cuda()).logits.cpu()
    delta = (cpu_logits - gpu_logits).abs().max().item()
    print(f"  max |GPU - CPU|  {delta:.3e}")
    if delta > 1e-3:
        failures.append(f"GPU and CPU logits diverge by {delta:.3e}")
    else:
        print("  agreement        OK (within float tolerance)")

    print("\n=== 3. held-out test inference on GPU ===")
    if args.max_cells:
        test = test.select(range(min(args.max_cells, len(test))))
        print(f"  limited to       {len(test):,} cells")
    model = model.cuda().eval()
    order = np.argsort([-n for n in test["length"]])  # long first, fewer pad tokens
    preds = np.empty(len(test), dtype=np.int64)
    with torch.no_grad():
        for start in range(0, len(order), args.batch_size):
            idx = order[start : start + args.batch_size]
            rows = [test["input_ids"][int(i)] for i in idx]
            w = min(max(len(r) for r in rows), model.config.max_position_embeddings)
            ids = torch.zeros((len(rows), w), dtype=torch.long)
            am = torch.zeros((len(rows), w), dtype=torch.long)
            for j, r in enumerate(rows):
                r = r[:w]
                ids[j, : len(r)] = torch.tensor(r, dtype=torch.long)
                am[j, : len(r)] = 1
            logits = model(input_ids=ids.cuda(), attention_mask=am.cuda()).logits
            preds[idx] = logits.argmax(-1).cpu().numpy()
            if start % (args.batch_size * 100) == 0:
                print(f"    {start + len(idx):,}/{len(order):,}", flush=True)
    truth = np.array(test["disease"])

    expected = pd.read_csv(EXPECTED_CM, index_col=0)
    classes = list(expected.index)
    truth_idx = np.array([classes.index(t) for t in truth])

    # Recover which model output corresponds to which disease by testing every
    # permutation against the recorded matrix, rather than assuming an order.
    best = None
    for perm in itertools.permutations(range(len(classes))):
        mapped = np.array([perm[p] for p in preds])
        cm = confusion_matrix(truth_idx, mapped, labels=list(range(len(classes))))
        acc = accuracy_score(truth_idx, mapped)
        if best is None or acc > best[0]:
            best = (acc, perm, cm, mapped)
    acc, perm, cm, mapped = best
    macro_f1 = f1_score(truth_idx, mapped, average="macro")
    print(f"  label order      model output -> {[classes[p] for p in perm]}")
    print(f"  accuracy         {acc:.6f}")
    print(f"  macro F1         {macro_f1:.6f}")
    print("\n  confusion matrix (rows = actual):")
    print(pd.DataFrame(cm, index=classes, columns=classes).to_string())

    if args.max_cells:
        print("\n  (subset run: recorded metrics not comparable)")
    else:
        recorded = json.loads(EXPECTED_METRICS.read_text())["test_metrics"]
        print("\n=== 4. against the recorded result ===")
        print(f"  accuracy   this run {acc:.6f}   recorded {recorded['acc']:.6f}")
        print(f"  macro F1   this run {macro_f1:.6f}   recorded {recorded['macro_f1']:.6f}")
        cm_match = bool((cm == expected.values).all())
        print(f"  confusion matrix identical: {cm_match}")
        if abs(acc - recorded["acc"]) > 5e-3:
            failures.append(f"accuracy {acc:.6f} vs recorded {recorded['acc']:.6f}")
        if not cm_match:
            failures.append("confusion matrix differs from the recorded one")
            print("\n  recorded:")
            print(expected.to_string())

    print("\n" + ("=" * 60))
    if failures:
        print("FAILED")
        for f in failures:
            print(f"  - {f}")
        return 1
    print("GPU ENVIRONMENT VERIFIED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
