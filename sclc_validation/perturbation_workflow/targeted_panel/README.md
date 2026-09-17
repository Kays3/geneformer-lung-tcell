# Targeted 50-gene panel perturbation

Scoped-down replacement for the full all-gene delete + overexpress sweep
(projected ~7-8 days; see `../METHODS.md`). Instead of every gene in every
held-out cell, this screen tests 50 genes individually, one gene at a time,
delete and overexpress, across all three source states.

## Target gene set (`target_gene_panel.json`)

- 21 genes: the audit's pre-registered signature panel.
- 29 genes: the top drivers by |shift| from the prior LUAD/LUSC/normal
  all-gene screen's `top_goal_shift_genes.csv`, deduplicated across the 6
  directional comparisons and with known ambient-RNA/contamination markers
  excluded (`SFTPC, SFTPB, NAPSA, MUC1, PIGR, FBLN1, DCN, ACKR1`, plus any
  `KRT*` keratin), per `evaluation/biology/README.md`'s own warning about
  the unfiltered rankings.

Rebuild with `build_target_gene_list.py` (reads `top_goal_shift_genes.csv`).

## Implementation note

Geneformer's `InSilicoPerturber.genes_to_perturb` treats a gene list as one
combined co-perturbation event, not an individual sweep. `run_targeted_panel.py`
therefore calls it once per gene (`genes_to_perturb=[single_ensembl_id]`),
50 x 3 sources x 2 types = 300 calls total, reusing the training-reference
and held-out-test datasets and the training-donor state embeddings already
computed for the (paused) all-gene sweep.

For `delete`, cells that don't detect the target gene are filtered out
automatically (cheap). For `overexpress` that filter does not apply -- every
held-out cell in the source is processed for every gene, since a gene can be
induced from an undetected state.

A specific-gene-list run also routes through a different internal code path
(`isp_perturb_set*`) that forks worker processes for `Dataset.map()` even at
`num_proc=1`; this crashes if CUDA is already initialized in the parent
process ("Cannot re-initialize CUDA in forked subprocess"). Forcing the
`spawn` multiprocessing start method before any CUDA-touching import does
**not** fix this on its own -- that only changes stdlib `multiprocessing`'s
default, but this fork comes from the separate `multiprocess` (dill-based)
package that `datasets.map()` actually uses, which has its own default and
ignores the stdlib setting. The real fix (2026-09-17, found via the bf16
precision-replication task): `NPROC` must be 1 for any GPU run, same as
`run_t4_overexpression.py` already does for the identical failure mode.

Check a live or completed run with:

```bash
./check_status.sh
```

Override the compute-output directory with `TARGETED_PANEL_RUN_DIR` when the
run is stored somewhere other than the default `~/workspace/KD/` location.
**Note (2026-09-17): this script does not currently read `TARGETED_PANEL_RUN_DIR`
and has never written under `~/workspace/KD/...` in its current form** --
that mismatch between this README/`check_status.sh` and the actual code
predates the bf16 task and is unrelated to it; flagged here, not fixed.
Separately, as of the bf16 precision-replication task, a run invoked with
`--dtype`/`--run-tag` (or with neither, which defaults `--run-tag` to
`fp32`) writes under
`sclc_validation/bf16_bench/runs/<run-tag>/targeted_panel/` instead --
check there for any run made after 2026-09-17.
