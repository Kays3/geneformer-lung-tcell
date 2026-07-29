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
process ("Cannot re-initialize CUDA in forked subprocess"). Fixed by forcing
the `spawn` multiprocessing start method before any CUDA-touching import.

Check a live or completed run with:

```bash
./check_status.sh
```

Override the compute-output directory with `TARGETED_PANEL_RUN_DIR` when the
run is stored somewhere other than the default `~/workspace/KD/` location.
