# Two-machine workflow

This project runs across a laptop (authoring, review, figures) and two GB10 DGX Spark
nodes, `thinkstation2` and `thinkstation1` (GPU compute). Large data and run artifacts
live outside Git under `KD/` — see [`RELATIONSHIP_TO_KD.md`](RELATIONSHIP_TO_KD.md).

## Sync

```bash
./tools/sync.sh
```

Rebases onto `origin/main`, pushes, then fast-forwards `thinkstation2` directly. Run it
**before starting work on a machine and after finishing**. `--check` reports divergence
without changing anything; `--no-ts2` skips the workstation when it is unreachable.

The laptop can push straight to `thinkstation2` because that repo sets
`receive.denyCurrentBranch=updateInstead`, which updates its working tree in place — no
GitHub round-trip for a one-line script change. Setup, already applied:

```bash
ssh ts2 'cd ~/workspace/geneformer-lung-tcell && git config receive.denyCurrentBranch updateInstead'
git remote add ts2 ts2:workspace/geneformer-lung-tcell
```

Make rebase the default so a forgotten pull cannot create a merge bubble:

```bash
git config pull.rebase true
```

## Lessons from the fork this replaced

On 2026-08-07 the laptop and `thinkstation2` had each committed on top of the same
parent, forking history. Untangling it needed a manual rebase. Nothing was lost, but the
following practices would have prevented it entirely.

### 1. Run analyses from the repo checkout, not a scratch directory

The denoised spatial validation was first run by `scp`-ing the script to
`~/scratch_denoised_spatial/` on the workstation. It worked, but the script that
produced the results was not the script under version control — if the two had drifted,
the committed results would have been unreproducible.

Prefer:

```bash
./tools/sync.sh                                   # laptop: publish the script
ssh ts2 'cd ~/workspace/geneformer-lung-tcell && \
  GSE263196_RAW_DIR=... python sclc_validation/spatial_validation/spatial_validation_denoised_programs.py'
```

Provenance is then automatic: the commit hash on the workstation *is* the code that ran.

### 2. Keep machine-specific paths out of tracked files

Commit `16c061e` had to hand-edit `/home/petadimensionlab/...` to
`/home/thinkstation2/...` across 8 files. That class of change recurs on every host move.

The newer scripts already do this correctly — `SCLC_PERTURBATION_ROOT`,
`GSE263196_RAW_DIR` and `DENOISED_CANDIDATES` are environment overrides with defaults.
Extend that pattern rather than hardcoding: put per-host values in an untracked
`config/paths.env`, source it, and keep absolute paths out of the repository.

### 3. Preflight the inputs before a long run

The GSE263196 matrices had never been extracted — only the 244 MB
`GSE263196_RAW.tar` was on disk, so the spatial pipeline could not have run at all. That
was discovered by inspection, not by an error message.

A long GPU job should assert its inputs exist and are the expected shape before it starts,
so a missing file fails in seconds instead of after the sample-loading phase.

### 4. Pull results back into the repo, not into a home directory

Results were retrieved with an ad-hoc `rsync` into
`sclc_validation/spatial_validation/results_denoised_programs/`. That destination was the
right one, but the command was typed by hand. Wrapping the retrieval in the run script —
or committing on the workstation and syncing — removes the chance of results landing
somewhere untracked.

## Recovery

`tools/sync.sh` never force-pushes and refuses to run against a dirty tree, so it cannot
discard work. Before any history rewrite, tag the pre-rewrite state:

```bash
git branch -f backup/pre-sync-$(git rev-parse --short HEAD) HEAD
```

`thinkstation2` currently holds `backup/pre-sync-9804944` from the fork described above.
It can be deleted once the rebased history is confirmed good.
