# Reproducible machine migration

This directory migrates the active NSCLC Geneformer experiment and its report
repository to another Linux machine. It deliberately rebuilds Python instead
of copying the existing virtual environment, whose compiled packages and
absolute symlinks are tied to the source machine.

## Pinned source state

- Monitor repository: current `main` branch of
  `https://github.com/Kays3/geneformer-nsclc.git`.
- Geneformer upstream: `https://huggingface.co/ctheodoris/Geneformer`, commit
  `f45a6c7`.
- Python: 3.12.
- Environment manager: `uv`, using the transferred `pyproject.toml` and
  `uv.lock`.
- Source hardware at inventory time: Ubuntu 24.04 ARM64, NVIDIA GB10, CUDA
  13.0, and PyTorch `2.12.0+cu130`.

The target may use a different CPU architecture or NVIDIA GPU, but its driver
and selected PyTorch wheel must be compatible. Never copy `.venv/`.

## Migration scopes

The default active-experiment scope is approximately 3.2 GB:

```text
Geneformer-V2-104M/
KD/tcell_luad_lusc_normal_10k_from_atlas/
KD/tcell_luad_lusc_normal_luscmax7000_finetune/
KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/
pyproject.toml
uv.lock
.python-version
```

Add `KD/data/nsclc/nsclc_integrated.h5ad`—approximately 13 GB—when the target
must reproduce cohort construction from the original atlas. The complete
105 GB source workspace contains unrelated experiments and is intentionally
outside this migration scope.

## 1. Configure the source

Copy the example without committing the resulting local file:

```bash
cd /home/petadimensionlab/workspace/geneformer-nsclc-monitor
cp migration/migration.env.example migration/migration.env
```

Edit `migration/migration.env` for the source and target machines. The file is
ignored by Git because it may contain a private SSH hostname.

## 2. Inventory and checksum the source

```bash
set -a
. migration/migration.env
set +a

migration/scripts/inventory_source.sh \
  "$SOURCE_GF" \
  "$MIGRATION_INVENTORY_DIR"
```

To include the original atlas in the checksum manifest:

```bash
migration/scripts/inventory_source.sh \
  "$SOURCE_GF" \
  "$MIGRATION_INVENTORY_DIR" \
  --include-atlas
```

The inventory records system/GPU information, Git state, installed Python
packages, sizes, and SHA-256 hashes. Review `git-state.txt`: the upstream
Geneformer checkout contains local and untracked work, so cloning upstream
alone is not a backup of this experiment.

## 3. Prepare clean repositories on the target

Run on the target machine:

```bash
mkdir -p /home/petadimensionlab/workspace
cd /home/petadimensionlab/workspace

git clone https://github.com/Kays3/geneformer-nsclc.git \
  geneformer-nsclc-monitor

git clone https://huggingface.co/ctheodoris/Geneformer
cd Geneformer
git checkout f45a6c7
```

Use the same absolute directory layout when possible. If it differs, update
the environment variables and cron command described below.

## 4. Transfer active assets

Run from the monitor repository on the source machine:

```bash
set -a
. migration/migration.env
set +a

migration/scripts/transfer_active_assets.sh \
  "$SOURCE_GF" \
  "$TARGET_HOST" \
  "$TARGET_GF"
```

Add `--include-atlas` for full-from-atlas reproduction. The transfer uses
`rsync --partial` and never uses `--delete`.

Copy the generated inventory separately so it can be checked on the target:

```bash
rsync -aH --partial --info=progress2 \
  "$MIGRATION_INVENTORY_DIR/" \
  "$TARGET_HOST:$TARGET_INVENTORY_DIR/"
```

## 5. Rebuild the environment on the target

Install Python 3.12 and `uv`, then run:

```bash
cd /home/petadimensionlab/workspace/geneformer-nsclc-monitor

migration/scripts/bootstrap_target.sh \
  /home/petadimensionlab/workspace/Geneformer
```

The bootstrap refuses a mismatched Geneformer source commit and executes
`uv sync --frozen`. If the target driver cannot run the CUDA 13.0 PyTorch
build, resolve that platform compatibility before executing perturbations.

## 6. Verify files, results, and runtime

Without a full checksum manifest:

```bash
cd /home/petadimensionlab/workspace/geneformer-nsclc-monitor

migration/scripts/verify_target.py \
  --geneformer-root /home/petadimensionlab/workspace/Geneformer \
  --monitor-root /home/petadimensionlab/workspace/geneformer-nsclc-monitor \
  --runtime
```

With the generated manifest:

```bash
migration/scripts/verify_target.py \
  --geneformer-root /home/petadimensionlab/workspace/Geneformer \
  --monitor-root /home/petadimensionlab/workspace/geneformer-nsclc-monitor \
  --manifest "$TARGET_INVENTORY_DIR/files.sha256" \
  --runtime
```

The verifier checks critical model hashes, the six expected statistical
tables and row counts, report artifacts, the optional full manifest, local
Geneformer import, package versions, and CUDA visibility.

After verification, the existing small perturbation smoke test can be run:

```bash
cd /home/petadimensionlab/workspace/Geneformer

.venv/bin/python \
  KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/scripts/run_heldout_allgene.py \
  smoke-test
```

## 7. Restore reporting only after validation

The completed run does not need active monitoring. For a future resumed or new
run, use target-specific values rather than copying the old crontab blindly:

```cron
*/30 * * * * PUBLISH_TO_GIT=1 ANALYSIS_ROOT=/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation PYTHON_BIN=/home/petadimensionlab/workspace/Geneformer/.venv/bin/python /home/petadimensionlab/workspace/geneformer-nsclc-monitor/current_workflow/monitoring/refresh_live_report.sh >> /home/petadimensionlab/workspace/geneformer-nsclc-monitor/current_workflow/monitoring/report_generation.log 2>&1
```

The monitor scripts also accept `PERTURBATION_STATS_DIR`. Historical JSON and
Markdown files retain source-machine paths as provenance; they do not control
the relative-path perturbation runner.

## 8. Cutover checklist

- Source monitor repository is clean and pushed.
- Source inventory and checksum manifest are retained outside the source disk.
- Target critical hashes and six statistical table counts pass.
- Target Geneformer import and CUDA smoke checks pass.
- Perturbation smoke test completes.
- Technical HTML report and executed notebook open on the target.
- GitHub and Hugging Face authentication are configured independently.
- Old scheduled jobs are disabled before enabling target jobs.
- Source data is retained until at least one independent target backup exists.

Do not transfer SSH keys, GitHub credentials, Hugging Face tokens, or other
secrets with the experiment directories.
