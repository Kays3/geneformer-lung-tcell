# Held-out all-gene deletion status

![Animated overall cell progress](progress_animation.gif)

> **15-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary, and the gallery below shows the recent single-cell snapshots. The
> thumbnails are orientation aids only; they are not measured ligand-receptor
> or pathology results.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Current snapshot

**What changed since the prior report:** The run advanced by 52 cells and 2 shards, lifting completion from 69.58% to 71.12%. LUSC remained complete; LUAD remained complete; NORMAL moved from 380 to 432 cells; GPU utilization fell from 96% to 91%.

| Metric | Value |
| --- | --- |
| Generated | 2026-07-18T15:16:02+09:00 |
| Run status | RUNNING |
| Overall cell progress | 2,403 / 3,379 (71.12%) |
| GPU | NVIDIA GB10 |
| GPU utilization | 91% |
| GPU temperature | 79 C |
| GPU power | 83.2 W |
| Perturbation GPU memory | 2,301 MiB |
| System memory used | 39.3 GiB |
| System memory available | 80.4 GiB |
| Swap used | 0.0 GiB |

### Progress by source

| Source | Cells | Shards | Raw files | Marker deletions |
| --- | --- | --- | --- | --- |
| LUSC | 560 / 560 (100.00%) | 23 / 23 | 1,120 | 348,313 / 348,313 |
| LUAD | 1,411 / 1,411 (100.00%) | 57 / 57 | 2,822 | 1,150,097 / 1,150,097 |
| NORMAL | 432 / 1,408 (30.68%) | 17 / 57 | 868 | 439,324 / 1,439,366 |

## Final statistical comparisons

**0 / 6 comparisons are complete. Final aggregation waits for the deletion screen, currently 97 / 137 shards.**

| Comparison | State | Result rows | Updated | Output |
| --- | --- | --- | --- | --- |
| LUSC → NORMAL | WAITING FOR PERTURBATION | — | — | `heldout_allgene_lusc_to_normal.csv` |
| LUSC → LUAD | WAITING FOR PERTURBATION | — | — | `heldout_allgene_lusc_to_luad.csv` |
| LUAD → NORMAL | WAITING FOR PERTURBATION | — | — | `heldout_allgene_luad_to_normal.csv` |
| LUAD → LUSC | WAITING FOR PERTURBATION | — | — | `heldout_allgene_luad_to_lusc.csv` |
| NORMAL → LUAD | WAITING FOR PERTURBATION | — | — | `heldout_allgene_normal_to_luad.csv` |
| NORMAL → LUSC | WAITING FOR PERTURBATION | — | — | `heldout_allgene_normal_to_lusc.csv` |

Result-row counts confirm artifact generation only; they do not establish
biological significance. Gene rankings should be interpreted only after all
six comparisons complete and coverage, FDR, and donor-consistency checks pass.

## Snapshot gallery

Each thumbnail is a single-cell diagram that encodes overall progress and the
LUAD, LUSC, and NORMAL pathology balance at that snapshot.

<table><thead><tr><th align="center">Single-cell snapshot</th><th align="center">Single-cell snapshot</th><th align="center">Single-cell snapshot</th></tr></thead><tbody><tr><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T151602+0900.svg" alt="Snapshot 2026-07-18T15:16:02+09:00" width="260"/><br/><sub>2026-07-18T15:16:02+09:00 · 71.12% · 2,403 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T143930+0900.svg" alt="Snapshot 2026-07-18T14:39:30+09:00" width="260"/><br/><sub>2026-07-18T14:39:30+09:00 · 69.58% · 2,351 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T133927+0900.svg" alt="Snapshot 2026-07-18T13:39:27+09:00" width="260"/><br/><sub>2026-07-18T13:39:27+09:00 · 66.82% · 2,258 cells</sub></td></tr><tr><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T123923+0900.svg" alt="Snapshot 2026-07-18T12:39:23+09:00" width="260"/><br/><sub>2026-07-18T12:39:23+09:00 · 64.28% · 2,172 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T113920+0900.svg" alt="Snapshot 2026-07-18T11:39:20+09:00" width="260"/><br/><sub>2026-07-18T11:39:20+09:00 · 61.35% · 2,073 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T103917+0900.svg" alt="Snapshot 2026-07-18T10:39:17+09:00" width="260"/><br/><sub>2026-07-18T10:39:17+09:00 · 58.42% · 1,974 cells</sub></td></tr></tbody></table>

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

| Timestamp | Cells | Progress | GPU util | Temp | Power | Shards |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-07-18T15:16:02+09:00 | 2,403 | 71.12% | 91% | 79 C | 83.2 W | 97 |
| 2026-07-18T14:39:30+09:00 | 2,351 | 69.58% | 96% | 77 C | 85.5 W | 95 |
| 2026-07-18T13:39:27+09:00 | 2,258 | 66.82% | 96% | 76 C | 84.1 W | 91 |
| 2026-07-18T12:39:23+09:00 | 2,172 | 64.28% | 96% | 84 C | 88.8 W | 88 |
| 2026-07-18T11:39:20+09:00 | 2,073 | 61.35% | 96% | 81 C | 90.1 W | 84 |
| 2026-07-18T10:39:17+09:00 | 1,974 | 58.42% | 96% | 85 C | 89.0 W | 80 |
| 2026-07-18T09:39:14+09:00 | 1,849 | 54.72% | 96% | 78 C | 89.6 W | 74 |
| 2026-07-18T08:39:11+09:00 | 1,717 | 50.81% | 91% | 80 C | 90.1 W | 69 |

## Job notes

- Job entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Statistics source: `/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/stats` (override with `PERTURBATION_STATS_DIR`)
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, `cell_interaction_diagram.svg`, and `snapshot_gallery/*.svg`
- Cadence: 15 minutes
