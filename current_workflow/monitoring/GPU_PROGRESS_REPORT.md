# Held-out all-gene deletion status

![Animated overall cell progress](progress_animation.gif)

> **15-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary, and the gallery below shows the recent single-cell snapshots. The
> thumbnails are orientation aids only; they are not measured ligand-receptor
> or pathology results.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Current snapshot

**What changed since the prior report:** The run advanced by 99 cells and 4 shards, lifting completion from 61.35% to 64.28%. LUSC remained complete; LUAD remained complete; NORMAL moved from 102 to 201 cells; GPU utilization stayed at 96%.

| Metric | Value |
| --- | --- |
| Generated | 2026-07-18T12:39:23+09:00 |
| Run status | RUNNING |
| Overall cell progress | 2,172 / 3,379 (64.28%) |
| GPU | NVIDIA GB10 |
| GPU utilization | 96% |
| GPU temperature | 84 C |
| GPU power | 88.8 W |
| Perturbation GPU memory | 3,151 MiB |
| System memory used | 39.1 GiB |
| System memory available | 80.6 GiB |
| Swap used | 0.0 GiB |

### Progress by source

| Source | Cells | Shards | Raw files | Marker deletions |
| --- | --- | --- | --- | --- |
| LUSC | 560 / 560 (100.00%) | 23 / 23 | 1,120 | 348,313 / 348,313 |
| LUAD | 1,411 / 1,411 (100.00%) | 57 / 57 | 2,822 | 1,150,097 / 1,150,097 |
| NORMAL | 201 / 1,408 (14.28%) | 8 / 57 | 404 | 204,266 / 1,439,366 |

## Snapshot gallery

Each thumbnail is a single-cell diagram that encodes overall progress and the
LUAD, LUSC, and NORMAL pathology balance at that snapshot.

<table><thead><tr><th align="center">Single-cell snapshot</th><th align="center">Single-cell snapshot</th><th align="center">Single-cell snapshot</th></tr></thead><tbody><tr><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T123923+0900.svg" alt="Snapshot 2026-07-18T12:39:23+09:00" width="260"/><br/><sub>2026-07-18T12:39:23+09:00 · 64.28% · 2,172 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T113920+0900.svg" alt="Snapshot 2026-07-18T11:39:20+09:00" width="260"/><br/><sub>2026-07-18T11:39:20+09:00 · 61.35% · 2,073 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T103917+0900.svg" alt="Snapshot 2026-07-18T10:39:17+09:00" width="260"/><br/><sub>2026-07-18T10:39:17+09:00 · 58.42% · 1,974 cells</sub></td></tr><tr><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T093914+0900.svg" alt="Snapshot 2026-07-18T09:39:14+09:00" width="260"/><br/><sub>2026-07-18T09:39:14+09:00 · 54.72% · 1,849 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T083911+0900.svg" alt="Snapshot 2026-07-18T08:39:11+09:00" width="260"/><br/><sub>2026-07-18T08:39:11+09:00 · 50.81% · 1,717 cells</sub></td><td align="center" valign="top"><img src="snapshot_gallery/snapshot_20260718T073907+0900.svg" alt="Snapshot 2026-07-18T07:39:07+09:00" width="260"/><br/><sub>2026-07-18T07:39:07+09:00 · 46.67% · 1,577 cells</sub></td></tr></tbody></table>

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

| Timestamp | Cells | Progress | GPU util | Temp | Power | Shards |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-07-18T12:39:23+09:00 | 2,172 | 64.28% | 96% | 84 C | 88.8 W | 88 |
| 2026-07-18T11:39:20+09:00 | 2,073 | 61.35% | 96% | 81 C | 90.1 W | 84 |
| 2026-07-18T10:39:17+09:00 | 1,974 | 58.42% | 96% | 85 C | 89.0 W | 80 |
| 2026-07-18T09:39:14+09:00 | 1,849 | 54.72% | 96% | 78 C | 89.6 W | 74 |
| 2026-07-18T08:39:11+09:00 | 1,717 | 50.81% | 91% | 80 C | 90.1 W | 69 |
| 2026-07-18T07:39:07+09:00 | 1,577 | 46.67% | 0% | 62 C | 14.9 W | 63 |
| 2026-07-18T06:39:04+09:00 | 1,436 | 42.50% | 0% | 63 C | 15.1 W | 58 |
| 2026-07-18T05:39:01+09:00 | 1,296 | 38.35% | 96% | 80 C | 90.7 W | 52 |

## Job notes

- Job entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, `cell_interaction_diagram.svg`, and `snapshot_gallery/*.svg`
- Cadence: 15 minutes
