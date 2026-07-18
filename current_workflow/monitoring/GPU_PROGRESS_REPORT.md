# Held-out all-gene deletion status

![Animated overall cell progress](progress_animation.gif)

> **15-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary plus the compact sketch below. The sketch is an orientation aid only;
> it is not a measured ligand-receptor or pathology result.

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

| Snapshot | Snapshot | Snapshot |
| --- | --- | --- |
| **2026-07-18T12:39:23+09:00**<br/>`████████░░░░` 64.28%<br/>Cells 2,172 · GPU 96% · 84 C · 88.8 W<br/>Shards 88 | **2026-07-18T11:39:20+09:00**<br/>`███████░░░░░` 61.35%<br/>Cells 2,073 · GPU 96% · 81 C · 90.1 W<br/>Shards 84 | **2026-07-18T10:39:17+09:00**<br/>`███████░░░░░` 58.42%<br/>Cells 1,974 · GPU 96% · 85 C · 89.0 W<br/>Shards 80 |
| **2026-07-18T09:39:14+09:00**<br/>`███████░░░░░` 54.72%<br/>Cells 1,849 · GPU 96% · 78 C · 89.6 W<br/>Shards 74 | **2026-07-18T08:39:11+09:00**<br/>`██████░░░░░░` 50.81%<br/>Cells 1,717 · GPU 91% · 80 C · 90.1 W<br/>Shards 69 | **2026-07-18T07:39:07+09:00**<br/>`██████░░░░░░` 46.67%<br/>Cells 1,577 · GPU 0% · 62 C · 14.9 W<br/>Shards 63 |

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
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, and `cell_interaction_diagram.svg`
- Cadence: 15 minutes
