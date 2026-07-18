# Held-out all-gene deletion status

![Animated overall cell progress](progress_animation.gif)

> **30-minute report job:** This report is generated from `latest_status.json`
> and `hourly_history.csv` on each run. Every snapshot includes a short delta
> summary. A single final diagram is appended for each disease source when its
> perturbation screen completes.

![Single-cell-inspired interaction sketch](cell_interaction_diagram.svg)

## Job run summaries

Newest refreshes are appended at the top and retained for the most recent
48 runs.

<!-- JOB_RUN_SUMMARIES_START -->
- **2026-07-18T18:00:01+09:00** — RUNNING; 2,646 / 3,379 cells (78.31%). The run advanced by 16 cells and 1 shards, lifting completion from 77.83% to 78.31%. LUSC remained complete; LUAD remained complete; NORMAL moved from 659 to 675 cells; GPU utilization rose from 5% to 75%.
- **2026-07-18T17:52:31+09:00** — RUNNING; 2,630 / 3,379 cells (77.83%). The run advanced by 5 cells and 0 shards, lifting completion from 77.69% to 77.83%. LUSC remained complete; LUAD remained complete; NORMAL moved from 654 to 659 cells; GPU utilization fell from 91% to 5%.
- **2026-07-18T17:48:04+09:00** — RUNNING; 2,625 / 3,379 cells (77.69%). The run advanced by 12 cells and 1 shards, lifting completion from 77.33% to 77.69%. LUSC remained complete; LUAD remained complete; NORMAL moved from 642 to 654 cells; GPU utilization rose from 4% to 91%.
- **2026-07-18T17:41:08+09:00** — RUNNING; 2,613 / 3,379 cells (77.33%). The run advanced by 210 cells and 8 shards, lifting completion from 71.12% to 77.33%. LUSC remained complete; LUAD remained complete; NORMAL moved from 432 to 642 cells; GPU utilization fell from 91% to 4%.
<!-- JOB_RUN_SUMMARIES_END -->

## Current snapshot

**What changed since the prior report:** The run advanced by 16 cells and 1 shards, lifting completion from 77.83% to 78.31%. LUSC remained complete; LUAD remained complete; NORMAL moved from 659 to 675 cells; GPU utilization rose from 5% to 75%.

| Metric | Value |
| --- | --- |
| Generated | 2026-07-18T18:00:01+09:00 |
| Run status | RUNNING |
| Overall cell progress | 2,646 / 3,379 (78.31%) |
| GPU | NVIDIA GB10 |
| GPU utilization | 75% |
| GPU temperature | 75 C |
| GPU power | 85.2 W |
| Perturbation GPU memory | 3,381 MiB |
| System memory used | 39.9 GiB |
| System memory available | 79.8 GiB |
| Swap used | 0.0 GiB |

### Progress by source

| Source | Cells | Shards | Raw files | Marker deletions |
| --- | --- | --- | --- | --- |
| LUSC | 560 / 560 (100.00%) | 23 / 23 | 1,120 | 348,313 / 348,313 |
| LUAD | 1,411 / 1,411 (100.00%) | 57 / 57 | 2,822 | 1,150,097 / 1,150,097 |
| NORMAL | 675 / 1,408 (47.94%) | 27 / 57 | 1,359 | 696,156 / 1,439,366 |

## Final statistical comparisons

**0 / 6 comparisons are complete. Final aggregation waits for the deletion screen, currently 107 / 137 shards.**

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

## Monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

The history table below shows the newest samples first.

| Timestamp | Cells | Progress | GPU util | Temp | Power | Shards |
| --- | --- | --- | --- | --- | --- | --- |
| 2026-07-18T18:00:01+09:00 | 2,646 | 78.31% | 75% | 75 C | 85.2 W | 107 |
| 2026-07-18T17:52:31+09:00 | 2,630 | 77.83% | 5% | 62 C | 15.7 W | 106 |
| 2026-07-18T17:48:04+09:00 | 2,625 | 77.69% | 91% | 78 C | 83.2 W | 106 |
| 2026-07-18T17:41:08+09:00 | 2,613 | 77.33% | 4% | 67 C | 16.5 W | 105 |
| 2026-07-18T15:16:02+09:00 | 2,403 | 71.12% | 91% | 79 C | 83.2 W | 97 |
| 2026-07-18T14:39:30+09:00 | 2,351 | 69.58% | 96% | 77 C | 85.5 W | 95 |
| 2026-07-18T13:39:27+09:00 | 2,258 | 66.82% | 96% | 76 C | 84.1 W | 91 |
| 2026-07-18T12:39:23+09:00 | 2,172 | 64.28% | 96% | 84 C | 88.8 W | 88 |

## Job notes

- Scheduled entrypoint: `current_workflow/monitoring/refresh_live_report.sh`
- Render entrypoint: `current_workflow/monitoring/generate_progress_report.py`
- Statistics source: `/home/petadimensionlab/workspace/Geneformer/KD/tcell_luad_lusc_normal_luscmax7000_heldout_allgene_perturbation/stats` (override with `PERTURBATION_STATS_DIR`)
- Output files: `GPU_PROGRESS_REPORT.md`, `progress_animation.gif`, `progress_animation.svg`, `cell_interaction_diagram.svg`, and `disease_completion/*.svg`
- Cadence: 30 minutes


## Disease completion diagrams

One final diagram is appended for each source after its perturbation screen completes. Cell totals are shown explicitly.

<table><tbody><tr><td align="center" valign="top"><img src="disease_completion/lusc_complete.svg" alt="LUSC perturbation complete" width="360"/></td><td align="center" valign="top"><img src="disease_completion/luad_complete.svg" alt="LUAD perturbation complete" width="360"/></td></tr></tbody></table>
