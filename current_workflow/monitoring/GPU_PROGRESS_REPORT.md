# Held-out all-gene perturbation status

![Animated overall cell progress](progress_animation.gif)

## GPU monitoring history

![GPU utilization, temperature, power, and memory over time](gpu_statistics.png)

Values are point-in-time monitor samples; brief compute and idle phases may occur between observations.

Newest snapshots appear first. Existing entries are retained below; the machine-readable history is in `hourly_history.csv`.

<!-- NEWEST_SNAPSHOTS_BELOW -->

---

## Snapshot — 2026-07-17T21:43:58+09:00

Run status: **RUNNING**

Overall cell progress: **144 / 3,379 (4.26%)**

### GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 4% |
| GPU temperature | 63 C |
| GPU power | 16.0 W |
| Perturbation GPU memory | 701 MiB |
| System memory used | 41.5 GiB |
| System memory available | 78.2 GiB |
| Swap used | 0.0 GiB |

### Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 144 / 560 (25.71%) | 5 / 23 | 288 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.


---

## Snapshot — 2026-07-17T21:42:53+09:00

Run status: **RUNNING**

Overall cell progress: **140 / 3,379 (4.14%)**

### GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 93% |
| GPU temperature | 72 C |
| GPU power | 86.9 W |
| Perturbation GPU memory | 1391 MiB |
| System memory used | 38.8 GiB |
| System memory available | 80.9 GiB |
| Swap used | 0.0 GiB |

### Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 140 / 560 (25.00%) | 5 / 23 | 281 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.


---

## Snapshot — 2026-07-17T21:37:30+09:00

Run status: **RUNNING**

Overall cell progress: **126 / 3,379 (3.73%)**

### GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 4% |
| GPU temperature | 67 C |
| GPU power | 16.4 W |
| Perturbation GPU memory | 701 MiB |
| System memory used | 41.3 GiB |
| System memory available | 78.3 GiB |
| Swap used | 0.0 GiB |

### Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 126 / 560 (22.50%) | 5 / 23 | 252 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.


<details>
<summary>Earlier report content</summary>

# Held-out all-gene perturbation status

Last updated: **2026-07-17T21:22:30+09:00**

Run status: **RUNNING**

## GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 3% |
| GPU temperature | 68 C |
| GPU power | 21.4 W |
| Perturbation GPU memory | 701 MiB |
| System memory used | 40.5 GiB |
| System memory available | 79.2 GiB |
| Swap used | 0.0 GiB |

## Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 82 / 560 (14.64%) | 3 / 23 | 164 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.

## Recent hourly history

```csv
timestamp,run_active,gpu_utilization_percent,gpu_temperature_c,gpu_power_w,gpu_process_memory_mib,system_used_gib,system_available_gib,lusc_cells,luad_cells,normal_cells,completed_shards
2026-07-17T21:04:09+09:00,True,4.0,67.0,16.64,701.0,41.488990783691406,78.20438766479492,30,0,0,1
2026-07-17T21:06:10+09:00,True,93.0,78.0,86.44,1591.0,38.868473052978516,80.82490539550781,35,0,0,1
2026-07-17T21:22:30+09:00,True,3.0,68.0,21.39,701.0,40.47722625732422,79.21615219116211,82,0,0,3
```


---

## Snapshot — 2026-07-17T21:31:12+09:00

Run status: **RUNNING**

### GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 91% |
| GPU temperature | 75 C |
| GPU power | 67.8 W |
| Perturbation GPU memory | 1571 MiB |
| System memory used | 38.8 GiB |
| System memory available | 80.9 GiB |
| Swap used | 0.0 GiB |

### Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 106 / 560 (18.93%) | 4 / 23 | 213 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.


</details>
