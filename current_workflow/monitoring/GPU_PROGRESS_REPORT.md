# Held-out all-gene perturbation status

Last updated: **2026-07-17T21:06:10+09:00**

Run status: **RUNNING**

## GPU and memory

| Metric | Current value |
|---|---:|
| GPU | NVIDIA GB10 |
| GPU utilization | 93% |
| GPU temperature | 78 C |
| GPU power | 86.4 W |
| Perturbation GPU memory | 1591 MiB |
| System memory used | 38.9 GiB |
| System memory available | 80.8 GiB |
| Swap used | 0.0 GiB |

## Perturbation progress

| Source | Cells written | Shards complete | Raw files |
|---|---:|---:|---:|
| LUSC | 35 / 560 (6.25%) | 1 / 23 | 70 |
| LUAD | 0 / 1,411 (0.00%) | 0 / 57 | 0 |
| NORMAL | 0 / 1,408 (0.00%) | 0 / 57 | 0 |

Statistics are generated after all shards for a source state complete. Cell
counts include only cells with both raw checkpoint batches present.

## Recent hourly history

```csv
timestamp,run_active,gpu_utilization_percent,gpu_temperature_c,gpu_power_w,gpu_process_memory_mib,system_used_gib,system_available_gib,lusc_cells,luad_cells,normal_cells,completed_shards
2026-07-17T21:04:09+09:00,True,4.0,67.0,16.64,701.0,41.488990783691406,78.20438766479492,30,0,0,1
2026-07-17T21:06:10+09:00,True,93.0,78.0,86.44,1591.0,38.868473052978516,80.82490539550781,35,0,0,1
```
