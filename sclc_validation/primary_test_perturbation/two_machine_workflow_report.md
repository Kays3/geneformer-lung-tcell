# Two-Machine Perturbation Workflow Report

**Snapshot time:** 2026-08-06 18:04 JST  
**Workflow:** SCLC/LUAD/normal HTAN T-cell held-out all-gene perturbation

## Executive Summary

The distributed workflow is running normally on two machines. The machines use the same prepared dataset and model assets, but split the perturbation types across GPUs:

- `thinkstation2` (`192.168.100.2`) runs **gene deletion** perturbations.
- `thinkstation1` (`192.168.100.1`) runs **gene overexpression** perturbations.

Both jobs are actively producing output files. Overall progress is approximately **733 of 752 shard markers complete (97.5%)**.

## Work Allocation

Each perturbation type processes the three source groups independently: `normal`, `sclc`, and `luad`. Each source is divided into 25-cell dataset shards.

| Machine | Perturbation | Normal | SCLC | LUAD | Current state |
|---|---|---:|---:|---:|---|
| `thinkstation2` | `delete` | 23/23 | 97/97 | 242/256 | LUAD shard 244 active; shard 243 was the latest completed marker |
| `thinkstation1` | `overexpress` | 23/23 | 97/97 | 251/256 | LUAD shard 252 active; shard 251 was the latest completed marker |

The shard counts are completion-marker counts. A shard being active may not yet have a corresponding `.complete.json` marker.

## Workflow Between the Machines

1. The orchestration script prepares and synchronizes the model and held-out shard assets to both machines.
2. `thinkstation2` launches `run_heldout_allgene.py` with `--perturb-types delete`.
3. `thinkstation1` launches the same runner with `--perturb-types overexpress` over SSH, using its local GPU.
4. Each process loads one source shard at a time, runs forward passes in batches, and writes raw perturbation outputs plus a completion marker under its local `raw/` directory.
5. A monitor records both process IDs, shard counts, and GPU readings in `progress.log`.
6. Once both jobs finish, the raw outputs can be synchronized or combined for downstream analysis.

## Compute State

- Both main Python processes have been running for approximately 4 days 23 hours.
- `thinkstation2` local process: PID `251301`, four multiprocessing workers.
- `thinkstation1` remote process: PID `249898`, four multiprocessing workers.
- `thinkstation2` GPU: utilization varies between idle at batch boundaries and approximately 96% during inference; latest direct sample was 693 MiB memory and 61°C.
- `thinkstation1` GPU: approximately 89% utilization, 2.33 GiB memory, and 81°C in the latest sample.

The temporary local GPU idle readings are expected between shards or batches; the continuously advancing logs and output files show that the job is not stalled.

## Network Use

The workflow uses the direct `enp1s0f0np0` link between `192.168.100.2` and `192.168.100.1` for SSH orchestration, monitoring, and asset synchronization. A two-second traffic sample showed very low active transfer:

- `thinkstation2`: approximately 0.001 MB/s received and near-zero transmitted.
- `thinkstation1`: approximately 0.001 MB/s received and 0.001 MB/s transmitted.

This indicates that the current workload is primarily GPU and local-storage bound. The network is not carrying a sustained data stream during inference.

## Logs and Monitoring

- Local progress monitor: `sclc_validation/primary_test_perturbation/progress.log`
- Local workflow log: `/home/thinkstation2/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/logs/2node/local_delete.log`
- Remote workflow log: `/home/thinkstation1/workspace/KD/sclc_luad_normal_htan_heldout_allgene_perturbation/logs/2node/remote_overexpress.log`

No fatal runtime errors were observed in the latest inspection. Resource-tracker semaphore cleanup warnings appear in the logs, but processing continues and output markers are advancing.

## Expected Completion

The remaining work is concentrated in the final LUAD shards on both machines. Based on the recent shard completion rate, several more hours may be required, with `thinkstation2` currently the slower side of the pair.
