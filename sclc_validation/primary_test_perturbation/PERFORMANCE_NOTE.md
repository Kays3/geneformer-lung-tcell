# Performance Note

## DGX Spark benchmark outcome

Benchmarks varying Geneformer's `forward_batch_size` and dataset `nproc` did
not materially change the speed of the all-gene perturbation calculations on
the DGX Spark GB10 machines.

The limiting factor is the inherent GPU memory and execution setup of this
class of DGX Spark systems, rather than only the nominal forward batch size or
the number of dataset worker processes. GPU utilization is high during model
forward passes, with intermittent idle periods caused by the per-perturbation
execution path and host-side orchestration.

## Operational conclusion

- Keep the validated `forward_batch_size=16` and `nproc=4` settings for this
  analysis unless a future hardware/software stack changes the benchmark.
- Do not assume that increasing batch size or dataset workers will improve
  throughput on DGX Spark GB10 nodes.
- Optimize future runs primarily through independent shard scheduling across
  additional GPUs, resumable execution, local NVMe placement, and staged
  targeted screens rather than aggressive per-GPU batch tuning.
- Preserve the same numerical settings for the primary held-out analysis to
  maintain comparability and reproducibility.

This note records the observed benchmark result for the current Geneformer V2
SCLC/LUAD/normal workflow; it is not a general claim about other GPU systems.
