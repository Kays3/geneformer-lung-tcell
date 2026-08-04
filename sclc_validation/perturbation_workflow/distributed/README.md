# Two-node all-gene execution

The SCLC all-gene screen is an independent, resumable workload rather than a
`torchrun`/DDP model. `InSilicoPerturber` loads its own model and processes
datasets locally, so the reliable multi-GPU design is to assign deletion to one
GPU and overexpression to the other. The direct 200G interface is used for
asset transfer and result synchronization.

On `thinkstation2` (the local node), after confirming SSH access to
`thinkstation1` over `192.168.100.1`:

```bash
cd ~/workspace/geneformer-lung-tcell
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh prepare
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh start
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh status
```

After both arms report all shards complete:

```bash
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh sync
./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh stats
```

The launcher defaults to the documented topology and can be overridden without
editing it:

```bash
REMOTE_USER=thinkstation1 REMOTE_IP=192.168.100.1 IFACE=enp1s0f0np0 \
  ./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh status
```

`prepare` transfers the fine-tune workspace and analysis inputs, then runs one
delete and one overexpression smoke test. `start` is resumable through the
existing completion markers. It does not launch NCCL or move tensors over the
network during inference.
