# Two-node all-gene execution

The SCLC all-gene screen is an independent, resumable workload rather than a
`torchrun`/DDP model. `InSilicoPerturber` loads its own model and processes
datasets locally, so the reliable multi-GPU design is to assign deletion to one
GPU and overexpression to the other. The direct 200G interface is used for
asset transfer and result synchronization.

On `thinkstation2` (the local node), after confirming SSH access to
`thinkstation1` over the direct link (`<NODE1_IP>`):

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

The launcher's network settings are placeholders. Set them for your machines
without editing the script:

```bash
REMOTE_USER=<SSH_USER> REMOTE_IP=<NODE1_IP> IFACE=<IFACE> \
  ./sclc_validation/perturbation_workflow/distributed/run_2node_allgene.sh status
```

Fill in the placeholders: `<SSH_USER>` is the SSH user on the remote node, `<NODE1_IP>` its address on the
direct link, and `<IFACE>` the local interface of that link (`ip -br addr` lists them). Export them once per
shell or prefix each command as above; `monitor_2node_progress.sh` reads the same `REMOTE_USER` and `REMOTE_IP`.

`prepare` transfers the fine-tune workspace and analysis inputs, then runs one
delete and one overexpression smoke test. `start` is resumable through the
existing completion markers. It does not launch NCCL or move tensors over the
network during inference.
