"""Two registered/planned checks that should have run BEFORE Phase 6 (run late, 2026-09-25).

1. NO-OP GATE (registration s.3, stop condition s.8): two independent forward passes
   over one donor's ISP input, no token edits, scored with Geneformer's own
   quant_cos_sims against that donor's goals.pkl, exactly as run_isp scores. Every
   per-cell shift must be exactly 0.0, for every state.
2. WRAPPER EQUIVALENCE (plan, Verification): the Phase 6 stack (dtype_cast +
   model_cache + inproc_map) must reproduce a NATIVE InSilicoPerturber call
   (dtype_cast only, which the bf16 model needs) byte for byte, both operations.

Usage: --mode noop | native   (native writes the pickles; comparison is done by sha256 outside)
"""
import argparse
import json
import os
import pickle
import sys

STATES = {"state_key": "origin", "start_state": "tumor_primary", "goal_state": "normal_adjacent",
          "alt_states": ["global_normal_adjacent"]}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--mode", choices=("noop", "native"), required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--donor", required=True)
    p.add_argument("--gene", required=True)
    p.add_argument("--bf16-bench", required=True)
    p.add_argument("--scripts", required=True)
    p.add_argument("--out", required=True)
    a = p.parse_args()
    sys.path[:0] = [a.bf16_bench, a.scripts]
    from dtype_cast import install_dtype_cast
    install_dtype_cast("bf16")                       # native path: NO model_cache, NO inproc_map
    import torch
    from datasets import load_from_disk
    from geneformer import InSilicoPerturber
    from geneformer import perturber_utils as pu
    from geneformer.emb_extractor import get_embs
    m = json.load(open(a.manifest))[a.donor]
    goals = pickle.load(open(m["goals"], "rb"))
    os.makedirs(a.out, exist_ok=True)
    isp_kw = dict(genes_to_perturb=[a.gene], combos=0, anchor_gene=None, model_type="CellClassifier",
                  num_classes=2, emb_mode="cls", filter_data=None, cell_states_to_model=STATES,
                  state_embs_dict=goals, max_ncells=None, emb_layer=0, forward_batch_size=128, nproc=1,
                  model_version="V2", clear_mem_ncells=1000)
    if a.mode == "native":
        for op in ("delete", "overexpress"):
            od = os.path.join(a.out, op, a.gene); os.makedirs(od, exist_ok=True)
            InSilicoPerturber(perturb_type=op, **isp_kw).perturb_data(
                model_directory=m["model_dir"], input_data_file=m["isp_input"], output_directory=od,
                output_prefix=a.donor.replace("/", "_"))
        print("NATIVE_DONE"); return
    isp = InSilicoPerturber(perturb_type="delete", **isp_kw)   # setup only; perturb_data never called
    model = pu.load_model("CellClassifier", 2, m["model_dir"], mode="eval")
    layer = pu.quant_layers(model) + 0
    ds = load_from_disk(m["isp_input"])
    ds = ds.filter(lambda x: isp.tokens_to_perturb[0] in x["input_ids"])
    acc = {}
    B = 128
    for i in range(0, len(ds), B):
        mb = ds.select(range(i, min(i + B, len(ds))))
        e1 = get_embs(model, mb, "cls", layer, isp.pad_token_id, B, token_gene_dict=isp.token_gene_dict,
                      summary_stat=None, silent=True)
        e2 = get_embs(model, mb, "cls", layer, isp.pad_token_id, B, token_gene_dict=isp.token_gene_dict,
                      summary_stat=None, silent=True)
        for s, v in pu.quant_cos_sims(e2, e1, STATES, goals, emb_mode="cell").items():
            acc.setdefault(s, []).extend(torch.as_tensor(v).float().flatten().tolist())
    res = {s: {"n": len(v), "max_abs": max(abs(x) for x in v), "n_nonzero": sum(x != 0 for x in v)}
           for s, v in acc.items()}
    res["PASS"] = all(r["max_abs"] == 0.0 for r in res.values())
    json.dump(res, open(os.path.join(a.out, "noop_gate.json"), "w"), indent=1)
    print(json.dumps(res))


if __name__ == "__main__":
    main()
