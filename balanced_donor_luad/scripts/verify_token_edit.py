"""Unregistered additional verification (run late, 2026-09-25): does InSilicoPerturber edit the
RIGHT token and nothing else? Records exactly what the library feeds the model: every get_embs()
call inside isp_perturb_set_special comes as a pair (original minibatch, perturbed minibatch).
Checks per cell:
  delete:      perturbed == original with every occurrence of the target token removed
  overexpress: perturbed == [original[0], target] + (original[1:] without target)
and that the target is present in the original. Any other difference fails."""
import argparse, json, os, pickle, sys


def main():
    p = argparse.ArgumentParser()
    for k in ("manifest", "donor", "gene", "token_dict", "bf16_bench", "scripts", "out"):
        p.add_argument("--" + k.replace("_", "-"), required=True)
    p.add_argument("--expect-token", type=int, default=None, help="NEGATIVE CONTROL ONLY: expect a different token")
    a = p.parse_args()
    sys.path[:0] = [a.bf16_bench, a.scripts]
    from dtype_cast import install_dtype_cast
    install_dtype_cast("bf16")
    import geneformer.in_silico_perturber as isp_mod
    from geneformer import InSilicoPerturber
    calls = []
    real = isp_mod.get_embs
    def recording_get_embs(model, data, *args, **kw):
        calls.append([list(x) for x in data["input_ids"]])
        return real(model, data, *args, **kw)
    isp_mod.get_embs = recording_get_embs
    m = json.load(open(a.manifest))[a.donor]
    goals = pickle.load(open(m["goals"], "rb"))
    target = pickle.load(open(a.token_dict, "rb"))[a.gene]
    states = {"state_key": "origin", "start_state": "tumor_primary", "goal_state": "normal_adjacent",
              "alt_states": ["global_normal_adjacent"]}
    edited = target
    if a.expect_token is not None:
        target = a.expect_token
    report = {"gene": a.gene, "target_token": int(target), "donor": a.donor, "negative_control": a.expect_token is not None}
    for op in ("delete", "overexpress"):
        calls.clear()
        od = os.path.join(a.out, op); os.makedirs(od, exist_ok=True)
        InSilicoPerturber(perturb_type=op, genes_to_perturb=[a.gene], combos=0, anchor_gene=None,
                          model_type="CellClassifier", num_classes=2, emb_mode="cls", filter_data=None,
                          cell_states_to_model=states, state_embs_dict=goals, max_ncells=None, emb_layer=0,
                          forward_batch_size=128, nproc=1, model_version="V2", clear_mem_ncells=1000
                          ).perturb_data(model_directory=m["model_dir"], input_data_file=m["isp_input"],
                                         output_directory=od, output_prefix="verify")
        assert len(calls) % 2 == 0, "get_embs calls not in (original, perturbed) pairs"
        n = bad = 0; first_bad = None; n_occ = []
        for orig_b, pert_b in zip(calls[0::2], calls[1::2]):
            assert len(orig_b) == len(pert_b)
            for o, q in zip(orig_b, pert_b):
                n += 1; n_occ.append(o.count(target))
                if op == "delete":
                    expected = [t for t in o if t != target]
                else:
                    expected = [o[0], target] + [t for t in o[1:] if t != target]
                ok = target in o and q == expected
                if not ok:
                    bad += 1
                    if first_bad is None:
                        first_bad = {"orig_head": o[:8], "pert_head": q[:8], "len_o": len(o), "len_q": len(q)}
        report[op] = {"cells_checked": n, "cells_wrong": bad, "first_wrong": first_bad,
                      "target_occurrences_per_cell": sorted(set(n_occ)), "PASS": n > 0 and bad == 0}
    report["PASS"] = report["delete"]["PASS"] and report["overexpress"]["PASS"]
    json.dump(report, open(os.path.join(a.out, "token_edit_verification.json"), "w"), indent=1)
    print(json.dumps(report, indent=1))


if __name__ == "__main__":
    main()
