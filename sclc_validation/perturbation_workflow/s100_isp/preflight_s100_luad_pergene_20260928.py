"""CPU-only rebuild of the S100 LUAD controls under per-gene matching
(human ruling 2026-09-28 ~13:00 JST; Amendment 8).

Run from sclc_validation/perturbation_workflow/targeted_panel/ with
S100_ISP_DIR set to the s100_isp directory and the same root environment
variables as the GPU run. Reuses the frozen LUAD gene stats and
eligibility tables from s100_luad_preflight_20260928/ (hashes checked, not
recomputed). No model is loaded; no GPU work.
"""
from __future__ import annotations

import hashlib
import json
import os
import pickle
import sys
from collections import defaultdict
from pathlib import Path

sys.argv = ["run_targeted_panel.py", "--run-tag", "s100_luad_preflight_pergene_20260928", "--dtype", "bf16"]
sys.path.insert(0, ".")
import run_targeted_panel as rtp  # noqa: E402

sys.path.insert(0, os.environ["S100_ISP_DIR"])
import matched_controls as mc  # noqa: E402

FROZEN = rtp.BF16_BENCH_ROOT / "runs" / "s100_luad_preflight_20260928"
OUT = rtp.BF16_BENCH_ROOT / "runs" / "s100_luad_preflight_pergene_20260928"
OUT.mkdir(parents=True, exist_ok=True)
PANEL_FILE = Path("../s100_isp/s100_gene_panel_20260922.json").resolve()
PRIMARY = ["S100A2", "S100B", "S100A13", "S100PBP", "S100A4", "S100A6", "S100A10", "S100A11"]
EXPECT = {"luad_gene_stats.csv": "cbad7c7bfde174c102f0c20fff0a9415e936574bd62994e76749e4acf4e35926",
          "luad_isp_eligibility_all_genes.csv": "552ee81da50b9dc5ab90962f575e2718e1c8e587e9c81d11b3e3ffdf1820fe53"}
LOAD_S, PER_CELL_S = 6.7, (40.1 - 6.7) / 300
LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def sha(p) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


say(f"runner sha256={sha('run_targeted_panel.py')}")
say(f"panel sha256={sha(PANEL_FILE)}")
for name, want in EXPECT.items():
    got = sha(FROZEN / name)
    say(f"{name} sha256={got} matches_frozen={got == want}")
    assert got == want, name

panel = json.loads(PANEL_FILE.read_text())
panel_genes = panel["genes"]
stats = rtp.pd.read_csv(FROZEN / "luad_gene_stats.csv")
elig = rtp.pd.read_csv(FROZEN / "luad_isp_eligibility_all_genes.csv").set_index("ensembl_id")
eligible_ids = set(elig.index[elig["isp_eligible"]])

name_file = Path(rtp.TOKEN_DICTIONARY_FILE).with_name(
    Path(rtp.TOKEN_DICTIONARY_FILE).name.replace("token_dictionary", "gene_name_id_dict"))
id_to_names = defaultdict(set)
for name, eid in pickle.load(open(name_file, "rb")).items():
    id_to_names[eid].add(str(name))
s100_family = {eid for eid, names in id_to_names.items() if any(n.upper().startswith("S100") for n in names)}
say(f"S100-family ids excluded: {len(s100_family)}")

per_gene = {g: (g,) for g in PRIMARY}
strata = mc.build_matched_control_table(panel_genes, stats, strata=per_gene,
                                        isp_eligible_ensembl_ids=eligible_ids,
                                        excluded_ensembl_ids=s100_family)
say("\nper-gene matched controls (seed 20260922 + gene symbol):")
for g in PRIMARY:
    r = strata[g]
    say(f"  {g:>8s} status={r['status']} candidates={r['n_candidates']} drawn={len(r['controls'] or ())}")
    assert r["status"] == "eligible", (g, r)
strata_path = OUT / "matched_control_sets_pergene.json"
strata_path.write_text(json.dumps({"rule": "per-gene (human ruling 2026-09-28)", "seed": mc.SEED,
                                   "min_common_controls": mc.MIN_COMMON_CONTROLS,
                                   "luad_gene_stats_sha256": EXPECT["luad_gene_stats.csv"],
                                   "eligibility_sha256": EXPECT["luad_isp_eligibility_all_genes.csv"],
                                   "sets": strata}, indent=2, default=list) + "\n")
say(f"matched_control_sets_pergene.json sha256={sha(strata_path)}")

used_by = defaultdict(list)
for g, r in strata.items():
    for c in r["controls"]:
        used_by[c].append(g)
control_ids = sorted(used_by)
say(f"\ndistinct controls: {len(control_ids)} (of {sum(len(r['controls']) for r in strata.values())} draws); "
    f"shared by >1 gene: {sum(1 for c in control_ids if len(used_by[c]) > 1)}")

control_rows = []
for eid in control_ids:
    sym = sorted(id_to_names.get(eid, {eid}))[0]
    gene = {"gene": sym, "ensembl_id": eid, "stratum": ",".join(used_by[eid]), "role": "matched_control"}
    _, info = rtp.paired_eligible_dataset("luad", gene)
    mine = elig.loc[eid]
    assert info["eligible"] and info["n_donors"] == mine["n_donors"] \
        and info["n_cells"] == mine["n_cells_after_cap"] \
        and sum(info["donor_cell_counts"].values()) == mine["n_token_positive_cells"], (gene, info)
    control_rows.append({**gene, "n_cells": int(mine["n_token_positive_cells"]), "n_donors": int(info["n_donors"]),
                         "n_cells_after_cap": int(mine["n_cells_after_cap"]), "eligible": True})
say(f"all {len(control_rows)} controls eligible via paired_eligible_dataset(); counts agree")

panel_rows = []
for g in panel_genes:
    _, info = rtp.paired_eligible_dataset("luad", g)
    mine = elig.loc[g["ensembl_id"]]
    panel_rows.append({**g, "eligible": bool(info["eligible"]), "n_cells": int(mine["n_token_positive_cells"]),
                       "n_donors": int(mine["n_donors"]), "n_cells_after_cap": int(mine["n_cells_after_cap"])})

run_panel = dict(panel)
run_panel["genes"] = panel_genes + [{k: r[k] for k in ("gene", "ensembl_id", "stratum", "role")} for r in control_rows]
run_panel["note"] = ("S100 LUAD run panel, per-gene matching (Amendment 8, 2026-09-28): the 12 registered genes "
                     "(s100_gene_panel_20260922.json, unchanged) plus the per-gene matched controls in "
                     "matched_control_sets_pergene.json. For controls, 'stratum' lists the primary gene(s) they serve.")
run_panel_path = OUT / "s100_luad_run_panel_pergene_20260928.json"
run_panel_path.write_text(json.dumps(run_panel, indent=2) + "\n")
say(f"run panel {run_panel_path.name}: {len(run_panel['genes'])} genes, sha256={sha(run_panel_path)}")

manifest = rtp.pd.DataFrame(panel_rows + control_rows)
manifest_path = OUT / "s100_luad_eligibility_manifest_pergene.csv"
manifest.to_csv(manifest_path, index=False)
say(f"eligibility manifest sha256={sha(manifest_path)}")

run_genes = [r for r in panel_rows if r["eligible"]] + control_rows
arm_s = [LOAD_S + PER_CELL_S * r["n_cells_after_cap"] for r in run_genes]
total = 2 * sum(arm_s)
say(f"\nprojection: {len(run_genes)} genes run ({sum(1 for r in panel_rows if r['eligible'])} panel + "
    f"{len(control_rows)} controls) x 2 operations, LUAD only; capped cells per arm total "
    f"{sum(r['n_cells_after_cap'] for r in run_genes)}; per arm = {LOAD_S} s + {PER_CELL_S:.4f} s/cell "
    f"-> {total:.0f} s = {total / 3600:.2f} GPU-h (hard stop 3.4)")
(OUT / "preflight_pergene_log.txt").write_text("\n".join(LOG) + "\n")
