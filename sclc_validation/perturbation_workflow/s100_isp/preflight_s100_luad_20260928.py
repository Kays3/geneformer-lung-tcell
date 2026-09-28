"""CPU-only gate-3 completion for the S100 LUAD run (2026-09-28).

Run from sclc_validation/perturbation_workflow/targeted_panel/ (it imports
run_targeted_panel.py from the working directory), with S100_ISP_DIR set to
this directory and the same root environment variables as the GPU run
(GENEFORMER_ROOT, SCLC_PERTURBATION_ROOT, HTAN_FINETUNE_ROOT,
STATE_EMB_FILE_OVERRIDE, TARGET_GENES_FILE_OVERRIDE). No model is loaded;
no GPU work. Outputs go to bf16_bench/runs/s100_luad_preflight_20260928/.

First run, 2026-09-28 on thinkstation1: two of the six strata (clean_high,
low_a9_pbp) have zero candidates under the registered all-members matching
rule, so the primary test is not estimable as registered. See
retained_rows_spec_20260922.md.
"""
from __future__ import annotations

import hashlib
import json
import os
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.argv = ["run_targeted_panel.py", "--run-tag", "s100_luad_preflight_20260928", "--dtype", "bf16"]
sys.path.insert(0, ".")
import run_targeted_panel as rtp  # noqa: E402

sys.path.insert(0, os.environ["S100_ISP_DIR"])
import matched_controls as mc  # noqa: E402
from datasets import load_from_disk  # noqa: E402

PANEL_FILE = Path("../s100_isp/s100_gene_panel_20260922.json").resolve()
OUT = rtp.BF16_BENCH_ROOT / "runs" / "s100_luad_preflight_20260928"
OUT.mkdir(parents=True, exist_ok=True)
LOG = []


def say(msg=""):
    print(msg, flush=True)
    LOG.append(str(msg))


def sha(p: Path) -> str:
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


panel = json.loads(PANEL_FILE.read_text())
panel_genes = panel["genes"]
say(f"panel {PANEL_FILE} sha256={sha(PANEL_FILE)}")
say(f"runner sha256={sha(Path('run_targeted_panel.py'))}")
say(f"TEST_DATASET={rtp.TEST_DATASET}")

tok = rtp.gene_token_dict()
say(f"token dictionary {rtp.TOKEN_DICTIONARY_FILE} sha256={sha(rtp.TOKEN_DICTIONARY_FILE)} n={len(tok)}")
name_file = Path(rtp.TOKEN_DICTIONARY_FILE).with_name(
    Path(rtp.TOKEN_DICTIONARY_FILE).name.replace("token_dictionary", "gene_name_id_dict"))
name_to_id = pickle.load(open(name_file, "rb"))
say(f"gene name dictionary {name_file} sha256={sha(name_file)} n={len(name_to_id)}")
id_to_names = defaultdict(set)
for name, eid in name_to_id.items():
    id_to_names[eid].add(str(name))
s100_family = {eid for eid, names in id_to_names.items() if any(n.upper().startswith("S100") for n in names)}
panel_ids = {g["ensembl_id"] for g in panel_genes}
say(f"S100-family ensembl ids by symbol prefix: {len(s100_family)}; panel ids covered: "
    f"{len(panel_ids & s100_family)}/12; family members not on panel: "
    f"{sorted(n for e in s100_family - panel_ids for n in id_to_names[e])}")
assert panel_ids <= s100_family, "a panel gene is missing from the S100-family set"

# ---- 1. frozen LUAD gene stats (the ruled source: heldout_test filtered to LUAD)
stats_path = OUT / "luad_gene_stats.csv"
stats, stats_sha = mc.load_and_freeze_luad_gene_stats(rtp.TEST_DATASET, tok, stats_path)
say(f"luad_gene_stats.csv sha256={stats_sha} rows={len(stats)}")

# ---- 2. eligibility counts in one pass over the runner's own LUAD source dataset
luad = load_from_disk(str(rtp.source_dataset_path("luad")))
say(f"LUAD cells: runner source dataset {len(luad)}")
cells_by_tok = defaultdict(int)
donors_by_tok = defaultdict(set)
cells_by_tok_donor = defaultdict(lambda: defaultdict(int))
for ids, ind in zip(luad["input_ids"], luad["individual"]):
    for t in set(ids):
        cells_by_tok[t] += 1
        donors_by_tok[t].add(ind)
        cells_by_tok_donor[t][ind] += 1
elig_rows = []
eligible_ids = set()
for eid, t in tok.items():
    n = cells_by_tok.get(t, 0)
    d = len(donors_by_tok.get(t, ()))
    ok = n >= rtp.MIN_CELLS_ELIGIBLE and d >= rtp.MIN_DONORS_ELIGIBLE
    capped = sum(min(c, rtp.DONOR_CELL_CAP) for c in cells_by_tok_donor[t].values()) if t in cells_by_tok_donor else 0
    if ok:
        eligible_ids.add(eid)
    elig_rows.append({"ensembl_id": eid, "n_token_positive_cells": n, "n_donors": d,
                      "n_cells_after_cap": capped, "isp_eligible": ok})
elig = rtp.pd.DataFrame(elig_rows)
elig_path = OUT / "luad_isp_eligibility_all_genes.csv"
elig.to_csv(elig_path, index=False)
say(f"luad_isp_eligibility_all_genes.csv sha256={sha(elig_path)}; eligible genes {len(eligible_ids)} / {len(tok)}")

# stats-table positive counts must agree with the eligibility pass (two code paths, same cells)
merged = stats.merge(elig, on="ensembl_id")
mismatch = merged[merged["n_positive_cells"] != merged["n_token_positive_cells"]]
say(f"cross-check n_positive_cells (stats path) vs n_token_positive_cells (eligibility path): "
    f"{len(mismatch)} mismatches over {len(merged)} genes")

def counts_agree(info, mine) -> bool:
    """paired_eligible_dataset() reports n_cells AFTER the per-donor cap when
    eligible (uncapped per-donor counts in donor_cell_counts), and the
    uncapped token-positive count when not eligible."""
    if bool(info["eligible"]) != bool(mine["isp_eligible"]) or info.get("n_donors") != mine["n_donors"]:
        return False
    if info["eligible"]:
        return info["n_cells"] == mine["n_cells_after_cap"] and \
            sum(info["donor_cell_counts"].values()) == mine["n_token_positive_cells"]
    return info["n_cells"] == mine["n_token_positive_cells"]


# ---- 3. panel genes: one-pass counts vs the runner's own paired_eligible_dataset()
say("\npanel genes, LUAD: one-pass count vs paired_eligible_dataset()")
by_id = elig.set_index("ensembl_id")
panel_rows = []
for g in panel_genes:
    _, info = rtp.paired_eligible_dataset("luad", g)
    mine = by_id.loc[g["ensembl_id"]]
    agree = counts_agree(info, mine)
    say(f"  {g['gene']:>8s} {g['stratum']:>14s} eligible={info['eligible']!s:5s} cells={int(mine['n_token_positive_cells'])} "
        f"donors={info.get('n_donors')} capped={mine['n_cells_after_cap']} agree={agree}")
    assert agree, (g, info, dict(mine))
    panel_rows.append({**g, "eligible": bool(info["eligible"]), "n_cells": int(mine["n_token_positive_cells"]),
                       "n_donors": info.get("n_donors"), "n_cells_after_cap": int(mine["n_cells_after_cap"]),
                       "donor_cell_counts": info.get("donor_cell_counts"), "reason": info.get("reason")})

# ---- 4. the six strata
strata = mc.build_matched_control_table(panel_genes, stats, isp_eligible_ensembl_ids=eligible_ids,
                                        excluded_ensembl_ids=s100_family)
say("\nmatched-control strata:")
for name, r in strata.items():
    say(f"  {name:>14s} status={r['status']} candidates={r['n_candidates']} "
        f"drawn={len(r['controls']) if r['controls'] else 0} reason={r.get('reason')}")
strata_path = OUT / "matched_control_strata.json"
strata_path.write_text(json.dumps({"seed": mc.SEED, "min_common_controls": mc.MIN_COMMON_CONTROLS,
                                   "luad_gene_stats_sha256": stats_sha,
                                   "eligibility_sha256": sha(elig_path),
                                   "strata": strata}, indent=2, default=list) + "\n")
say(f"matched_control_strata.json sha256={sha(strata_path)}")

# ---- 5. every drawn control through the runner's own eligibility / paired-list path
control_ids = sorted({c for r in strata.values() if r["controls"] for c in r["controls"]})
control_strata = defaultdict(list)
for name, r in strata.items():
    for c in (r["controls"] or ()):
        control_strata[c].append(name)
say(f"\ndistinct controls drawn: {len(control_ids)}; appearing in >1 stratum: "
    f"{sum(1 for c in control_ids if len(control_strata[c]) > 1)}")
control_rows = []
for eid in control_ids:
    sym = sorted(id_to_names.get(eid, {eid}))[0]
    gene = {"gene": sym, "ensembl_id": eid, "stratum": ",".join(control_strata[eid]), "role": "matched_control"}
    _, info = rtp.paired_eligible_dataset("luad", gene)
    assert info["eligible"], (gene, info)
    assert counts_agree(info, by_id.loc[eid]), (gene, info)
    control_rows.append({**gene, "n_cells": int(by_id.loc[eid, "n_token_positive_cells"]), "n_donors": info["n_donors"],
                         "n_cells_after_cap": int(by_id.loc[eid, "n_cells_after_cap"])})
say(f"all {len(control_rows)} drawn controls eligible via paired_eligible_dataset(); counts agree")

# ---- 6. run panel: registered panel genes + controls
run_panel = dict(panel)
run_panel["genes"] = panel_genes + [{k: r[k] for k in ("gene", "ensembl_id", "stratum", "role")} for r in control_rows]
run_panel["note"] = ("S100 LUAD run panel, 2026-09-28: the 12 registered genes (s100_gene_panel_20260922.json, "
                     "unchanged) plus the matched controls drawn in matched_control_strata.json.")
run_panel_path = OUT / "s100_luad_run_panel_20260928.json"
run_panel_path.write_text(json.dumps(run_panel, indent=2) + "\n")
say(f"run panel {run_panel_path.name}: {len(run_panel['genes'])} genes, sha256={sha(run_panel_path)}")

manifest = rtp.pd.DataFrame(panel_rows + [{**r, "eligible": True} for r in control_rows])
manifest_path = OUT / "s100_luad_eligibility_manifest.csv"
manifest.to_csv(manifest_path, index=False)
say(f"eligibility manifest sha256={sha(manifest_path)}")

# ---- 7. GPU projection
runs = rtp.BF16_BENCH_ROOT / "runs" / "316m_bf16" / "targeted_panel" / "raw" / "delete" / "luad"
xs, ys = [], []
for m in sorted(runs.glob("*.complete.json")):
    j = json.loads(m.read_text())
    if j.get("dtype") != "bf16" or "elapsed_seconds" not in j:
        continue
    t = tok.get(j["ensembl_id"])
    n = cells_by_tok.get(t, 0)
    if n:
        xs.append(n)
        ys.append(j["elapsed_seconds"])
xs, ys = np.array(xs, float), np.array(ys, float)
b, a = np.polyfit(xs, ys, 1)
say(f"\n316M bf16 LUAD delete calibration: {len(xs)} genes, cells {int(xs.min())}-{int(xs.max())}; "
    f"elapsed = {a:.2f} s + {b:.4f} s/cell")
noop_s_per_cell = 18.74 / 73
say(f"cross-check: today's no-op 18.74 s / 73 cells ({noop_s_per_cell:.3f} s/cell incl. load, 2 passes)")
run_genes = [r for r in panel_rows if r["eligible"]] + control_rows
per_arm = np.array([a + b * r["n_cells_after_cap"] for r in run_genes])
total_s = 2 * per_arm.sum()
say(f"projection: {len(run_genes)} eligible genes x 2 operations (delete + overexpress), LUAD only, "
    f"{sum(r['n_cells_after_cap'] for r in run_genes)} capped cells/arm total -> "
    f"{total_s:.0f} s = {total_s / 3600:.2f} GPU-h (budget 2.98)")

(OUT / "preflight_log.txt").write_text("\n".join(LOG) + "\n")
