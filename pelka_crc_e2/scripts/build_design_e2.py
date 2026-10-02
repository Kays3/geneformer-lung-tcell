"""E2: build design.json for analyse.py (unchanged) from COMMITTED PRE-GPU inputs only.

Copy of balanced_donor_luad/scripts/build_design.py. Changes: no Panel A and no ambient table (Panel A is
not part of E2; ambient labels enter only Panel A statuses), so panel_a, ambient_flagged and ambient_label
are empty; the not-run set is DERIVED (panel genes in a stratum with fewer than 20 controls) instead of the
lung constant (BTG1); the panel B list is the balanced-donor file, read from the lung directory.

Reads (relative to pelka_crc_e2/ unless noted, all committed before any ISP call):
  ../balanced_donor_luad/registration/panel_B.json          panel, registered order
  controls/pre_gpu_eligibility.csv                         estimable donors, eligibility, stratum
  controls/pre_gpu_strata.json                             strata -> status, 20 controls, members
  controls/pre_gpu_host_map.csv                            genes Phase 6 was told to run
  cohort/folds.json                                        donor -> held-out fold, donor -> study
  phase4_results/classifier_gate.json                      per-donor held-out balanced accuracy (3f)

No Phase 6 output is read. Every cross-check below refuses (exit 1) rather than guessing.
"""
import argparse
import csv
import hashlib
import json
import os
import sys

INPUTS = {
    "panel_b": "../balanced_donor_luad/registration/panel_B.json",
    "eligibility": "controls/pre_gpu_eligibility.csv",
    "strata": "controls/pre_gpu_strata.json",
    "host_map": "controls/pre_gpu_host_map.csv",
    "folds": "cohort/folds.json",
    "classifier_gate": "phase4_results/classifier_gate.json",
}
D_MIN = {"B": 11}
N_CONTROLS = 20


class DesignError(ValueError):
    pass


def check(cond, msg):
    if not cond:
        raise DesignError(msg)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def read_csv(path):
    with open(path, newline="") as f:
        return list(csv.DictReader(f))


def build(base):
    p = {k: os.path.join(base, v) for k, v in INPUTS.items()}
    for k, v in p.items():
        check(os.path.isfile(v), f"missing input {k}: {v}")
    pb = json.load(open(p["panel_b"]))
    strata = json.load(open(p["strata"]))
    folds = json.load(open(p["folds"]))
    gate = json.load(open(p["classifier_gate"]))
    elig = {r["ensembl_id"]: r for r in read_csv(p["eligibility"])}
    host_genes = {r["ensembl_id"] for r in read_csv(p["host_map"])}

    # ---- panels
    panel_a = []
    panel_b = [{"ensembl_id": g["ensembl_id"], "symbol": g["symbol"]} for g in pb["genes"] if g["perturbable"]]
    panel_b_not_run = [{"ensembl_id": g["ensembl_id"], "symbol": g["symbol"], "reason": g["status_if_not"]}
                       for g in pb["genes"] if not g["perturbable"]]
    check(len(panel_b) == pb["n_perturbable"], f"Panel B: {len(panel_b)} perturbable != {pb['n_perturbable']}")
    check(len(panel_b) + len(panel_b_not_run) == pb["n_listed"], "Panel B listed count mismatch")
    check(all(nr["reason"] for nr in panel_b_not_run), "a Panel B NOT_RUN gene has no reason")
    ids = [g["ensembl_id"] for g in panel_a + panel_b]
    check(len(ids) == len(set(ids)), "a gene appears twice across panels")

    # ---- eligibility table must cover exactly the perturbable panel genes
    check(set(elig) == set(ids), f"eligibility table genes differ from panels: "
                                 f"{sorted(set(elig) ^ set(ids))[:5]}")
    for P, genes in (("B", panel_b),):
        for g in genes:
            r = elig[g["ensembl_id"]]
            check(r["panel"] == P, f"{g['ensembl_id']}: panel {r['panel']} != {P}")
            check(int(r["d_min"]) == D_MIN[P], f"{g['ensembl_id']}: d_min {r['d_min']} != registered {D_MIN[P]}")
            is_elig = r["eligible"] == "True"
            check(is_elig == (int(r["estimable_donors"]) >= D_MIN[P]),
                  f"{g['ensembl_id']}: eligible flag disagrees with estimable_donors >= d_min")
            check(is_elig == bool(r["stratum"]), f"{g['ensembl_id']}: stratum present iff eligible violated")
    expected_estimable = {g: int(elig[g]["estimable_donors"]) for g in ids}
    gene_stratum = {g: elig[g]["stratum"] for g in ids if elig[g]["stratum"]}

    # ---- strata: membership and control counts
    member_list = [(m, s) for s, v in strata.items() for m in v["members"]]
    members = dict(member_list)
    check(len(members) == len(member_list), "strata members: a gene is listed in more than one stratum")
    check(members == gene_stratum, "strata members disagree with the eligibility table's stratum column")
    for s, v in strata.items():
        ok = v["status"] == "eligible"
        check(ok == (len(v.get("controls") or []) == N_CONTROLS), f"stratum {s}: status {v['status']} vs "
                                                                   f"{len(v.get('controls') or [])} controls")
        for m in v["members"]:
            check((elig[m]["controls_status"] or None) == v["status"], f"{m}: controls_status != stratum status")
        check(not set(v.get("controls") or []) & set(ids), f"stratum {s}: a panel gene is used as a control")
    no_ctrl = {m for s, v in strata.items() if v["status"] != "eligible" for m in v["members"]}
    not_run_genes = no_ctrl

    # ---- the genes ISP is told to run == eligible panel genes + all controls (host map includes control-less
    #      members, which the E2 driver excludes)
    controls = {c for v in strata.values() for c in (v.get("controls") or [])}
    check(host_genes == set(gene_stratum) | controls, "host map != eligible panel genes + controls")

    ambient_flagged, ambient_label = {}, {}

    # ---- donors: folds, study, held-out BA must agree on the same 19
    donors = sorted(folds["donor_test_fold"])
    check(set(folds["donor_study"]) == set(donors), "folds.json donor_study keys != donor_test_fold keys")
    check(set(gate["per_donor_balanced_accuracy"]) == set(donors), "classifier gate donors != fold donors")
    check(len(donors) == gate["n_donors"], f"{len(donors)} donors != gate n_donors {gate['n_donors']}")
    check(gate["PASS"] is True, "classifier gate did not pass")
    for f in folds["folds"]:
        test = [d for d, k in folds["donor_test_fold"].items() if k == f["fold"]]
        check(not set(test) & set(f["train"]), f"fold {f['fold']}: a held-out donor is in its training set")
    check(max(expected_estimable.values()) <= len(donors), "estimable donors exceeds cohort size")

    design = dict(
        panel_a=panel_a, panel_b=panel_b, panel_b_not_run=panel_b_not_run,
        strata={s: {"status": v["status"], "controls": list(v.get("controls") or []), "members": v["members"]}
                for s, v in strata.items()},
        gene_stratum=gene_stratum, ambient_flagged=ambient_flagged, ambient_label=ambient_label,
        donors=donors, donor_fold={d: int(folds["donor_test_fold"][d]) for d in donors},
        donor_study={d: folds["donor_study"][d] for d in donors},
        donor_ba={d: float(gate["per_donor_balanced_accuracy"][d]) for d in donors},
        expected_estimable=expected_estimable, not_run_genes=sorted(not_run_genes),
    )
    provenance = {k: {"path": INPUTS[k], "sha256": sha256(v)} for k, v in p.items()}
    return design, provenance


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--base", required=True, help="pelka_crc_e2/ directory of the committed branch")
    ap.add_argument("--out", required=True, help="design.json to write")
    a = ap.parse_args(argv)
    try:
        design, prov = build(a.base)
    except DesignError as e:
        print(f"REFUSED: {e}", file=sys.stderr)
        return 1
    with open(a.out, "w") as f:
        json.dump(design, f, indent=1, sort_keys=True)
    with open(os.path.splitext(a.out)[0] + "_inputs.json", "w") as f:
        json.dump(prov, f, indent=1, sort_keys=True)
    n_run = len(design["gene_stratum"]) - len(design["not_run_genes"])
    print(json.dumps({"panel_a": len(design["panel_a"]), "panel_b": len(design["panel_b"]),
                      "panel_b_not_run": len(design["panel_b_not_run"]), "donors": len(design["donors"]),
                      "panel_genes_run": n_run, "strata": len(design["strata"])}))
    return 0


if __name__ == "__main__":
    sys.exit(main())
