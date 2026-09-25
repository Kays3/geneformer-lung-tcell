"""Synthetic Phase 6 output with PLANTED answers, in the exact on-disk format run_isp.py writes:
<root>/<op>/<gene>/<donor>.complete.json and
<root>/<op>/<gene>/in_silico_<op>_<donor>_cell_embs_dict_[<token>]_raw.pickle
holding {state: {(token, "cell_emb"): [per-cell shifts]}}.

No real Phase 6 output is read anywhere in the tests.
"""
import json
import os
import pickle
from collections import defaultdict

import numpy as np

STATES = ("tumor_primary", "normal_adjacent", "global_normal_adjacent")
N_FOLDS, PER_FOLD = 5, 8
DONORS = [f"{'Leader_Merad_2021' if i % 2 == 0 else 'Kim_Lee_2020'}_D{i:02d}" for i in range(N_FOLDS * PER_FOLD)]
FOLD = {d: i // PER_FOLD for i, d in enumerate(DONORS)}
CONTROLS = [f"CTRL{i:02d}" for i in range(20)]

# gene -> (panel, delete effect per fold (list of 5), overexpress effect, stratum, expected status)
PLANT = {
    "GA_REP":   ("A", [0.03] * 5, -0.03, "S0", "REPLICATED"),
    "GA_AMB":   ("A", [0.03] * 5, -0.03, "S0", "REPLICATED_AMBIENT"),
    "GA_OPEN":  ("A", [0.0] * 5, 0.0, "S0", "OPEN"),
    "GA_REV":   ("A", [-0.03] * 5, 0.03, "S0", "REVERSED"),
    "GA_LOWD":  ("A", None, None, None, "NOT_RUN"),               # ineligible before GPU: no stratum
    "GA_NOCTL": ("A", None, None, "SBAD", "NOT_ESTIMABLE_CONTROLS"),
    "GB_TOW":   ("B", [0.03] * 5, -0.03, "S0", "T_CELL_SIGNAL_TOWARD"),
    "GB_AWAY":  ("B", [-0.03] * 5, 0.03, "S0", "T_CELL_SIGNAL_AWAY"),
    "GB_INC":   ("B", [0.03] * 5, 0.03, "S0", "DOSE_INCOHERENT"),
    "GB_DEL":   ("B", [0.03] * 5, 0.0, "S0", "DELETION_ONLY"),
    "GB_FOLD":  ("B", [0.03, 0.03, 0.03, -0.001, -0.001], -0.03, "S0", "OPEN"),   # S3 downgrade
}
AMBIENT = {"GA_AMB"}


def token_of(gene):
    return 10000 + sum(ord(c) for c in gene)


def write_call(root, op, gene, donor, shifts_by_state, n):
    od = os.path.join(root, op, gene)
    os.makedirs(od, exist_ok=True)
    tok = token_of(gene)
    if n == 0:
        json.dump({"gene": gene, "op": op, "donor": donor, "status": "no_token_cells", "n_token_cells": 0},
                  open(os.path.join(od, f"{donor}.complete.json"), "w"))
        return
    d = {s: defaultdict(list, {(tok, "cell_emb"): list(map(float, v))}) for s, v in shifts_by_state.items()}
    pickle.dump(d, open(os.path.join(od, f"in_silico_{op}_{donor}_cell_embs_dict_[{tok}]_raw.pickle"), "wb"))
    json.dump({"gene": gene, "op": op, "donor": donor, "status": "done", "n_token_cells": n},
              open(os.path.join(od, f"{donor}.complete.json"), "w"))


def build(root, seed=7, n_cells=30, ctrl_missing=None):
    """Write the synthetic tree. ctrl_missing: {donor: k} makes k controls non-estimable in that donor."""
    rng = np.random.default_rng(seed)
    ctrl_missing = ctrl_missing or {}
    for op in ("delete", "overexpress"):
        for gene in CONTROLS:
            for i, d in enumerate(DONORS):
                miss = CONTROLS.index(gene) < ctrl_missing.get(d, 0)
                n = 3 if miss else n_cells                 # 3 < MIN_CELLS: present but not estimable
                sh = {s: rng.normal(0.0, 0.002, n) for s in STATES}
                write_call(root, op, gene, d, sh, n)
        for gene, (panel, dels, ovx, stratum, _) in PLANT.items():
            if dels is None or stratum == "SBAD":
                continue                                   # never perturbed in Phase 6
            for d in DONORS:
                eff = dels[FOLD[d]] if op == "delete" else ovx
                sh = {s: rng.normal(eff, 0.002, n_cells) for s in STATES}
                write_call(root, op, gene, d, sh, n_cells)


def design():
    pa = [{"ensembl_id": g, "symbol": g} for g, v in PLANT.items() if v[0] == "A"]
    pb = [{"ensembl_id": g, "symbol": g} for g, v in PLANT.items() if v[0] == "B"]
    return dict(
        panel_a=pa, panel_b=pb,
        panel_b_not_run=[{"ensembl_id": "GB_TRAC", "symbol": "TRAC", "reason": "NOT_RUN: not in V2 token dictionary"}],
        strata={"S0": {"status": "eligible", "controls": CONTROLS, "members": []},
                "SBAD": {"status": "not_estimable_control_stratum", "controls": [], "members": []}},
        gene_stratum={g: v[3] for g, v in PLANT.items() if v[3]},
        ambient_flagged={g: g in AMBIENT for g, v in PLANT.items() if v[0] == "A"},
        ambient_label={g: ("AMBIENT_FLAGGED" if g in AMBIENT else "NOT_FLAGGED") for g, v in PLANT.items() if v[0] == "A"},
        donors=DONORS, donor_fold=FOLD, donor_study={d: d.rsplit("_", 1)[0] for d in DONORS},
        donor_ba={d: 0.6 + 0.3 * (i / len(DONORS)) for i, d in enumerate(DONORS)},
        expected_estimable={**{g: len(DONORS) for g, v in PLANT.items() if v[1] is not None and v[3] == "S0"},
                            "GA_LOWD": 4},
        not_run_genes=set(),
    )
