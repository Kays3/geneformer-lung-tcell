"""Phase 3a: build the balanced paired-donor cohort from the verified LuCA file.

Selection is the registered one (PHASE0_SELECTION_RULE.md, Phase 1 survey):
LUAD, CD4/CD8 T cells, singlets, 10x assays, origin tumor_primary vs
normal_adjacent, donors qualifying (>= 100 cells) on BOTH sides.

The draw is deterministic and independent of file row order: for each
(donor, origin) the cell IDs are sorted, then permuted by a generator seeded
from sha256("<seed>|<donor>|<origin>"). A cell's rank in that permutation fixes
both cohorts, which are nested:
  - analysis cohort: rank < C_MIN (100), used for ISP statistics;
  - training pool:   rank < C_TRAIN (300), or all cells if a side has fewer.

Outputs (committed; small and determinative):
  cohort/cohort_cells.csv (training-pool cells; analysis cells are flagged),
  cohort/COHORT_DEFINITION.json
Outputs (not committed; regenerable):
  <slim_dir>/balanced_donor_luad_pool300.h5ad  raw integer counts for tokenisation

Balance assertions A1-A4 and F2-F4 are evaluated on the drawn cohort by
qc_cohort.py, which exits non-zero if any fails.
"""
import argparse
import hashlib
import json
import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

SEED = 20260924
C_MIN = 100
C_TRAIN = 300
ORIGINS = ("tumor_primary", "normal_adjacent")
T_TYPES = ("T cell CD4", "T cell CD8")
DISEASE = "lung adenocarcinoma"
OBS_COLS = ["disease", "origin", "study", "assay", "donor_id", "cell_type_major",
            "doublet_status", "sample", "platform"]


def sha256_file(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def obs_frame(f):
    index = f["obs"][f["obs"].attrs["_index"]][:].astype(str)
    out = {}
    for name in OBS_COLS:
        g = f["obs"][name]
        if isinstance(g, h5py.Group):
            cats = g["categories"][:].astype(str)
            codes = g["codes"][:]
            out[name] = np.where(codes >= 0, cats[np.clip(codes, 0, None)], None)
        else:
            out[name] = g[:].astype(str)
    df = pd.DataFrame(out)
    df["cell_id"] = index
    df["row"] = np.arange(len(df))
    return df


def eligible_cells(obs, paired_donors):
    m = ((obs.disease == DISEASE) & obs.cell_type_major.isin(T_TYPES)
         & (obs.doublet_status == "singlet") & obs.origin.isin(ORIGINS)
         & obs.assay.astype(str).str.startswith("10x") & obs.donor_id.notna())
    e = obs[m].copy()
    counts = e.groupby(["donor_id", "origin"]).size().unstack(fill_value=0)
    qualifying = counts[(counts[list(ORIGINS)] >= C_MIN).all(axis=1)].index
    if set(qualifying) != set(paired_donors):
        raise SystemExit(
            f"Qualifying donors in this file ({len(qualifying)}) differ from the "
            f"registered Phase 1 list ({len(paired_donors)}): "
            f"extra={sorted(set(qualifying) - set(paired_donors))} "
            f"missing={sorted(set(paired_donors) - set(qualifying))}")
    return e[e.donor_id.isin(paired_donors)]


def draw(e, seed=SEED):
    parts = []
    for (donor, origin), g in e.groupby(["donor_id", "origin"], sort=True):
        g = g.sort_values("cell_id")
        key = hashlib.sha256(f"{seed}|{donor}|{origin}".encode()).hexdigest()
        rng = np.random.default_rng(int(key[:16], 16))
        rank = np.empty(len(g), dtype=int)
        rank[rng.permutation(len(g))] = np.arange(len(g))
        g = g.assign(draw_rank=rank)
        parts.append(g)
    c = pd.concat(parts)
    c["in_analysis100"] = c.draw_rank < C_MIN
    c["in_pool300"] = c.draw_rank < C_TRAIN
    return c


def write_slim(f, cohort, out_path):
    pool = cohort[cohort.in_pool300].sort_values("row")
    rows = pool.row.to_numpy()
    X = f["raw"]["X"]
    indptr = X["indptr"][:]
    data, indices, ptr = [], [], [0]
    for r in rows:
        a, b = indptr[r], indptr[r + 1]
        data.append(X["data"][a:b])
        indices.append(X["indices"][a:b])
        ptr.append(ptr[-1] + (b - a))
    n_genes = int(X.attrs["shape"][1])
    m = sp.csr_matrix((np.concatenate(data), np.concatenate(indices), np.array(ptr)),
                      shape=(len(rows), n_genes))
    if not np.array_equal(m.data, np.round(m.data)):
        raise SystemExit("raw/X values for the drawn cells are not integer-valued")
    m = m.astype(np.float32)
    var_index = f["raw"]["var"][f["raw"]["var"].attrs["_index"]][:].astype(str)
    obs = pd.DataFrame({
        "cell_id": pool.cell_id.to_numpy(),
        "individual": pool.donor_id.to_numpy(),
        "origin": pool.origin.to_numpy(),
        "celltype": pool.cell_type_major.to_numpy(),
        "study": pool.study.to_numpy(),
        "in_analysis100": pool.in_analysis100.astype(int).to_numpy(),
        "n_counts": np.asarray(m.sum(axis=1)).ravel().astype(float),
    }, index=pool.cell_id.to_numpy())
    var = pd.DataFrame({"ensembl_id": var_index}, index=var_index)
    ad.AnnData(X=m, obs=obs, var=var).write_h5ad(out_path)
    return len(rows), int(m.nnz)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--h5ad", required=True)
    p.add_argument("--manifest", required=True, help="MANIFEST.json from verify_download.py")
    p.add_argument("--paired-donors", required=True, help="Phase 1 luca_luad_paired_10x_donors.csv")
    p.add_argument("--out-dir", required=True, help="balanced_donor_luad/ in the repo")
    p.add_argument("--slim-dir", required=True)
    a = p.parse_args()

    manifest = json.load(open(a.manifest))
    if not manifest.get("all_passed"):
        raise SystemExit("MANIFEST.json does not record a passing verification")
    sha = sha256_file(a.h5ad)
    if sha != manifest["sha256"]:
        raise SystemExit(f"h5ad sha256 {sha} != manifest {manifest['sha256']}")

    donors = sorted(pd.read_csv(a.paired_donors).donor_id.astype(str))
    with h5py.File(a.h5ad, "r") as f:
        obs = obs_frame(f)
        e = eligible_cells(obs, donors)
        cohort = draw(e)
        os.makedirs(a.slim_dir, exist_ok=True)
        slim = os.path.join(a.slim_dir, "balanced_donor_luad_pool300.h5ad")
        n_rows, nnz = write_slim(f, cohort, slim)

    cols = ["cell_id", "donor_id", "study", "origin", "cell_type_major", "assay",
            "sample", "draw_rank", "in_analysis100", "in_pool300"]
    eligible_per_side = cohort.groupby(["donor_id", "origin"]).size()
    cohort = cohort[cohort.in_pool300].sort_values(["donor_id", "origin", "draw_rank"])
    csv = os.path.join(a.out_dir, "cohort", "cohort_cells.csv")
    cohort[cols].to_csv(csv, index=False)

    definition = {
        "registered_rule": "PHASE0_SELECTION_RULE.md (sha256 58e9b153..., amendments 1, 1a)",
        "source_h5ad_sha256": sha,
        "source_dataset_version_id": manifest["dataset_version_id"],
        "seed": SEED,
        "draw": "per (donor, origin): sort cell_id, permute with default_rng(int(sha256(f'{seed}|{donor}|{origin}')[:16],16)); rank fixes both cohorts",
        "filters": {"disease": DISEASE, "cell_type_major": list(T_TYPES),
                    "doublet_status": "singlet", "assay": "startswith 10x",
                    "origin": list(ORIGINS), "qualifying": f">= {C_MIN} cells on both origins"},
        "c_min_analysis": C_MIN,
        "c_train_pool": C_TRAIN,
        "n_paired_donors": len(donors),
        "n_eligible_cells": int(eligible_per_side.sum()),
        "eligible_cells_per_donor_origin": {f"{d}|{o}": int(n) for (d, o), n in eligible_per_side.items()},
        "n_analysis100": int(cohort.in_analysis100.sum()),
        "n_pool300": int(cohort.in_pool300.sum()),
        "cohort_cells_csv_sha256": sha256_file(csv),
        "slim_h5ad": {"path": slim, "n_cells": n_rows, "nnz": nnz,
                      "sha256": sha256_file(slim), "committed": False},
    }
    with open(os.path.join(a.out_dir, "cohort", "COHORT_DEFINITION.json"), "w") as fh:
        json.dump(definition, fh, indent=1)
    print(json.dumps(definition, indent=1))


if __name__ == "__main__":
    main()
