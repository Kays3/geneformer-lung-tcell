"""E2 step 1: build the balanced paired-donor colorectal cohort from Pelka et al. 2021 (GSE178341).

Selection (registration s.2, fixed before any expression value is read):
  cells     clMidwayPr in {TCD4, TCD8} (author annotation, cluster table)
            PROCESSING_TYPE == "unsorted" (both tissues; the 19-donor unsorted-only cohort)
            SPECIMEN_TYPE T -> origin "tumor_primary", N -> origin "normal_adjacent"
  donors    PID with >= 100 such cells in BOTH tissues. The recomputed set must equal the
            19 registered PIDs exactly, else the script refuses.
The draw is the balanced-donor one, unchanged (build_cohort.py): per (donor, origin), sort cell IDs,
permute with default_rng(int(sha256(f"{seed}|{donor}|{origin}")[:16], 16)); rank < 100 is the
analysis cohort, rank < 300 the training pool.

Genes: GSE178341 feature IDs carry a version and a suffix (e.g. ENSG00000243485.5_4); 28 pseudoautosomal
copies also end in _PAR_Y. Version and suffixes are stripped; IDs that collide after stripping (exactly the
28 PAR_Y copies onto their X-chromosome IDs) are summed (counted and reported).
The h5 also holds 35 non-gene features labelled "Gene Expression" (4 cell hashtags HASH-S1..S4 and
31 antibody-derived tags, e.g. CD4-SK3). Every feature whose ID does not start with ENSG is dropped
before anything else; their names are recorded in COHORT_DEFINITION.json.

Outputs: cohort/cohort_cells.csv, cohort/COHORT_DEFINITION.json (committed);
<slim_dir>/pelka_crc_e2_pool300.h5ad raw integer counts (not committed, hash recorded).
"""
import argparse
import hashlib
import json
import os
import re

import numpy as np
import pandas as pd

SEED = 20260924
C_MIN = 100
C_TRAIN = 300
T_TYPES = ("TCD4", "TCD8")
PROCESSING = "unsorted"
ORIGIN_OF = {"T": "tumor_primary", "N": "normal_adjacent"}
STUDY = "Pelka_2021"
PREFIX = "pelka_crc_e2_pool300"
REGISTERED_DONORS = ("C107", "C110", "C111", "C123", "C125", "C126", "C129", "C130", "C132", "C134",
                     "C135", "C137", "C140", "C142", "C143", "C155", "C157", "C162", "C170")
SOURCE_SHA256 = {
    "GSE178341_crc10x_full_c295v4_submit.h5": "f435bb2651ff5297d0c24a99daf58850ed67ae1ed6c5ef05fad48fa3f0186670",
    "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz": "42ed840feb62aec786f092a18ea8c0cb88d07f919df207cc8f12dc2b08e96699",
    "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz": "0452149b0bcff9660951de87913b82a4d72784c08d0c782e49395f756a4b9853",
}
_VERSION = re.compile(r"\.\d+(_\d+)?(_PAR_Y)?$")


def sha256_file(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def strip_version(ids):
    out = [_VERSION.sub("", i) for i in ids]
    bad = [i for i in out if not re.fullmatch(r"ENSG\d{11}", i)]
    if bad:
        raise SystemExit(f"{len(bad)} feature IDs are not plain Ensembl gene IDs after stripping: {bad[:5]}")
    return out


def cell_table(cluster_csv, meta_csv):
    c = pd.read_csv(cluster_csv)
    m = pd.read_csv(meta_csv)
    if not (c.sampleID.to_numpy() == m.cellID.to_numpy()).all():
        raise SystemExit("cluster and metatables rows are not in the same cell order")
    x = pd.concat([c, m.drop(columns="cellID")], axis=1).rename(columns={"sampleID": "cell_id"})
    return x


def eligible_cells(x):
    m = x.clMidwayPr.isin(T_TYPES) & (x.PROCESSING_TYPE == PROCESSING) & x.SPECIMEN_TYPE.isin(list(ORIGIN_OF))
    e = x[m].copy()
    e["origin"] = e.SPECIMEN_TYPE.map(ORIGIN_OF)
    e["donor_id"] = e.PID.astype(str)
    counts = e.groupby(["donor_id", "origin"]).size().unstack(fill_value=0)
    for o in ORIGIN_OF.values():
        if o not in counts:
            counts[o] = 0
    qualifying = sorted(counts[(counts[list(ORIGIN_OF.values())] >= C_MIN).all(axis=1)].index)
    if tuple(qualifying) != REGISTERED_DONORS:
        raise SystemExit(f"qualifying donors {qualifying} != registered {list(REGISTERED_DONORS)}")
    return e[e.donor_id.isin(qualifying)]


def draw(e, seed=SEED):
    parts = []
    for (donor, origin), g in e.groupby(["donor_id", "origin"], sort=True):
        g = g.sort_values("cell_id")
        key = hashlib.sha256(f"{seed}|{donor}|{origin}".encode()).hexdigest()
        rng = np.random.default_rng(int(key[:16], 16))
        rank = np.empty(len(g), dtype=int)
        rank[rng.permutation(len(g))] = np.arange(len(g))
        parts.append(g.assign(draw_rank=rank))
    c = pd.concat(parts)
    c["in_analysis100"] = c.draw_rank < C_MIN
    c["in_pool300"] = c.draw_rank < C_TRAIN
    return c


def collapse_genes(ids):
    """Return (unique_ids, column index of each original feature into unique_ids, n_collapsed)."""
    uniq, inv = np.unique(np.asarray(ids), return_inverse=True)
    return uniq, inv, int(len(ids) - len(uniq))


def slim_matrix(h5_path, cells):
    """Raw counts for `cells` (DataFrame with column cell_id, in output order) as CSR cells x genes."""
    import h5py
    import scipy.sparse as sp
    with h5py.File(h5_path, "r") as f:
        g = f["matrix"]
        barcodes = g["barcodes"][:].astype(str)
        ftype = set(g["features"]["feature_type"][:].astype(str))
        if ftype != {"Gene Expression"}:
            raise SystemExit(f"unexpected feature types {ftype}")
        raw_ids = g["features"]["id"][:].astype(str)
        n_genes, n_cells = (int(v) for v in g["shape"][:])
        if len(barcodes) != n_cells or len(raw_ids) != n_genes:
            raise SystemExit("h5 shape disagrees with barcodes/features")
        keep = np.char.startswith(raw_ids, "ENSG")
        dropped = sorted(raw_ids[~keep].tolist())
        ids = strip_version(raw_ids[keep].tolist())
        pos = pd.Series(np.arange(n_cells), index=barcodes)
        if not pos.index.is_unique:
            raise SystemExit("h5 barcodes are not unique")
        missing = set(cells.cell_id) - set(pos.index)
        if missing:
            raise SystemExit(f"{len(missing)} selected cells absent from the h5, e.g. {sorted(missing)[:3]}")
        cols = pos.loc[cells.cell_id].to_numpy()
        indptr = g["indptr"][:]
        # Read each run of nearby selected cells as one contiguous slice (per-cell reads decompress the same
        # h5 chunks thousands of times); values are identical, only the read pattern changes.
        order = np.sort(np.unique(cols))
        runs, start = [], 0
        for i in range(1, len(order) + 1):
            if i == len(order) or order[i] - order[i - 1] > 2000:
                runs.append(order[start:i]); start = i
        block = {}
        for r in runs:
            lo, hi = indptr[r[0]], indptr[r[-1] + 1]
            dat, ind = g["data"][lo:hi], g["indices"][lo:hi]
            for c in r:
                a, b = indptr[c] - lo, indptr[c + 1] - lo
                block[c] = (dat[a:b], ind[a:b])
        data, indices, ptr = [], [], [0]
        for c in cols:
            dv, iv = block[c]
            data.append(dv)
            indices.append(iv)
            ptr.append(ptr[-1] + len(dv))
    m = sp.csr_matrix((np.concatenate(data), np.concatenate(indices), np.array(ptr)), shape=(len(cols), n_genes))
    m = m[:, np.flatnonzero(keep)].tocsr()
    if not np.array_equal(m.data, np.round(m.data)) or (m.data < 0).any():
        raise SystemExit("selected counts are not non-negative integers")
    uniq, inv, n_coll = collapse_genes(ids)
    n_keep = len(ids)
    agg = sp.csr_matrix((np.ones(n_keep), (np.arange(n_keep), inv)), shape=(n_keep, len(uniq)))
    out = (m @ agg).tocsr().astype(np.float32)
    out.sum_duplicates()
    if out.sum(dtype=np.float64) != m.sum(dtype=np.float64):
        raise SystemExit("gene collapse changed the total count")
    return out, uniq, n_coll, n_genes, dropped


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--raw-dir", required=True)
    p.add_argument("--out-dir", required=True, help="pelka_crc_e2/ in the repo")
    p.add_argument("--slim-dir", required=True)
    a = p.parse_args()

    for fn, want in SOURCE_SHA256.items():
        got = sha256_file(os.path.join(a.raw_dir, fn))
        if got != want:
            raise SystemExit(f"{fn}: sha256 {got} != registered {want}")
    x = cell_table(os.path.join(a.raw_dir, "GSE178341_crc10x_full_c295v4_submit_cluster.csv.gz"),
                   os.path.join(a.raw_dir, "GSE178341_crc10x_full_c295v4_submit_metatables.csv.gz"))
    e = eligible_cells(x)
    cohort = draw(e)
    pool = cohort[cohort.in_pool300].sort_values(["donor_id", "origin", "draw_rank"])

    import anndata as ad
    m, genes, n_coll, n_feat, dropped = slim_matrix(os.path.join(a.raw_dir, "GSE178341_crc10x_full_c295v4_submit.h5"), pool)
    obs = pd.DataFrame({
        "cell_id": pool.cell_id.to_numpy(), "individual": pool.donor_id.to_numpy(),
        "origin": pool.origin.to_numpy(), "celltype": pool.clMidwayPr.to_numpy(),
        "study": STUDY, "in_analysis100": pool.in_analysis100.astype(int).to_numpy(),
        "n_counts": np.asarray(m.sum(axis=1)).ravel().astype(float),
    }, index=pool.cell_id.to_numpy())
    var = pd.DataFrame({"ensembl_id": genes}, index=genes)
    os.makedirs(a.slim_dir, exist_ok=True)
    slim = os.path.join(a.slim_dir, f"{PREFIX}.h5ad")
    ad.AnnData(X=m, obs=obs, var=var).write_h5ad(slim)

    pool = pool.assign(study=STUDY)
    cols = ["cell_id", "donor_id", "study", "origin", "clMidwayPr", "cl295v11SubShort", "batchID", "SINGLECELL_TYPE",
            "PROCESSING_TYPE", "MMRStatus", "draw_rank", "in_analysis100", "in_pool300"]
    os.makedirs(os.path.join(a.out_dir, "cohort"), exist_ok=True)
    csv = os.path.join(a.out_dir, "cohort", "cohort_cells.csv")
    pool[cols].to_csv(csv, index=False)
    per_side = cohort.groupby(["donor_id", "origin"]).size()
    donor_meta = pool.groupby("donor_id").agg(chemistry=("SINGLECELL_TYPE", lambda s: "/".join(sorted(set(s)))),
                                              mmr=("MMRStatus", "first"))
    definition = {
        "registration": "pelka_crc_e2/registration/E2_REGISTRATION.md s.2",
        "source_sha256": SOURCE_SHA256,
        "source": "GEO GSE178341 (Pelka et al. 2021, Cell), downloaded 2026-09-30 15:06-15:07 UTC",
        "seed": SEED,
        "draw": "per (donor, origin): sort cell_id, permute with default_rng(int(sha256(f'{seed}|{donor}|{origin}')[:16],16)); rank fixes both cohorts",
        "filters": {"clMidwayPr": list(T_TYPES), "PROCESSING_TYPE": PROCESSING, "origin_map": ORIGIN_OF,
                    "qualifying": f">= {C_MIN} cells on both origins"},
        "c_min_analysis": C_MIN, "c_train_pool": C_TRAIN,
        "donors": list(REGISTERED_DONORS), "n_paired_donors": len(REGISTERED_DONORS),
        "donor_chemistry": donor_meta.chemistry.to_dict(), "donor_mmr": donor_meta.mmr.to_dict(),
        "n_eligible_cells": int(per_side.sum()),
        "eligible_cells_per_donor_origin": {f"{d}|{o}": int(n) for (d, o), n in per_side.items()},
        "n_analysis100": int(pool.in_analysis100.sum()), "n_pool300": int(len(pool)),
        "genes": {"n_features": n_feat, "n_non_gene_features_dropped": len(dropped),
                  "non_gene_features_dropped": dropped,
                  "n_after_version_strip_and_collapse": int(len(genes)), "n_collapsed": n_coll},
        "cohort_cells_csv_sha256": sha256_file(csv),
        "slim_h5ad": {"path": slim, "n_cells": int(m.shape[0]), "nnz": int(m.nnz), "sha256": sha256_file(slim),
                      "committed": False},
    }
    json.dump(definition, open(os.path.join(a.out_dir, "cohort", "COHORT_DEFINITION.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in definition.items() if k != "eligible_cells_per_donor_origin"}, indent=1))


if __name__ == "__main__":
    main()
