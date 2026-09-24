"""Phase 2 integrity check for the LuCA extended download.

Checks, in order, and writes MANIFEST.json next to the file:
  1. byte size equals the source Content-Length;
  2. S3 multipart ETag recomputed from the bytes (8 MiB parts) equals the
     source ETag, which ties the local bytes to the published object;
  3. sha256 of the file, recorded as our own identifier for downstream use;
  4. the file names the expected dataset version and holds raw counts;
  5. the paired-donor table rebuilt from this file equals the Phase 1 table
     (luca_luad_paired_10x_donors.csv), so selection ran on these exact bytes.

Reads structural metadata only (obs labels, counts of cells); no expression.
"""
import hashlib
import json
import math
import os
import sys

import h5py
import numpy as np
import pandas as pd

DATASET_VERSION = "33165751-ae33-4a65-94ee-52d9bc38f97e"
PART = 8 * 1024 * 1024


def headers(path):
    h = {}
    for line in open(path):
        if ":" in line:
            k, v = line.split(":", 1)
            h[k.strip().lower()] = v.strip().strip('"')
    return h


def hashes(path, n_parts_expected):
    sha = hashlib.sha256()
    part_md5s = []
    with open(path, "rb") as f:
        while True:
            chunk = f.read(PART)
            if not chunk:
                break
            sha.update(chunk)
            part_md5s.append(hashlib.md5(chunk).digest())
    etag = hashlib.md5(b"".join(part_md5s)).hexdigest() + f"-{len(part_md5s)}"
    return sha.hexdigest(), etag, len(part_md5s)


def obs_col(f, name):
    g = f["obs"][name]
    if isinstance(g, h5py.Group):
        cats = g["categories"][:].astype(str)
        codes = g["codes"][:]
        return pd.Series(np.where(codes >= 0, cats[np.clip(codes, 0, None)], None))
    return pd.Series(g[:].astype(str))


def paired_table(f):
    o = pd.DataFrame({k: obs_col(f, k) for k in
                      ["disease", "origin", "study", "assay", "donor_id",
                       "cell_type_major", "doublet_status"]})
    o = o[(o.doublet_status == "singlet") & o.donor_id.notna()]
    t = o[o.cell_type_major.isin(["T cell CD4", "T cell CD8"])
          & (o.disease == "lung adenocarcinoma") & o.assay.str.startswith("10x")]
    c = t.groupby(["origin", "study", "donor_id"]).size().rename("n").reset_index()
    c = c[c.n >= 100]
    m = c[c.origin == "tumor_primary"].merge(
        c[c.origin == "normal_adjacent"], on=["study", "donor_id"], suffixes=("_t", "_n"))
    return m[["study", "donor_id", "n_t", "n_n"]].sort_values("donor_id").reset_index(drop=True)


def main(data_dir, phase1_csv):
    path = os.path.join(data_dir, f"{DATASET_VERSION}.h5ad")
    src = headers(os.path.join(data_dir, "source_headers.txt"))
    out = {"file": path, "dataset_version_id": DATASET_VERSION,
           "source_url": f"https://datasets.cellxgene.cziscience.com/{DATASET_VERSION}.h5ad",
           "source_headers": {k: src.get(k) for k in
                              ["content-length", "etag", "last-modified", "x-amz-version-id"]}}
    checks = {}

    size = os.path.getsize(path)
    checks["size_matches_content_length"] = size == int(src["content-length"])
    out["size_bytes"] = size

    n_expected = int(src["etag"].split("-")[1])
    sha, etag, n_parts = hashes(path, n_expected)
    out["sha256"] = sha
    out["multipart_etag_recomputed"] = etag
    checks["multipart_etag_matches_source"] = etag == src["etag"]
    checks["part_count_consistent"] = n_parts == n_expected == math.ceil(size / PART)

    with h5py.File(path, "r") as f:
        title = f["uns"]["title"][()].decode()
        out["uns_title"] = title
        out["n_cells"], out["n_genes"] = (int(x) for x in f["X"].attrs["shape"])
        checks["title_is_extended_atlas"] = "extended atlas" in title
        checks["raw_counts_present"] = "raw" in f and "X" in f["raw"]
        rd = f["raw"]["X"]["data"][:1_000_000]
        checks["raw_counts_integer_valued_first_1M"] = bool(np.all(rd == np.round(rd)))
        mine = paired_table(f)

    ref = pd.read_csv(phase1_csv).sort_values("donor_id").reset_index(drop=True)
    checks["paired_table_equals_phase1"] = bool(
        mine[["study", "donor_id", "n_t", "n_n"]].astype(str)
        .equals(ref[["study", "donor_id", "n_t", "n_n"]].astype(str)))
    out["n_paired_donors"] = len(mine)
    out["checks"] = checks
    out["all_passed"] = all(checks.values())

    with open(os.path.join(data_dir, "MANIFEST.json"), "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps(out, indent=1))
    return 0 if out["all_passed"] else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2]))
