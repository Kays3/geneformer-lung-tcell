"""Phase 3c: tokenise the training-pool cohort with Geneformer V2, then QC it.

Settings mirror current_workflow (July runner): TranscriptomeTokenizer,
model_version="V2" (model_input_size 4096, special tokens on), raw integer
counts in X, var "ensembl_id", obs "n_counts". Nothing is filtered here: the
cohort was fixed before tokenisation, so every drawn cell must come out.

Input goes through Geneformer's LOOM path, not its h5ad path: the installed
tokenizer's h5ad branch indexes a string-indexed pandas Series with integer
positions, which pandas 3 rejects (KeyError in ensembl_id_collapsed). The slim
h5ad is converted to loom here and the conversion is checked exactly (same
matrix, gene IDs and cell attributes) before tokenising.

QC (written to qc/QC_TOKENIZATION.md, exit non-zero on failure):
  - the slim input matches the sha256 recorded in COHORT_DEFINITION.json;
  - tokenised cell_id set == cohort_cells.csv cell_id set (no cell lost or added);
  - no empty cell; no cell longer than the model input size;
  - share of the input genes that map to the V2 token dictionary.
"""
import argparse
import hashlib
import json
import os
import pickle
import sys

import numpy as np
import pandas as pd

ATTRS = {"cell_id": "cell_id", "individual": "individual", "origin": "origin",
         "celltype": "celltype", "study": "study", "in_analysis100": "in_analysis100"}
PREFIX = "balanced_donor_luad_pool300"


def sha256_file(path, block=1 << 24):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(block), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--slim-dir", required=True)
    p.add_argument("--out-dir", required=True)
    p.add_argument("--repo-dir", required=True, help="balanced_donor_luad/ in the repo")
    p.add_argument("--nproc", type=int, default=8)
    a = p.parse_args()

    definition = json.load(open(os.path.join(a.repo_dir, "cohort", "COHORT_DEFINITION.json")))
    slim = os.path.join(a.slim_dir, f"{PREFIX}.h5ad")
    sha = sha256_file(slim)
    if sha != definition["slim_h5ad"]["sha256"]:
        raise SystemExit(f"slim h5ad sha256 {sha} != COHORT_DEFINITION {definition['slim_h5ad']['sha256']}")

    import anndata as ad
    from geneformer import TranscriptomeTokenizer
    from geneformer.tokenizer import TOKEN_DICTIONARY_FILE
    from datasets import load_from_disk

    out_ds = os.path.join(a.out_dir, f"{PREFIX}.dataset")
    if os.path.exists(out_ds):
        raise SystemExit(f"{out_ds} exists; remove it deliberately before re-tokenising")

    import loompy
    loom_dir = os.path.join(a.out_dir, "loom_input")
    os.makedirs(loom_dir, exist_ok=True)
    loom = os.path.join(loom_dir, f"{PREFIX}.loom")
    src = ad.read_h5ad(slim)
    col_attrs = {k: src.obs[k].to_numpy() for k in ATTRS}
    col_attrs["n_counts"] = src.obs["n_counts"].to_numpy()
    loompy.create(loom, src.X.T.tocoo(), row_attrs={"ensembl_id": src.var["ensembl_id"].to_numpy()},
                  col_attrs=col_attrs)
    with loompy.connect(loom, "r") as ds_l:
        conv_ok = (ds_l.shape == (src.n_vars, src.n_obs)
                   and np.array_equal(ds_l.ra["ensembl_id"], src.var["ensembl_id"].to_numpy())
                   and np.array_equal(ds_l.ca["cell_id"], src.obs["cell_id"].to_numpy())
                   and ds_l.sparse().sum(dtype=np.float64) == src.X.sum(dtype=np.float64)
                   and (ds_l.sparse().T.tocsr() != src.X.tocsr()).nnz == 0)
    if not conv_ok:
        raise SystemExit("h5ad -> loom conversion changed the data")

    tk = TranscriptomeTokenizer(ATTRS, nproc=a.nproc, chunk_size=512, model_version="V2")
    tk.tokenize_data(data_directory=loom_dir, output_directory=a.out_dir,
                     output_prefix=PREFIX, file_format="loom")

    ds = load_from_disk(out_ds)
    cohort = pd.read_csv(os.path.join(a.repo_dir, "cohort", "cohort_cells.csv"))
    tok_ids = set(ds["cell_id"])
    want = set(cohort.cell_id)
    lengths = np.array(ds["length"])
    token_dict = pickle.load(open(TOKEN_DICTIONARY_FILE, "rb"))
    var = ad.read_h5ad(slim, backed="r").var
    mapped = float(var.ensembl_id.isin(token_dict.keys()).mean())
    n_analysis = int(np.sum(np.array(ds["in_analysis100"]).astype(int)))

    rows = [
        ("slim input sha256 matches COHORT_DEFINITION", sha[:16] + "...", "match", True),
        ("h5ad -> loom conversion exact (matrix, gene IDs, cell IDs)", conv_ok, "True", conv_ok),
        ("tokenised cells", len(ds), f"== {len(want)}", len(ds) == len(want)),
        ("tokenised cell_id set == cohort cell_id set", tok_ids == want, "True", tok_ids == want),
        ("analysis cells flagged in tokenised data", n_analysis, f"== {int(cohort.in_analysis100.sum())}",
         n_analysis == int(cohort.in_analysis100.sum())),
        ("min tokens per cell", int(lengths.min()), "> 2 (more than the two special tokens)",
         int(lengths.min()) > 2),
        ("max tokens per cell", int(lengths.max()), "<= 4096", int(lengths.max()) <= 4096),
        ("median tokens per cell", float(np.median(lengths)), "report only", True),
        ("input genes present in V2 token dictionary", f"{mapped:.4f}", "report only", True),
    ]
    lines = ["# Tokenisation QC", "", "Computed values. Generated by `scripts/tokenize_cohort.py`.", "",
             "| Check | Computed | Threshold | Pass |", "|---|---|---|---|"]
    lines += [f"| {n} | {v} | {t} | {'PASS' if ok else '**FAIL**'} |" for n, v, t, ok in rows]
    ok = all(r[3] for r in rows)
    lines += ["", f"All passed: **{ok}**", "",
              f"Tokenizer: geneformer TranscriptomeTokenizer, model_version=V2, "
              f"token dictionary `{os.path.basename(TOKEN_DICTIONARY_FILE)}`.", ""]
    text = "\n".join(lines)
    open(os.path.join(a.repo_dir, "qc", "QC_TOKENIZATION.md"), "w").write(text)
    print(text)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
