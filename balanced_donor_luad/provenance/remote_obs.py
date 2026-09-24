import sys, fsspec, h5py, anndata, pandas as pd
from anndata.io import read_elem
url, out = sys.argv[1], sys.argv[2]
fs = fsspec.filesystem("http", client_kwargs={"trust_env": True}, block_size=8*2**20)
with fs.open(url, "rb", cache_type="blockcache") as f, h5py.File(f, "r") as h:
    print("obs cols:", list(h["obs"].keys()))
    obs = read_elem(h["obs"])
obs.to_parquet(out); print(obs.shape)
