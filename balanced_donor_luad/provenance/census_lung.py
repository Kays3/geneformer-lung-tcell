import cellxgene_census, pandas as pd, os
from urllib.parse import urlparse
u=urlparse(os.environ["HTTPS_PROXY"])
cfg={"vfs.s3.proxy_scheme":"http","vfs.s3.proxy_host":u.hostname,"vfs.s3.proxy_port":str(u.port),"vfs.s3.proxy_username":u.username,"vfs.s3.proxy_password":u.password}
with cellxgene_census.open_soma(census_version="stable", tiledb_config=cfg) as c:
    print("census", c["census_info"]["summary"].read().concat().to_pandas().set_index("label").loc["census_build_date","value"])
    obs = cellxgene_census.get_obs(c, "homo_sapiens",
        value_filter="tissue_general == 'lung' and is_primary_data == True",
        column_names=["dataset_id","donor_id","disease","tissue","tissue_type","cell_type","assay","suspension_type"])
obs.to_parquet("lung_obs.parquet")
print(obs.shape)
