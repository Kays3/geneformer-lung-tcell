import pandas as pd, json
o=pd.read_parquet("lung_obs.parquet")
T=o[o.cell_type.str.contains(r"\bT cell|T-cell|thymocyte",regex=True) & ~o.cell_type.str.contains("natural killer T|NK T",regex=True)]
print("T cells in lung:",len(T), "celltypes:",T.cell_type.nunique())
ds={}
try:
    for d in json.load(open("datasets.json")):
        ds[d["dataset_id"]]=(d.get("title","")[:70], d.get("collection_id",""))
except Exception as e: print("ds meta err",e)
g=T.groupby(["dataset_id","disease","donor_id"],observed=True).size().rename("n").reset_index()
q=g[g.n>=100]
tum=q[q.disease!="normal"]; nor=q[q.disease=="normal"]
rows=[]
for did,sub in q.groupby("dataset_id"):
    dis=sub.groupby("disease").donor_id.nunique().sort_values(ascending=False)
    normal_don=set(sub[sub.disease=="normal"].donor_id)
    tumor_dis=[x for x in dis.index if x!="normal"]
    paired=len(normal_don & set(sub[sub.disease!="normal"].donor_id))
    rows.append(dict(dataset=did[:8],title=ds.get(did,("?",))[0],n_normal=dis.get("normal",0),
       top_disease=(tumor_dis[0] if tumor_dis else ""),n_top=(dis[tumor_dis[0]] if tumor_dis else 0),paired_same_donorid=paired,
       assays=";".join(sorted(T[T.dataset_id==did].assay.unique()))[:40]))
r=pd.DataFrame(rows).sort_values("n_top",ascending=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",70)
print(r[r.n_top>=6].to_string(index=False))
g.to_csv("lung_tcell_donor_counts.csv",index=False)
