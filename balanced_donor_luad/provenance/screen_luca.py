import pandas as pd, numpy as np
C_MIN=100
o=pd.read_parquet("luca_obs.parquet")
o=o[(o.doublet_status=="singlet")&o.donor_id.notna()]
Tdef={"primary(CD4+CD8)":["T cell CD4","T cell CD8"],"sens(+Treg)":["T cell CD4","T cell CD8","T cell regulatory"]}
for tname,types in Tdef.items():
  T=o[o.cell_type_major.isin(types)]
  print(f"\n######## T definition: {tname}")
  for dis in ["lung adenocarcinoma","squamous cell lung carcinoma","non-small cell lung carcinoma"]:
    d=T[T.disease==dis]
    cnt=d.groupby(["study","platform","origin","donor_id"],observed=True).size().rename("n").reset_index()
    q=cnt[cnt.n>=C_MIN]
    tp=q[q.origin=="tumor_primary"]; na=q[q.origin=="normal_adjacent"]
    per=[]
    for s,g in q.groupby("study",observed=True):
        t=set(g[g.origin=="tumor_primary"].donor_id); n=set(g[g.origin=="normal_adjacent"].donor_id)
        per.append((s,";".join(sorted(g.platform.astype(str).unique())),len(t),len(n),len(t&n)))
    per=pd.DataFrame(per,columns=["study","platform","n_tumor","n_normadj","paired"])
    both=per[(per.n_tumor>0)&(per.n_normadj>0)]
    print(f"\n== {dis}: qualifying donors (>= {C_MIN} T cells) tumor_primary={tp.donor_id.nunique()}  normal_adjacent={na.donor_id.nunique()}  paired(same donor both)={int(per.paired.sum())}")
    print(per.to_string(index=False))
    # F3-compliant subset: studies contributing to both groups
    print(f"   F3 subset (studies with both groups): tumor={both.n_tumor.sum()} normadj={both.n_normadj.sum()} paired={both.paired.sum()}")
    # raw imbalance before capping, for the record: max share / max-min within paired donors
    pd_=set(q[q.origin=="tumor_primary"].donor_id)&set(q[q.origin=="normal_adjacent"].donor_id)
    for org in ["tumor_primary","normal_adjacent"]:
        x=cnt[(cnt.origin==org)&cnt.donor_id.isin(pd_)].n
        if len(x): print(f"   raw {org} over paired donors: cells min={x.min()} median={int(x.median())} max={x.max()} max/min={x.max()/x.min():.1f} max share={x.max()/x.sum():.1%}")
