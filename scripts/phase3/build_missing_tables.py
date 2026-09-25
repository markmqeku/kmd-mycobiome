#!/usr/bin/env python
"""Build the three outstanding supplementary tables: S2 Bloemfontein inventory,
S3 taxonomy threshold comparison, S4 DADA2 statistics for all 45 samples."""
import csv, zipfile, os
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; A = f"{O}/assembly"; os.makedirs(A, exist_ok=True)
RANKS=[("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
       ("f__","family"),("g__","genus"),("s__","species")]
CUT={"species":97.,"genus":95.,"family":90.,"order":85.,"class":80.,"phylum":80.,"kingdom":80.}
PH={"unidentified","unclassified","incertae_sedis","unknown",""}
def ranks_of(t):
    out={r:"" for _,r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part=part.strip()
        for pre,rk in RANKS:
            if part.startswith(pre):
                v=part[len(pre):].strip(); b=v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk=="species" and (b.endswith("_sp") or b.endswith("_sp.")))): v=""
                out[rk]=v
    return out
def load_tax(p):
    d={}
    for r in csv.reader(open(p),delimiter="\t"):
        if r and r[0]!="Feature ID" and not r[0].startswith("#"): d[r[0]]=r[1]
    return d
pid={}
z=zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn=[n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p=line.split("\t")
    if len(p)<3: continue
    try: v=float(p[2])
    except ValueError: continue
    if p[0] not in pid or v>pid[p[0]]: pid[p[0]]=v

meta={}
rd=csv.reader(open(f"{O}/metadata_v2_27-07-2026.tsv"),delimiter="\t"); h=next(rd); next(rd)
di={c:i for i,c in enumerate(h)}
for r in rd:
    if r and r[0].strip(): meta[r[0]]={c:r[di[c]] for c in h}
def load_table(p):
    lines=[l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hd=lines[0].lstrip("#").split("\t"); cols=[x for x in hd[1:] if x.strip()]
    per={s:{} for s in cols}
    for l in lines[1:]:
        q=l.split("\t")
        for s,v in zip(cols,q[1:1+len(cols)]):
            v=int(float(v))
            if v: per[s][q[0]]=v
    return cols,per
cols_f,per_f=load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
cols_a,per_a=load_table(f"{O}/table_exp/feature-table.tsv")
host={l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
final=load_tax(f"{O}/final_taxonomy.tsv")
yeast={r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"),delimiter="\t")
       if r and r[0]!="ASV_ID_md5" and r[2]=="1"}
def label(a):
    if a in yeast: return "Tremellomycetes yeast lineage (ITS2-unresolved)"
    rk=ranks_of(final.get(a,"")); last=("none","")
    for _,r in RANKS:
        if rk[r]: last=(r,rk[r])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"

# ---- S2 Bloemfontein inventory ----
bl=[s for s in cols_a if meta.get(s) and meta[s]["location"]=="Bloemfontein"]
sel=[]
for s in bl:
    tot=sum(per_a[s].values()); hs=sum(v for a,v in per_a[s].items() if a in host)
    hp=100*hs/tot if tot else 0
    if hp<=10: sel.append((s,hp,sum(per_f.get(s,{}).values())))
present=Counter()
for s,_,_ in sel:
    for a in per_f.get(s,{}): present[label(a)]+=1
with open(f"{A}/TableS2_bloemfontein_inventory.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["# Bloemfontein presence-absence inventory. Low-host samples only (host <=10%)."])
    w.writerow(["# No abundances, no diversity metrics, no cross-site comparison."])
    w.writerow(["sample","tissue","age","host_pct","fungal_reads"])
    for s,hp,fr in sorted(sel,key=lambda x:x[1]):
        w.writerow([s,meta[s]["tissue"],meta[s]["age"],f"{hp:.2f}",fr])
    w.writerow([]); w.writerow(["taxon_present","n_ASVs"])
    for k,v in present.most_common(): w.writerow([k,v])
print(f"S2: {len(sel)} samples, {len(present)} taxa")

# ---- S3 taxonomy thresholds ----
t82=load_tax(f"{B}/phase2/taxonomy/exp/taxonomy.tsv")
t10=load_tax(f"{O}/taxonomy_unite10/exp/taxonomy.tsv")
fungi=set(final)
with open(f"{A}/TableS3_taxonomy_thresholds.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["# ASVs assigned (%) at each rank, fungi-only table (n=1173)"])
    w.writerow(["rank","UNITE_8.2_permissive","UNITE_10.0_permissive","UNITE_10.0_rank_threshold"])
    for _,rk in RANKS:
        a=sum(1 for x in fungi if ranks_of(t82.get(x,""))[rk])
        b=sum(1 for x in fungi if ranks_of(t10.get(x,""))[rk])
        c=sum(1 for x in fungi if ranks_of(final.get(x,""))[rk])
        w.writerow([rk,f"{100*a/len(fungi):.1f}",f"{100*b/len(fungi):.1f}",f"{100*c/len(fungi):.1f}"])
print("S3: written")

# ---- S4 DADA2 stats ----
rows=[]
for grp,p in (("A_Bloemfontein","groupA2"),("B_Christiana_Pretoria","groupB"),("B_OT-P","otp")):
    fp=f"{B}/phase2/{p}/stats_exp/stats.tsv"
    if not os.path.exists(fp): continue
    with open(fp) as fh:
        r=csv.reader(fh,delimiter="\t"); hh=next(r); next(r)
        for row in r:
            if row and row[0].strip(): rows.append([grp]+row)
with open(f"{A}/TableS4_dada2_stats.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["# DADA2 statistics, all denoised samples. Group A trunc 228/200; Group B trunc 280/230."])
    w.writerow(["denoise_group","sample-id","input","filtered","pct_filtered","denoised",
                "merged","pct_merged","non-chimeric","pct_non_chimeric"])
    w.writerows(rows)
print(f"S4: {len(rows)} samples")
