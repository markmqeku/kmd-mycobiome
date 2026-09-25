#!/usr/bin/env python
"""Item 9 figures 1,2,4,5,6,7 + supplementary. 300 dpi, colourblind-safe (Okabe-Ito),
individual sample points everywhere (n = 2 to 6: no boxplots, no violins, no distribution summaries).
Group central tendency shown as a short horizontal line or a marked point only."""
import csv, zipfile, re, math, os
from collections import defaultdict, Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
F = f"{O}/figures"; os.makedirs(F, exist_ok=True)
DPI = 300
# Okabe-Ito colourblind-safe palette
CB = {"orange":"#E69F00","skyblue":"#56B4E9","green":"#009E73","yellow":"#F0E442",
      "blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7","black":"#000000","grey":"#999999"}
SITE_C = {"Bloemfontein":CB["orange"],"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M = {"reference_site":"s","asymptomatic":"o","symptomatic":"^"}
TISSUE_M = {"Leaves":"o","Twigs":"^","Inflorescence":"s","Seeds":"D","Malformation":"P"}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.dpi":DPI,"savefig.dpi":DPI,"savefig.bbox":"tight"})

# ---------------- data ----------------
def load_table(p):
    lines=[l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hdr=lines[0].lstrip("#").split("\t"); cols=[h for h in hdr[1:] if h.strip()]
    per={s:{} for s in cols}
    for l in lines[1:]:
        q=l.split("\t")
        for s,v in zip(cols,q[1:1+len(cols)]):
            v=int(float(v))
            if v: per[s][q[0]]=v
    return cols,per
cols_all,per_all = load_table(f"{O}/table_exp/feature-table.tsv")
cols_f,per_f     = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
host={l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
yeast={r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"),delimiter="\t")
       if r and r[0]!="ASV_ID_md5" and r[2]=="1"}
meta={}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd=csv.reader(fh,delimiter="\t"); h=next(rd); next(rd); di={c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]]={c:r[di[c]] for c in h}
core=[s for s in cols_f if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]
SRC={}

# ---------------- Figure 1: host DNA by site ----------------
tot={s:sum(per_all[s].values()) for s in cols_all}
hst={s:sum(v for a,v in per_all[s].items() if a in host) for s in cols_all}
pct={s:100*hst[s]/tot[s] if tot[s] else 0 for s in cols_all}
fig,ax=plt.subplots(figsize=(5.2,3.6))
sites=["Bloemfontein","Christiana","Pretoria"]
for i,site in enumerate(sites):
    ss=[s for s in cols_all if meta.get(s) and meta[s]["location"]==site]
    gm=100*sum(hst[s] for s in ss)/sum(tot[s] for s in ss)
    ax.bar(i,gm,width=.55,color=SITE_C[site],alpha=.35,edgecolor=SITE_C[site],linewidth=1.2,zorder=1)
    xs=[i+(hash(s)%100-50)/380 for s in ss]
    ax.scatter(xs,[pct[s] for s in ss],s=22,facecolor="white",edgecolor=SITE_C[site],
               linewidth=1.1,zorder=3)
    ax.text(i,gm+3,f"{gm:.1f}%",ha="center",fontsize=8.5,fontweight="bold")
ax.set_xticks(range(3)); ax.set_xticklabels([f"{s}\n(n={sum(1 for x in cols_all if meta.get(x) and meta[x]['location']==s)})" for s in sites])
ax.set_ylabel("Host (Searsia lancea) reads (% of sample)"); ax.set_ylim(-4,104)
ax.set_title("Figure 1. Host plant DNA content by site",fontsize=10,loc="left")
fig.savefig(f"{F}/Figure1_host_by_site.png"); plt.close(fig)
SRC["Figure 1"]="table_exp/feature-table.tsv; nontarget_screen/host_asv_ids_ALL.txt; metadata_v2_27-07-2026.tsv"

# ---------------- Figure 2: Bloemfontein host load vs fungal ASV count ----------------
bl=[s for s in cols_all if meta.get(s) and meta[s]["location"]=="Bloemfontein"]
nasv={s:sum(1 for a in per_f.get(s,{}) if per_f[s][a]>0) for s in cols_all}
fig,ax=plt.subplots(figsize=(6.6,4.8))
for t,m in TISSUE_M.items():
    ss=[s for s in bl if meta[s]["tissue"]==t]
    if not ss: continue
    ax.scatter([pct[s] for s in ss],[nasv.get(s,0) for s in ss],marker=m,s=52,
               facecolor="white",edgecolor=CB["black"],linewidth=1.1,label=t,zorder=3)
# stagger labels so the dense high-host cluster stays legible
_bl_sorted=sorted(bl,key=lambda s:(-pct[s],-nasv.get(s,0)))
_offsets=[(5,4),(5,-9),(-16,5),(-16,-9),(5,12),(-16,12)]
for i,s in enumerate(_bl_sorted):
    off=_offsets[i%len(_offsets)] if pct[s]>90 else (5,4)
    ax.annotate(s,(pct[s],nasv.get(s,0)),textcoords="offset points",xytext=off,
                fontsize=6.8,color=CB["grey"])
ax.set_xlabel("Host DNA in sample (%)"); ax.set_ylabel("Fungal ASVs recovered")
ax.set_title("Figure 2. Bloemfontein: apparent fungal richness tracks host load, not tissue",
             fontsize=10,loc="left")
ax.legend(title="Tissue",frameon=False,fontsize=8,title_fontsize=8,loc="upper right")
fig.savefig(f"{F}/Figure2_bloem_hostload_vs_ASVs.png"); plt.close(fig)
SRC["Figure 2"]="table_exp/feature-table.tsv; fungi_only/table_exp/feature-table.tsv; host_asv_ids_ALL.txt"

# ---------------- Figure 4: composition ----------------
RANKS=[("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
       ("f__","family"),("g__","genus"),("s__","species")]
CUT={"species":97.,"genus":95.,"family":90.,"order":85.,"class":80.,"phylum":80.,"kingdom":80.}
PH={"unidentified","unclassified","incertae_sedis","unknown",""}
tax={}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"),delimiter="\t"):
    if r and r[0]!="Feature ID" and not r[0].startswith("#"): tax[r[0]]=r[1]
z=zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn=[n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid={}
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p=line.split("\t")
    if len(p)<3: continue
    try: v=float(p[2])
    except ValueError: continue
    if p[0] not in pid or v>pid[p[0]]: pid[p[0]]=v
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
YL="Unresolved Tremellomycetes\nlineage"   # KMD PROMPT 4, R3a
def lab(a):
    if a in yeast: return YL
    rk=ranks_of(tax.get(a,"")); p=pid.get(a,0.); keep={}
    for _,r in RANKS:
        if rk[r] and p>=CUT[r]: keep[r]=rk[r]
        else: break
    last=(None,None)
    for _,r in RANKS:
        if keep.get(r): last=(r,keep[r])
    return f"{last[1]} ({last[0]})" if last[0] else "unassigned"
groups=[("Bloemfontein","reference_site"),("Christiana","asymptomatic"),("Christiana","symptomatic"),
        ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
comp={}; ns={}
for g in groups:
    ss=[s for s in cols_f if meta.get(s) and (meta[s]["location"],meta[s]["condition_std"])==g]
    ns[g]=len(ss); c=Counter(); t=0
    for s in ss:
        for a,v in per_f[s].items(): c[lab(a)]+=v; t+=v
    comp[g]={k:100*v/t for k,v in c.items()} if t else {}
top=set()
for g in groups: top|={k for k,_ in Counter(comp[g]).most_common(7)}
top=sorted(top,key=lambda k:-sum(comp[g].get(k,0) for g in groups))
palette=[CB["vermillion"],CB["blue"],CB["green"],CB["orange"],CB["purple"],CB["skyblue"],
         CB["yellow"],"#7F7F7F","#4D4D4D","#B2DF8A","#FB9A99","#CAB2D6"]
colmap={k:palette[i%len(palette)] for i,k in enumerate(top)}
fig,ax=plt.subplots(figsize=(7.6,4.6))
for i,g in enumerate(groups):
    bottom=0
    for k in top:
        v=comp[g].get(k,0)
        if v<=0: continue
        ax.bar(i,v,bottom=bottom,width=.62,color=colmap[k],edgecolor="white",linewidth=.4)
        bottom+=v
    other=100-bottom
    if other>0.01: ax.bar(i,other,bottom=bottom,width=.62,color="#DDDDDD",edgecolor="white",linewidth=.4)
ax.set_xticks(range(len(groups)))
ax.set_xticklabels([f"{g[0]}\n{g[1]}\n(n={ns[g]})" for g in groups],fontsize=8)
ax.set_ylabel("Relative abundance (% of fungal reads)"); ax.set_ylim(0,100)
ax.set_title("Figure 4. Taxonomic composition by site and condition",fontsize=10,loc="left")
h=[plt.Rectangle((0,0),1,1,color=colmap[k]) for k in top]+[plt.Rectangle((0,0),1,1,color="#DDDDDD")]
ax.legend(h,[k.replace("\n"," ") for k in top]+["other"],frameon=False,fontsize=7,
          bbox_to_anchor=(1.01,1),loc="upper left")
fig.savefig(f"{F}/Figure4_composition.png"); plt.close(fig)
SRC["Figure 4"]="fungi_only/table_exp/feature-table.tsv; taxonomy_unite10/exp/taxonomy.tsv; taxonomy_unite10/search.qza"

# ---------------- Figure 5: Hill numbers (parsed from H4.log) ----------------
hill=defaultdict(dict); cur=None
for line in open(f"{B}/phase3/H4.log",encoding="utf-8",errors="replace"):
    m=re.match(r"\s*--\s*q\s*=\s*(\d)",line)
    if m: cur=int(m.group(1)); continue
    m=re.match(r"\s*(\S+)\s+(Christiana|Pretoria)\s+(asymptomatic|symptomatic)\s+(\S+)\s+qD=\s*([\d.]+)",line)
    if m and cur is not None: hill[cur][m.group(1)]=float(m.group(5))
fig,axes=plt.subplots(1,3,figsize=(9.4,3.6),sharex=True)
titles={0:"q = 0  (richness)",1:"q = 1  (exp Shannon)",2:"q = 2  (inverse Simpson)"}
gg=[("Christiana","asymptomatic"),("Christiana","symptomatic"),("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
for qi,ax in zip((0,1,2),axes):
    for i,g in enumerate(gg):
        ss=[s for s in core if (meta[s]["location"],meta[s]["condition_std"])==g and s in hill[qi]]
        vals=[hill[qi][s] for s in ss]
        if not vals: continue
        xs=[i+(hash(s)%100-50)/420 for s in ss]
        ax.scatter(xs,vals,s=30,facecolor="white",edgecolor=SITE_C[g[0]],
                   marker=COND_M[g[1]],linewidth=1.2,zorder=3)
        mu=sum(vals)/len(vals)
        ax.plot([i-.22,i+.22],[mu,mu],color=CB["black"],linewidth=1.6,zorder=4)
    ax.set_xticks(range(len(gg)))
    ax.set_xticklabels([f"{g[0][:4]}\n{g[1][:5]}\nn={sum(1 for s in core if (meta[s]['location'],meta[s]['condition_std'])==g)}" for g in gg],fontsize=7.5)
    ax.set_title(titles[qi],fontsize=9)
axes[0].set_ylabel("Effective number of ASVs")
fig.suptitle("Figure 5. Hill numbers at depth 57,141. Points are samples; bars are group means.",
             fontsize=10,x=.02,ha="left")
fig.savefig(f"{F}/Figure5_hill_numbers.png"); plt.close(fig)
SRC["Figure 5"]="phase3/H4.log (iNEXT, per-sample qD at depth 57,141); metadata_v2_27-07-2026.tsv"

# ---------------- Figure 6: ordination ----------------
import random
rng=random.Random(42); DEPTH=57141
def rarefy(c,d):
    pool=[]
    for a,v in c.items(): pool.extend([a]*v)
    return Counter(pool if len(pool)<d else rng.sample(pool,d))
rar={s:rarefy(per_f[s],DEPTH) for s in core}
def bray(a,b):
    ks=set(a)|set(b); num=sum(abs(a.get(k,0)-b.get(k,0)) for k in ks); den=sum(a.values())+sum(b.values())
    return num/den if den else 0.
def jac(a,b):
    A,Bs=set(a),set(b); return 1-len(A&Bs)/len(A|Bs) if (A|Bs) else 0.
def pcoa(D,n):
    A=[[-0.5*D[i][j]**2 for j in range(n)] for i in range(n)]
    rm=[sum(r)/n for r in A]; cm=[sum(A[i][j] for i in range(n))/n for j in range(n)]; gm=sum(rm)/n
    G=[[A[i][j]-rm[i]-cm[j]+gm for j in range(n)] for i in range(n)]
    def ax_(G):
        v=[rng.random() for _ in range(n)]
        for _ in range(600):
            w=[sum(G[i][j]*v[j] for j in range(n)) for i in range(n)]
            nn=math.sqrt(sum(x*x for x in w)) or 1; v=[x/nn for x in w]
        lam=sum(v[i]*sum(G[i][j]*v[j] for j in range(n)) for i in range(n)); return lam,v
    l1,v1=ax_(G); G2=[[G[i][j]-l1*v1[i]*v1[j] for j in range(n)] for i in range(n)]; l2,v2=ax_(G2)
    tot=sum(G[i][i] for i in range(n))
    return ([v1[i]*math.sqrt(abs(l1)) for i in range(n)],[v2[i]*math.sqrt(abs(l2)) for i in range(n)],
            100*l1/tot,100*l2/tot)
fig,axes=plt.subplots(1,2,figsize=(9.2,4.2))
for ax,(nm,fnc) in zip(axes,[("Bray-Curtis",bray),("Jaccard",jac)]):
    n=len(core); D=[[fnc(rar[core[i]],rar[core[j]]) for j in range(n)] for i in range(n)]
    x,y,p1,p2=pcoa(D,n)
    for i,s in enumerate(core):
        ax.scatter(x[i],y[i],s=52,facecolor="white",edgecolor=SITE_C[meta[s]["location"]],
                   marker=COND_M[meta[s]["condition_std"]],linewidth=1.4,zorder=3)
    ax.axhline(0,color="#EEEEEE",zorder=0); ax.axvline(0,color="#EEEEEE",zorder=0)
    ax.set_xlabel(f"PCo1 ({p1:.1f}%)"); ax.set_ylabel(f"PCo2 ({p2:.1f}%)"); ax.set_title(nm,fontsize=9)
leg=[Line2D([0],[0],marker="o",color="w",markerfacecolor="w",markeredgecolor=SITE_C[s],
            markeredgewidth=1.4,markersize=7,label=s) for s in ("Christiana","Pretoria")]
leg+=[Line2D([0],[0],marker=COND_M[c],color="w",markerfacecolor="w",markeredgecolor=CB["black"],
             markeredgewidth=1.2,markersize=7,label=c) for c in ("asymptomatic","symptomatic")]
axes[1].legend(handles=leg,frameon=False,fontsize=8,bbox_to_anchor=(1.02,1),loc="upper left")
fig.suptitle("Figure 6. Ordination of the 16 core samples. Colour = site, symbol = condition.",
             fontsize=10,x=.02,ha="left")
fig.savefig(f"{F}/Figure6_ordination.png"); plt.close(fig)
SRC["Figure 6"]="fungi_only/table_exp/feature-table.tsv rarefied to 57,141; metadata_v2_27-07-2026.tsv"

# ---------------- Figure 7: Curvibasidium ----------------
def genus_of(a):
    t=tax.get(a,"")
    if not t or pid.get(a,0)<95: return ""
    for part in t.split(";"):
        part=part.strip()
        if part.startswith("g__"):
            v=part[3:].strip()
            return "" if (v.lower() in PH or v.lower().endswith("_incertae_sedis")) else v
    return ""
curv={a for a in tax if genus_of(a)=="Curvibasidium"}
fig,ax=plt.subplots(figsize=(6.0,4.0))
for i,g in enumerate(gg):
    ss=[s for s in core if (meta[s]["location"],meta[s]["condition_std"])==g]
    vals=[]
    for s in ss:
        t=sum(per_f[s].values()); c=sum(v for a,v in per_f[s].items() if a in curv)
        vals.append(100*c/t if t else 0)
    xs=[i+(hash(s)%100-50)/420 for s in ss]
    ax.scatter(xs,vals,s=40,facecolor="white",edgecolor=SITE_C[g[0]],marker=COND_M[g[1]],
               linewidth=1.3,zorder=3)
    sv=sorted(vals); med=(sv[len(sv)//2] if len(sv)%2 else (sv[len(sv)//2-1]+sv[len(sv)//2])/2)
    ax.plot([i-.22,i+.22],[med,med],color=CB["black"],linewidth=1.8,zorder=4)
    ax.text(i,-1.6,f"median {med:.2f}%",ha="center",fontsize=7.2)
ax.set_xticks(range(len(gg)))
ax.set_xticklabels([f"{g[0]}\n{g[1]}\n(n={sum(1 for s in core if (meta[s]['location'],meta[s]['condition_std'])==g)})" for g in gg],fontsize=8)
ax.set_ylabel("Curvibasidium (% of fungal reads)"); ax.set_ylim(-3,22)
ax.set_title("Figure 7. Per-sample Curvibasidium relative abundance",fontsize=10,loc="left")
fig.savefig(f"{F}/Figure7_curvibasidium.png"); plt.close(fig)
SRC["Figure 7"]="fungi_only/table_exp/feature-table.tsv; taxonomy_unite10 (genus at >=95% identity)"

# ---------------- Supplementary ----------------
# S2 DADA2 retention
stages=["input","filtered","denoised","merged","non-chimeric"]
acc=Counter()
for p in ("groupA2","groupB","otp"):
    fp=f"{B}/phase2/{p}/stats_exp/stats.tsv"
    if not os.path.exists(fp): continue
    with open(fp) as fh:
        rd=csv.reader(fh,delimiter="\t"); hh=next(rd); next(rd)
        idx={c:hh.index(c) for c in stages if c in hh}
        for r in rd:
            if r and r[0].strip():
                for c in stages: acc[c]+=int(float(r[idx[c]]))
fig,ax=plt.subplots(figsize=(5.4,3.4))
ax.plot(range(len(stages)),[acc[c] for c in stages],marker="o",color=CB["blue"])
for i,c in enumerate(stages):
    ax.annotate(f"{acc[c]:,}",(i,acc[c]),textcoords="offset points",xytext=(0,7),ha="center",fontsize=7.5)
ax.set_xticks(range(len(stages))); ax.set_xticklabels(stages,rotation=20)
ax.set_ylabel("Reads (all 45 samples)"); ax.set_title("Figure S2. DADA2 read retention by stage",fontsize=10,loc="left")
fig.savefig(f"{F}/FigureS2_dada2_retention.png"); plt.close(fig)
SRC["Figure S2"]="phase2/{groupA2,groupB,otp}/stats_exp/stats.tsv"

# S3 taxonomy assignment rate by rank
def rate(tp,thr=None):
    d={}
    for r in csv.reader(open(tp),delimiter="\t"):
        if r and r[0]!="Feature ID" and not r[0].startswith("#"): d[r[0]]=r[1]
    out=[]
    for _,rk in RANKS:
        n=0
        for a,t in d.items():
            rr=ranks_of(t)
            if thr is None:
                if rr[rk]: n+=1
            else:
                if rr[rk] and pid.get(a,0)>=CUT[rk]: n+=1
        out.append(100*n/len(d))
    return out
r82=rate(f"{O}/../outputs_27-07-2026/taxonomy_unite10/exp/taxonomy.tsv") if False else None
series={"UNITE 8.2":rate(f"{B}/phase2/taxonomy/exp/taxonomy.tsv"),
        "UNITE 10.0 permissive":rate(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"),
        "UNITE 10.0 rank-threshold":rate(f"{O}/taxonomy_unite10/exp/taxonomy.tsv",thr=True)}
fig,ax=plt.subplots(figsize=(6.2,3.8))
for (k,v),c in zip(series.items(),[CB["grey"],CB["orange"],CB["blue"]]):
    ax.plot(range(len(RANKS)),v,marker="o",label=k,color=c)
ax.set_xticks(range(len(RANKS))); ax.set_xticklabels([r for _,r in RANKS],rotation=25)
ax.set_ylabel("ASVs assigned (%)"); ax.legend(frameon=False,fontsize=8)
ax.set_title("Figure S3. Taxonomic assignment rate by rank",fontsize=10,loc="left")
fig.savefig(f"{F}/FigureS3_assignment_by_rank.png"); plt.close(fig)
SRC["Figure S3"]="phase2/taxonomy/exp/taxonomy.tsv; taxonomy_unite10/exp/taxonomy.tsv; search.qza"

# S4 guild coverage
gd=list(csv.DictReader(open(f"{O}/guilds/guild_assignments_per_ASV.tsv"),delimiter="\t"))
grand=sum(int(r["total_reads"]) for r in gd)
ft_a=sum(1 for r in gd if r["FungalTraits_primary_lifestyle"]); ft_r=sum(int(r["total_reads"]) for r in gd if r["FungalTraits_primary_lifestyle"])
fg_a=sum(1 for r in gd if r["FUNGuild_guild"]); fg_r=sum(int(r["total_reads"]) for r in gd if r["FUNGuild_guild"])
fg_r_corr=sum(int(r["total_reads"]) for r in gd if r["FUNGuild_guild"] and r["ASV_ID_md5"] not in yeast)
fig,ax=plt.subplots(figsize=(6.0,3.6))
labels=["FungalTraits\n% ASVs","FungalTraits\n% reads","FUNGuild\n% ASVs",
        "FUNGuild\n% reads\n(as reported)","FUNGuild\n% reads\n(excl. misassigned\nTremellomycetes lineage)"]
vals=[100*ft_a/len(gd),100*ft_r/grand,100*fg_a/len(gd),100*fg_r/grand,100*fg_r_corr/grand]
cols_=[CB["green"],CB["green"],CB["purple"],CB["purple"],CB["blue"]]
ax.bar(range(5),vals,color=cols_,width=.6)
for i,v in enumerate(vals): ax.text(i,v+1.5,f"{v:.1f}",ha="center",fontsize=8)
ax.set_xticks(range(5)); ax.set_xticklabels(labels,fontsize=7)
ax.set_ylabel("Coverage (%)"); ax.set_ylim(0,100)
ax.set_title("Figure S4. Guild assignment coverage by tool",fontsize=10,loc="left")
fig.savefig(f"{F}/FigureS4_guild_coverage.png"); plt.close(fig)
SRC["Figure S4"]="guilds/guild_assignments_per_ASV.tsv"

print("WROTE:")
for f in sorted(os.listdir(F)):
    if f.endswith(".png"): print(f"  {F}/{f}")
print("\nSOURCE FILES PER FIGURE:")
for k in sorted(SRC): print(f"  {k}: {SRC[k]}")
