#!/usr/bin/env python
"""Item 7 prep: write the 28 target ASVs grouped by neighbourhood, and a FASTA of the
4 that need an NCBI triage first (3 Ascomycota-only, 1 fully unassigned)."""
import csv, zipfile, os
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
P = f"{O}/placement"; os.makedirs(P, exist_ok=True)
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUT = {"species":97.,"genus":95.,"family":90.,"order":85.,"class":80.,"phylum":80.,"kingdom":80.}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

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
yeast={r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"),delimiter="\t")
       if r and r[0]!="ASV_ID_md5" and r[2]=="1"}

final={}
for a in tax:
    rk=ranks_of(tax[a]); p=pid.get(a,0.); keep={}
    for _,r in RANKS:
        if rk[r] and p>=CUT[r]: keep[r]=rk[r]
        else: break
    final[a]=keep

lines=[l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
       if l.strip() and not l.startswith("# Constructed")]
hdr=lines[0].lstrip("#").split("\t"); cols=[h for h in hdr[1:] if h.strip()]
per={s:{} for s in cols}
for l in lines[1:]:
    q=l.split("\t")
    for s,v in zip(cols,q[1:1+len(cols)]):
        v=int(float(v))
        if v: per[s][q[0]]=v
meta={}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd=csv.reader(fh,delimiter="\t"); h=next(rd); next(rd); di={c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]]={c:r[di[c]] for c in h}
groups=defaultdict(list)
for s in cols:
    if meta.get(s): groups[(meta[s]["location"],meta[s]["condition_std"])].append(s)
maxrel=defaultdict(float)
for g,ss in groups.items():
    tot=sum(per[s][a] for s in ss for a in per[s])
    if not tot: continue
    acc=Counter()
    for s in ss:
        for a,v in per[s].items(): acc[a]+=v
    for a,v in acc.items():
        rel=100*v/tot
        if rel>maxrel[a]: maxrel[a]=rel

target=[a for a in final if not final[a].get("genus") and maxrel.get(a,0)>=0.5 and a not in yeast]
seqs={}; sid=None; buf=[]
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line=line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid]="".join(buf)
        sid,buf=line[1:].split()[0],[]
    else: buf.append(line.strip())
if sid: seqs[sid]="".join(buf)

def deepest(d):
    last=(None,None)
    for _,r in RANKS:
        if d.get(r): last=(r,d[r])
    return last

nb=defaultdict(list)
for a in target:
    rk,val = deepest(final[a])
    nb[val if val else "UNASSIGNED"].append(a)

print(f"target ASVs to place: {len(target)}")
with open(f"{P}/targets.tsv","w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["ASV_ID_md5","neighbourhood","rank","max_group_rel_pct","length_bp"])
    for k in sorted(nb, key=lambda x:-len(nb[x])):
        for a in nb[k]:
            rk,_=deepest(final[a])
            w.writerow([a,k,rk or "none",f"{maxrel[a]:.2f}",len(seqs.get(a,''))])
        print(f"  {k:22} {len(nb[k])} ASV(s)")

# neighbourhood FASTAs
for k,al in nb.items():
    if k in ("Ascomycota","UNASSIGNED"): continue
    with open(f"{P}/{k}_asvs.fasta","w") as fh:
        for a in al: fh.write(f">ASV_{a[:8]}_rel{maxrel[a]:.1f}\n{seqs[a]}\n")
# triage set
triage=[a for k in ("Ascomycota","UNASSIGNED") for a in nb.get(k,[])]
with open(f"{P}/triage_query.fasta","w") as fh:
    for a in triage: fh.write(f">{a}\n{seqs[a]}\n")
print(f"\ntriage (NCBI first): {len(triage)} ASVs -> {P}/triage_query.fasta")
print(f"neighbourhood FASTAs written to {P}")
