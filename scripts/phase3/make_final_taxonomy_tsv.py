#!/usr/bin/env python
"""Write the FINAL reported taxonomy (rank-specific identity thresholds applied to the
permissive UNITE 10.0 run, yeast lineage relabelled) as a QIIME-style taxonomy TSV,
so the ASV catalogue is built on exactly the taxonomy the manuscript reports."""
import csv, zipfile
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
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
fungi=set()
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    if line.startswith(">"): fungi.add(line[1:].split()[0])

out=f"{O}/final_taxonomy.tsv"
n_yeast=0
with open(out,"w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["Feature ID","Taxon","Confidence"])
    for a in sorted(fungi):
        if a in yeast:
            # placement-based label; the erroneous Agaricales chain is NOT propagated
            w.writerow([a,"k__Fungi;p__Basidiomycota;c__Tremellomycetes","placement-based (ITS2 unresolved below class)"])
            n_yeast+=1; continue
        rk=ranks_of(tax.get(a,"")); p=pid.get(a,0.); keep=[]
        for pre,r in RANKS:
            if rk[r] and p>=CUT[r]: keep.append(f"{pre}{rk[r]}")
            else: break
        w.writerow([a,";".join(keep) if keep else "Unassigned",f"{p:.1f}" if p else ""])
print(f"wrote {out}")
print(f"  ASVs: {len(fungi)}   yeast lineage written as Tremellomycetes (class): {n_yeast}")
