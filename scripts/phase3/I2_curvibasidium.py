#!/usr/bin/env python
"""I-2: per-sample relative abundance of Curvibasidium across all 16 core samples.
Determines whether the enrichment is consistent within groups or driven by one or two samples."""
import csv, zipfile
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
f = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(f).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v

def genus(a):
    t = tax.get(a, "")
    if not t or t.lower().startswith("unassigned"): return ""
    for part in t.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            v = part[3:].strip()
            if v.lower() in PH or v.lower().endswith("_incertae_sedis"): return ""
            return v
    return ""

curv = {a for a in tax if genus(a) == "Curvibasidium" and pid.get(a, 0) >= 95}
print(f"ASVs assigned to Curvibasidium at >=95% identity: {len(curv)}")

lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
per = {s: {} for s in cols}
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]

print(f"\n{'sample':7} {'site':11} {'condition':13} {'tissue':14} {'total':>9} {'Curvi':>8} {'rel%':>7} {'nASV':>5}")
rows = []
for site in ("Christiana","Pretoria"):
    for cond in ("asymptomatic","symptomatic"):
        ss = [s for s in core if meta[s]["location"]==site and meta[s]["condition_std"]==cond]
        for s in sorted(ss):
            tot = sum(per[s].values())
            cv = sum(v for a,v in per[s].items() if a in curv)
            n = sum(1 for a in per[s] if a in curv)
            rel = 100*cv/tot if tot else 0
            rows.append((site,cond,s,rel,cv,tot,n))
            print(f"{s:7} {site:11} {cond:13} {meta[s]['tissue']:14} {tot:>9,} {cv:>8,} {rel:>6.2f}% {n:>5}")
        print()

print("=== consistency within each group ===")
for site in ("Christiana","Pretoria"):
    for cond in ("asymptomatic","symptomatic"):
        v = [r[3] for r in rows if r[0]==site and r[1]==cond]
        if not v: continue
        nz = sum(1 for x in v if x > 0)
        above1 = sum(1 for x in v if x >= 1.0)
        print(f"  {site:11} {cond:13} n={len(v)}  values=[{', '.join(f'{x:.2f}' for x in sorted(v, reverse=True))}]")
        print(f"      detected in {nz}/{len(v)}; >=1% in {above1}/{len(v)}; "
              f"min={min(v):.2f}% max={max(v):.2f}% median={sorted(v)[len(v)//2]:.2f}%")
print("\n=== verdict ===")
for site in ("Christiana","Pretoria"):
    a = [r[3] for r in rows if r[0]==site and r[1]=="asymptomatic"]
    s = [r[3] for r in rows if r[0]==site and r[1]=="symptomatic"]
    overlap = max(a) >= min(s)
    print(f"  {site}: asymptomatic max={max(a):.2f}%  symptomatic min={min(s):.2f}%  "
          f"{'RANGES OVERLAP' if overlap else 'ranges separate'}")
