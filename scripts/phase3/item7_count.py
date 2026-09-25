#!/usr/bin/env python
"""Item 7 scope: ASVs unassigned at GENUS level with >=0.5% relative abundance in ANY
site-by-condition group. Count first; if >40, stop and report."""
import csv, zipfile
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUT = {"species":97.0,"genus":95.0,"family":90.0,"order":85.0,"class":80.0,"phylum":80.0,"kingdom":80.0}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

def ranks_of(t):
    out = {r:"" for _,r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part = part.strip()
        for pre,rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk=="species" and (b.endswith("_sp") or b.endswith("_sp.")))): v = ""
                out[rk] = v
    return out

tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}

final = {}
for a in tax:
    rk = ranks_of(tax[a]); p = pid.get(a, 0.0); keep = {}
    for _, r in RANKS:
        if rk[r] and p >= CUT[r]: keep[r] = rk[r]
        else: break
    final[a] = keep

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

groups = defaultdict(list)
for s in cols:
    if meta.get(s): groups[(meta[s]["location"], meta[s]["condition_std"])].append(s)

maxrel = defaultdict(float); grp_of = {}
for g, ss in groups.items():
    tot = sum(per[s][a] for s in ss for a in per[s])
    if not tot: continue
    acc = Counter()
    for s in ss:
        for a, v in per[s].items(): acc[a] += v
    for a, v in acc.items():
        rel = 100*v/tot
        if rel > maxrel[a]: maxrel[a] = rel; grp_of[a] = g

target = [a for a in final if not final[a].get("genus") and maxrel.get(a, 0) >= 0.5]
print(f"ASVs unassigned at genus level with >=0.5% in any site x condition group: {len(target)}")
print(f"   of these, already-placed Tremellomycetes lineage members: {sum(1 for a in target if a in yeast)}")
print(f"   REMAINING to place: {len([a for a in target if a not in yeast])}")
print(f"\n{'ASV':10} {'max rel%':>9} {'in group':28} {'deepest confident assignment':40}")
for a in sorted(target, key=lambda x: -maxrel[x]):
    d = final[a]; last = ("none","")
    for _, r in RANKS:
        if d.get(r): last = (r, d[r])
    tagy = "  [Tremellomycetes lineage, already placed]" if a in yeast else ""
    print(f"{a[:9]:10} {maxrel[a]:8.2f}% {str(grp_of[a]):28} {last[1]+' ('+last[0]+')':40}{tagy}")
print("\n--- grouping for per-neighbourhood trees (excluding the already-placed lineage) ---")
c = Counter()
for a in target:
    if a in yeast: continue
    d = final[a]; last = ("none","")
    for _, r in RANKS:
        if d.get(r): last = (r, d[r])
    c[f"{last[1]} ({last[0]})"] += 1
for k, v in c.most_common(): print(f"   {k:44} {v} ASV(s)")
