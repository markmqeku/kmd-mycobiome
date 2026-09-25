#!/usr/bin/env python
"""G-5 corrected: the 78 yeast ASVs are relabelled from UNITE's erroneous 'Agaricales' to
'Tremellomycetes yeast lineage (ITS2-unresolved)' on the phylogenetic + identity evidence from G-3.
Rank-truncation cannot fix a WRONG-LINEAGE match, only an over-deep one; this correction is evidence-based."""
import csv, zipfile
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUT = {"species":97.0,"genus":95.0,"family":90.0,"order":85.0,"class":80.0,"phylum":80.0,"kingdom":80.0}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}
YEAST_LABEL = "Tremellomycetes yeast lineage (ITS2-unresolved)"

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

tax80 = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax80[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
f = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(f).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v

yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
fungi = set()
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    if line.startswith(">"): fungi.add(line[1:].split()[0])

def label(a):
    if a in yeast: return YEAST_LABEL
    rk = ranks_of(tax80.get(a,"")); p = pid.get(a, 0.0); keep = {}
    for _, r in RANKS:
        if rk[r] and p >= CUT[r]: keep[r] = rk[r]
        else: break
    last = (None,None)
    for _, r in RANKS:
        if keep.get(r): last = (r, keep[r])
    return f"{last[1]} ({last[0]})" if last[0] else "unassigned"

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

lab = {a: label(a) for a in fungi}
print("=== G-5 CORRECTED COMPOSITION (fungi-only table; yeast lineage relabelled on G-3 evidence) ===")
print("ALL earlier composition figures (Didymella 64->10%, Curvibasidium, Mycosphaerella) are OBSOLETE.\n")
groups = defaultdict(Counter); gtot = Counter(); gn = Counter()
for s in cols:
    m = meta.get(s)
    if not m: continue
    k = (m["location"], m["condition_std"]); gn[k] += 1
    for a, v in per[s].items():
        groups[k][lab[a]] += v; gtot[k] += v
for k in sorted(groups):
    print(f"--- {k[0]} / {k[1]}   ({gtot[k]:,} reads, n={gn[k]}) ---")
    for l2, v in groups[k].most_common(10):
        print(f"     {l2:48} {100*v/gtot[k]:5.1f}%")
    print()

print("=== the yeast lineage, per site ===")
for site in ("Bloemfontein","Christiana","Pretoria"):
    ss = [s for s in cols if meta.get(s) and meta[s]["location"] == site]
    t = sum(sum(per[s].values()) for s in ss)
    y = sum(v for s in ss for a, v in per[s].items() if a in yeast)
    print(f"   {site:14} {y:>9,} / {t:>9,} fungal reads = {100*y/t if t else 0:5.1f}%")
