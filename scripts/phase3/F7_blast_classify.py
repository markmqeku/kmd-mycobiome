#!/usr/bin/env python
"""F-7: classify remote-BLAST hits for the unknown ASVs; plus F-1 site distribution."""
import csv, re
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
D = f"{B}/phase3/outputs_27-07-2026/nontarget_screen"

rows = [r for r in csv.reader(open(f"{D}/query_abundance.tsv"), delimiter="\t") if r and r[0] != "ASV_ID_md5"]
ab = {r[0]: int(r[1]) for r in rows}
unknown = {r[0] for r in rows if r[2] == "1"}

best = {}
try:
    for r in csv.reader(open(f"{D}/blast_priority.tsv"), delimiter="\t"):
        if len(r) < 8: continue
        q, pid, qcov, ev, title = r[0], float(r[2]), float(r[4]), r[5], r[7]
        if q not in best or pid > best[q][0]:
            best[q] = (pid, qcov, ev, title)
except FileNotFoundError:
    print("no remote BLAST output yet"); raise SystemExit

PLANT = ("Searsia","Rhus","Anacard","chloroplast","Arabidopsis","Oryza","Vitis","Populus","Quercus",
         "Solanum","Nicotiana","Citrus","Mangifera","Pistacia","Schinus","Toxicodendron","plastid","Viridiplantae")
FUNGI = ("fungal","fungus","Ascomycota","Basidiomycota","Mucor","Fusarium","Penicillium","Aspergillus",
         "Cladosporium","yeast","mycor","Chytrid","Rhizophydi","Spizellomyc","Mortierell","uncultured fungus",
         "Agaric","Homophron","Psathyrella","Orbilia","Didymella","Alternaria","Aureobasidium")
BACT  = ("bacteri","Pseudomonas","Bacillus","Streptomyces","16S ribosomal RNA","Escherichia")

def klass(t):
    tl = t.lower()
    if any(k.lower() in tl for k in PLANT): return "Viridiplantae/plant"
    if any(k.lower() in tl for k in FUNGI): return "Fungi"
    if any(k.lower() in tl for k in BACT):  return "Bacteria"
    if "uncultured eukaryote" in tl or "environmental" in tl: return "other/uncertain eukaryote"
    return "other/uncertain"

print(f"=== F-7: remote NCBI nt BLAST — the {len(unknown)} unknown ('Homophron') ASVs ===")
hit = [a for a in unknown if a in best]
print(f"unknown ASVs with a remote hit: {len(hit)}/{len(unknown)}   "
      f"(reads represented: {sum(ab.get(a,0) for a in hit):,})")
cls = Counter(); reads = Counter()
for a in hit:
    k = klass(best[a][3]); cls[k] += 1; reads[k] += ab.get(a, 0)
nohit = [a for a in unknown if a not in best]
cls["no confident hit"] += len(nohit); reads["no confident hit"] += sum(ab.get(a,0) for a in nohit)
print(f"\n{'class':28} {'ASVs':>6} {'reads':>12}")
for k, v in cls.most_common():
    print(f"{k:28} {v:>6} {reads[k]:>12,}")

print("\ntop unknown ASVs by reads (identity / coverage / e-value / top hit):")
for a in sorted(hit, key=lambda x: -ab.get(x, 0))[:15]:
    pid, qcov, ev, title = best[a]
    print(f"  {a[:9]} reads={ab.get(a,0):>8,} pid={pid:5.1f} qcov={qcov:3.0f} e={ev:<9} {klass(title):22} {title[:52]}")

# ---- F-1 site distribution of unknowns ----
lines = [l.rstrip("\n") for l in open(f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); samples = [h for h in hdr[1:] if h.strip()]
counts = {}
for l in lines[1:]:
    p = l.split("\t"); counts[p[0]] = [int(float(x)) for x in p[1:1+len(samples)]]
meta = {}
with open(f"{B}/phase3/outputs_27-07-2026/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
print("\n=== F-1: the 78 unknown ASVs, per site ===")
print(f"{'site':14} {'total reads':>12} {'unknown reads':>14} {'%':>8} {'samples containing':>20}")
for loc in ("Bloemfontein", "Christiana", "Pretoria"):
    idx = [j for j, s in enumerate(samples) if meta[s]["location"] == loc]
    t = sum(counts[a][j] for a in counts for j in idx)
    u = sum(counts[a][j] for a in unknown if a in counts for j in idx)
    pos = sum(1 for j in idx if sum(counts[a][j] for a in unknown if a in counts) > 0)
    print(f"{loc:14} {t:>12,} {u:>14,} {100*u/t if t else 0:>7.2f}% {pos:>12}/{len(idx)}")
