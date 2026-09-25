#!/usr/bin/env python
"""F-1: per-site distribution of the 78 unknown ('Homophron') ASVs.
   F-5: Christiana/Pretoria retention after host removal at 1000/2500/5000/10000, per cell.
   F-6: per-Bloemfontein-sample totals."""
import csv
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
D = f"{B}/phase3/outputs_27-07-2026/nontarget_screen"

host = {l.strip() for l in open(f"{D}/host_asv_ids_ALL.txt") if l.strip()}
unknown = {r[0] for r in csv.reader(open(f"{D}/query_abundance.tsv"), delimiter="\t")
           if r and r[0] != "ASV_ID_md5" and r[2] == "1"}

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

# ---------- F-1 ----------
print("=== F-1: the 78 unknown ASVs, per site ===")
site_tot = defaultdict(int); site_unk = defaultdict(int); site_pos = defaultdict(int); site_n = defaultdict(int)
for j, s in enumerate(samples):
    loc = meta[s]["location"]; site_n[loc] += 1
    t = sum(counts[a][j] for a in counts)
    u = sum(counts[a][j] for a in unknown if a in counts)
    site_tot[loc] += t; site_unk[loc] += u
    if u > 0: site_pos[loc] += 1
print(f"{'site':14} {'total reads':>12} {'unknown reads':>14} {'%':>7} {'samples with them':>19}")
for loc in ("Bloemfontein", "Christiana", "Pretoria"):
    print(f"{loc:14} {site_tot[loc]:>12,} {site_unk[loc]:>14,} "
          f"{100*site_unk[loc]/site_tot[loc] if site_tot[loc] else 0:>6.2f}% {site_pos[loc]:>10}/{site_n[loc]}")

# ---------- F-6 ----------
print("\n=== F-6: per-Bloemfontein-sample ===")
print(f"{'sample':7} {'tissue':14} {'age':6} {'total':>9} {'host':>9} {'fungal':>9} {'fungalASVs':>11}")
bl = [s for s in samples if meta[s]["location"] == "Bloemfontein"]
for s in sorted(bl, key=lambda x: -sum(counts[a][samples.index(x)] for a in counts if a not in host)):
    j = samples.index(s)
    t = sum(counts[a][j] for a in counts)
    hh = sum(counts[a][j] for a in host if a in counts)
    f = t - hh
    na = sum(1 for a in counts if a not in host and counts[a][j] > 0)
    print(f"{s:7} {meta[s]['tissue']:14} {meta[s]['age']:6} {t:>9,} {hh:>9,} {f:>9,} {na:>11}")

# ---------- F-5 ----------
print("\n=== F-5: Christiana / Pretoria retention after host removal ===")
fung = {}
for j, s in enumerate(samples):
    fung[s] = sum(counts[a][j] for a in counts if a not in host)
for depth in (1000, 2500, 5000, 10000):
    print(f"\n--- depth {depth} ---")
    for site in ("Christiana", "Pretoria"):
        sub = [s for s in samples if meta[s]["location"] == site]
        kept = [s for s in sub if fung[s] >= depth]
        cells = Counter((meta[s]["condition_std"], meta[s]["tissue"]) for s in kept)
        print(f"   {site:11} retained {len(kept)}/{len(sub)}")
        for k in sorted(cells):
            print(f"      {k[0]:13} {k[1]:14} n={cells[k]}{' *' if cells[k] < 3 else ''}")
        lost = [s for s in sub if s not in kept]
        if lost: print(f"      dropped: {[(s, fung[s]) for s in lost]}")
