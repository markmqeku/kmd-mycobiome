#!/usr/bin/env python
"""T-4c/T-7 preview: design retention after removing confirmed host ASVs, at several depths."""
import csv
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
D = f"{B}/phase3/outputs_27-07-2026/nontarget_screen"
host = {l.strip() for l in open(f"{D}/host_asv_ids.txt") if l.strip()}

lines = [l.rstrip("\n") for l in open(f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); samples = [h for h in hdr[1:] if h.strip()]
fung = defaultdict(int); nasv = defaultdict(int); kept_asv = set()
for l in lines[1:]:
    p = l.split("\t"); a = p[0]
    if a in host: continue
    vals = [int(float(x)) for x in p[1:1+len(samples)]]
    if sum(vals) > 0: kept_asv.add(a)
    for s, v in zip(samples, vals):
        fung[s] += v
        if v > 0: nasv[s] += 1

meta = {}
with open(f"{B}/phase3/outputs_27-07-2026/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

print(f"ASVs after host removal: {len(kept_asv)} (was 1322; removed {len(host)} host ASVs)")
print(f"fungal reads total     : {sum(fung.values()):,} (was 4,010,617)\n")

for DEPTH in (5000, 2000, 1000, 500):
    kept = [s for s in samples if fung[s] >= DEPTH]
    print(f"=== depth {DEPTH}: retained {len(kept)}/45 ===")
    b = Counter((meta[s]["tissue"], meta[s]["age"]) for s in kept if meta[s]["location"] == "Bloemfontein")
    core = {k: v for k, v in b.items() if k[0] in ("Leaves", "Twigs")}
    print(f"   ARTICLE 1 Bloemfontein retained: {sum(b.values())}/29   core(Leaves+Twigs) cells: "
          f"{dict(sorted(core.items())) if core else 'NONE'}")
    for site in ("Christiana", "Pretoria"):
        c = Counter((meta[s]["condition_std"], meta[s]["tissue"]) for s in kept if meta[s]["location"] == site)
        print(f"   ARTICLE 2 {site:11}: {sum(c.values())}/{sum(1 for s in samples if meta[s]['location']==site)} retained")
    print()

print("Bloemfontein samples ranked by FUNGAL reads (post-host):")
bl = sorted([s for s in samples if meta[s]["location"] == "Bloemfontein"], key=lambda s: -fung[s])
for s in bl:
    print(f"   {s:5} {meta[s]['tissue']:14} {meta[s]['age']:6} fungal={fung[s]:>8,}  ASVs={nasv[s]:>4}")
