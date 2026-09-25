#!/usr/bin/env python
"""13a supplementary. The NCBI screen found a Searsia lancea ASV still in the fungal table.
It was not a reference-set gap: the original host screen DID hit it, at 93.2% identity, and
the >=95% identity cut discarded the hit. Quantify how many other ASVs sit in that same
band, using the host BLAST output already on disk. No network needed."""
import csv
from collections import defaultdict
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
best = {}
for r in csv.reader(open(f"{O}/nontarget_screen/blast_host_ALL.tsv"), delimiter="\t"):
    if len(r) < 4: continue
    try: pid, qcov = float(r[2]), float(r[3])
    except ValueError: continue
    if r[0] not in best or pid > best[r[0]][0]: best[r[0]] = (pid, qcov)

lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
reads, samples = {}, {}
for l in lines[1:]:
    q = l.split("\t"); vals = [int(float(v)) for v in q[1:1+len(cols)]]
    reads[q[0]] = sum(vals); samples[q[0]] = [s for s, v in zip(cols, vals) if v]
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
core = {s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")}

print("ASVs RETAINED in the fungal table that hit a Searsia/Rhus reference,")
print("binned by the identity the original screen measured (cut was >=95% id and >=80% qcov):\n")
print(f"{'identity band':16} {'ASVs':>5} {'reads':>9} {'in core-16':>11}")
bands = [(95,101,">=95 (removed as host)"),(93,95,"93.0-94.9"),(90,93,"90.0-92.9"),
         (85,90,"85.0-89.9"),(0,85,"<85")]
rows = []
for lo, hi, name in bands:
    grp = [a for a in reads if a in best and lo <= best[a][0] < hi and best[a][1] >= 80]
    nc = sum(1 for a in grp if any(s in core for s in samples[a]))
    print(f"{name:16} {len(grp):>5} {sum(reads[a] for a in grp):>9,} {nc:>11}")
    rows.append((name, grp))
print("\nThe >=95 row should be empty: those ASVs were removed as host before this table was built.")
print("\nindividual ASVs in the 90-95% band (the band that produced the missed Searsia ASV):")
print(f"{'asv':14} {'%id':>6} {'qcov':>5} {'reads':>8} {'n_samp':>7} sites")
for name, grp in rows:
    if not name.startswith(("93.0", "90.0")): continue
    for a in sorted(grp, key=lambda x: -reads[x]):
        sites = sorted({meta.get(s, {}).get("location", "?") for s in samples[a]})
        print(f"{a[:12]:14} {best[a][0]:>6.1f} {best[a][1]:>5.0f} {reads[a]:>8,} "
              f"{len(samples[a]):>7} {','.join(sites)}")
