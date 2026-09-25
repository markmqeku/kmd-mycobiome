#!/usr/bin/env python
"""T-4d: host (Searsia lancea) read fraction per sample, and by tissue / site / condition_std."""
import csv, statistics as st
from collections import defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
D = f"{B}/phase3/outputs_27-07-2026/nontarget_screen"

# confirmed host ASVs: >=95% identity AND >=80% query coverage vs Anacardiaceae/Searsia reference
best = {}
for r in csv.reader(open(f"{D}/blast_host.tsv"), delimiter="\t"):
    if len(r) < 8: continue
    q, pid, qcov = r[0], float(r[2]), float(r[4])
    if q not in best or pid > best[q][0]:
        best[q] = (pid, qcov)
host = {q for q, (pid, qcov) in best.items() if pid >= 95 and qcov >= 80}
with open(f"{D}/host_asv_ids.txt", "w") as fh:
    fh.write("\n".join(sorted(host)) + "\n")

lines = [l.rstrip("\n") for l in open(f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); samples = [h for h in hdr[1:] if h.strip()]
tot = defaultdict(int); hst = defaultdict(int)
for l in lines[1:]:
    p = l.split("\t"); a = p[0]
    vals = [int(float(x)) for x in p[1:1+len(samples)]]
    for s, v in zip(samples, vals):
        tot[s] += v
        if a in host: hst[s] += v

meta = {}
with open(f"{B}/phase3/outputs_27-07-2026/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

print(f"confirmed host ASVs: {len(host)}   host reads: {sum(hst.values()):,} / {sum(tot.values()):,} "
      f"= {100*sum(hst.values())/sum(tot.values()):.1f}% of dataset\n")
print(f"{'sample':7} {'site':13} {'tissue':14} {'cond_std':14} {'total':>9} {'host':>9} {'host%':>7} {'fungal_after':>13}")
rows = []
for s in sorted(samples, key=lambda x: -(hst[x]/tot[x] if tot[x] else 0)):
    m = meta.get(s, {})
    pct = 100*hst[s]/tot[s] if tot[s] else 0
    rows.append((s, m.get("tissue",""), m.get("location",""), m.get("condition_std",""), tot[s], hst[s], pct))
    print(f"{s:7} {m.get('location',''):13} {m.get('tissue',''):14} {m.get('condition_std',''):14} "
          f"{tot[s]:>9,} {hst[s]:>9,} {pct:>6.1f}% {tot[s]-hst[s]:>13,}")

def grp(idx, label):
    d = defaultdict(lambda: [0, 0])
    for r in rows:
        k = r[idx]; d[k][0] += r[4]; d[k][1] += r[5]
    print(f"\n--- host load by {label} ---")
    for k in sorted(d):
        t, hh = d[k]
        print(f"  {k:16} total={t:>10,}  host={hh:>10,}  {100*hh/t if t else 0:5.1f}%")
grp(1, "tissue"); grp(2, "site"); grp(3, "condition_std")

below = [(r[0], r[4]-r[5]) for r in rows if (r[4]-r[5]) < 5000]
print(f"\nsamples that would fall BELOW depth 5,000 after host removal: {sorted(below, key=lambda x:x[1])}")
