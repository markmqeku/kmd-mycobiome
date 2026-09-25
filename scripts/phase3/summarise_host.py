import csv
D = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen"
best = {}
for r in csv.reader(open(f"{D}/blast_host.tsv"), delimiter="\t"):
    if len(r) < 8: continue
    q, pid, qcov, title = r[0], float(r[2]), float(r[4]), r[7]
    if q not in best or pid > best[q][0]:
        best[q] = (pid, qcov, title)

rows = [r for r in csv.reader(open(f"{D}/query_abundance.tsv"), delimiter="\t") if r and r[0] != "ASV_ID_md5"]
ab = {r[0]: int(r[1]) for r in rows}
hom = {r[0] for r in rows if r[2] == "1"}
qtot = sum(ab.values())

hit_reads = sum(ab.get(q, 0) for q in best)
print(f"\n=== HOST-FAMILY (Anacardiaceae/Searsia/Rhus) BLAST OF NON-TARGET CANDIDATES ===")
print(f"query ASVs                 : {len(ab)}   ({qtot:,} reads)")
print(f"ASVs with any host hit     : {len(best)}   ({hit_reads:,} reads = {100*hit_reads/qtot:.1f}% of query set)")
strong = {q: v for q, v in best.items() if v[0] >= 95 and v[1] >= 80}
sreads = sum(ab.get(q, 0) for q in strong)
print(f"STRONG hits (>=95% id, >=80% qcov): {len(strong)} ASVs, {sreads:,} reads ({100*sreads/qtot:.1f}% of query set)")
hs = [q for q in strong if q in hom]
print(f"  of the 78 Homophron ASVs, strong host hits: {len(hs)}")
print(f"  of the 78 Homophron ASVs, any host hit    : {len([q for q in best if q in hom])}")

print("\n  top host-matching ASVs by reads:")
for q, v in sorted(best.items(), key=lambda x: -ab.get(x[0], 0))[:12]:
    tag = "HOMOPHRON" if q in hom else ""
    print(f"    {q[:10]} reads={ab.get(q,0):>9,} pid={v[0]:5.1f} qcov={v[1]:3.0f}  {v[2][:60]} {tag}")

import collections
org = collections.Counter()
for q, v in strong.items():
    t = v[2]
    parts = t.split()
    name = " ".join(parts[1:3]) if len(parts) > 2 else t[:30]
    org[name] += ab.get(q, 0)
print("\n  strong-hit reads by top matching organism:")
for k, v in org.most_common(8):
    print(f"    {k:35} {v:>9,} reads")
