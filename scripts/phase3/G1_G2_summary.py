#!/usr/bin/env python
"""G-1 fungi-only table summary + retention; G-2 full Bloemfontein per-sample table."""
import csv
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}

def load_table(p):
    """Return {sample: {asv: count}} keyed by NAME, never by column position,
    because the fungi-only and full tables have different sample sets and orders."""
    lines = [l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
    per = {s: {} for s in cols}
    for l in lines[1:]:
        q = l.split("\t")
        for s, v in zip(cols, q[1:1+len(cols)]):
            v = int(float(v))
            if v: per[s][q[0]] = v
    return cols, per

samples_f, fung_s = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
samples_a, full_s = load_table(f"{O}/table_exp/feature-table.tsv")
samples = samples_a                      # report over ALL original samples
missing = [s for s in samples_a if s not in samples_f]
print(f"NOTE: samples absent from the fungi-only table (zero fungal reads): {missing}\n")
fung = {}                                # asv -> per-sample list aligned to `samples`
for a in {a for d in fung_s.values() for a in d}:
    fung[a] = [fung_s.get(s, {}).get(a, 0) for s in samples]
full = {}
for a in {a for d in full_s.values() for a in d}:
    full[a] = [full_s.get(s, {}).get(a, 0) for s in samples]

meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

fs = {s: sum(fung[a][j] for a in fung) for j, s in enumerate(samples)}
fa = {s: sum(1 for a in fung if fung[a][j] > 0) for j, s in enumerate(samples)}
hs = {s: sum(full[a][j] for a in full if a in host) for j, s in enumerate(samples)}
ts = {s: sum(full[a][j] for a in full) for j, s in enumerate(samples)}

print("=== G-1: FUNGI-ONLY TABLE ===")
print(f"ASVs                 : {len(fung)}   (1,322 minus 149 host)")
print(f"total fungal reads   : {sum(fs.values()):,}   (was 4,010,617 including host)")
print(f"samples              : {len(samples)}")
srt = sorted(fs.items(), key=lambda x: x[1])
print(f"minimum frequency    : {srt[0][0]} = {srt[0][1]:,}")
print(f"maximum frequency    : {srt[-1][0]} = {srt[-1][1]:,}")
print(f"median frequency     : {sorted(fs.values())[len(fs)//2]:,}")
print("\nfive lowest:", [(s, v) for s, v in srt[:5]])

print("\n--- retention by depth and site ---")
print(f"{'depth':>7} {'TOTAL':>9} {'Bloemfontein':>14} {'Christiana':>12} {'Pretoria':>10}")
for d in (1000, 2500, 5000, 10000):
    row = [d, sum(1 for s in samples if fs[s] >= d)]
    for site in ("Bloemfontein", "Christiana", "Pretoria"):
        sub = [s for s in samples if meta[s]["location"] == site]
        row.append(f"{sum(1 for s in sub if fs[s] >= d)}/{len(sub)}")
    print(f"{row[0]:>7} {row[1]:>6}/45 {row[2]:>14} {row[3]:>12} {row[4]:>10}")

print("\n\n=== G-2: BLOEMFONTEIN, ALL 29 SAMPLES (sorted by fungal ASV count) ===")
bl = [s for s in samples if meta[s]["location"] == "Bloemfontein"]
print(f"{'sample':7} {'tissue':14} {'age':6} {'total':>9} {'host':>9} {'host%':>7} {'fungal':>9} {'fungalASVs':>11}")
for s in sorted(bl, key=lambda x: -fa[x]):
    hp = 100*hs[s]/ts[s] if ts[s] else 0
    print(f"{s:7} {meta[s]['tissue']:14} {meta[s]['age']:6} {ts[s]:>9,} {hs[s]:>9,} {hp:>6.1f}% {fs[s]:>9,} {fa[s]:>11}")

print("\n--- do the fungal-rich samples group by tissue or age? (top 10 by fungal ASV count) ---")
top = sorted(bl, key=lambda x: -fa[x])[:10]
print("  tissue:", dict(Counter(meta[s]["tissue"] for s in top)))
print("  age   :", dict(Counter(meta[s]["age"] for s in top)))
bot = sorted(bl, key=lambda x: fa[x])[:10]
print("  bottom-10 tissue:", dict(Counter(meta[s]["tissue"] for s in bot)))
print("  bottom-10 age   :", dict(Counter(meta[s]["age"] for s in bot)))
print("\n  mean fungal ASVs by tissue:")
for t in sorted({meta[s]["tissue"] for s in bl}):
    v = [fa[s] for s in bl if meta[s]["tissue"] == t]
    print(f"    {t:14} n={len(v):2}  mean={sum(v)/len(v):6.1f}  range={min(v)}-{max(v)}")
print("  mean fungal ASVs by age:")
for a in sorted({meta[s]["age"] for s in bl}):
    v = [fa[s] for s in bl if meta[s]["age"] == a]
    print(f"    {a:14} n={len(v):2}  mean={sum(v)/len(v):6.1f}  range={min(v)}-{max(v)}")
