#!/usr/bin/env python
"""Item 71 supporting checks. Verify the specific numeric and source claims the reasoning
review turns on, so the review rests on recomputation rather than impression."""
import csv, re, zipfile
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
C = f"{O}/clean"
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
lines = [l.rstrip("\n") for l in open(f"{C}/feature-table-clean.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [x for x in hdr[1:] if x.strip()]
per = {s: {} for s in cols}
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
def gen(a):
    for part in tax.get(a, "").split(";"):
        if part.strip().startswith("g__"): return part.strip()[3:]
    return ""
curv = {a for a in tax if gen(a) == "Curvibasidium" and pid.get(a, 0) >= 95}
def pc(s):
    t = sum(per[s].values()); return 100*sum(v for a, v in per[s].items() if a in curv)/t if t else 0

print("=== Results 3.7 Curvibasidium values, recomputed on the corrected table ===")
def med(v):
    v = sorted(v); n = len(v)
    return v[n//2] if n % 2 else (v[n//2-1]+v[n//2])/2
for site in ("Christiana", "Pretoria"):
    for cond in ("asymptomatic", "symptomatic"):
        ss = [s for s in cols if meta.get(s) and meta[s]["location"] == site
              and meta[s]["condition_std"] == cond]
        v = [pc(s) for s in ss]
        print(f"   {site:11} {cond:13} n={len(v)}  values {[f'{x:.2f}' for x in sorted(v, reverse=True)]}")
        print(f"      mean {sum(v)/len(v):.2f}   median {med(v):.2f}")
print("\n   The draft states: 'Averaged across each group, it rose from 0.24 to 5.82 percent at")
print("   Christiana and from 4.88 to 11.28 percent at Pretoria.' Compare against the above.")

print("\n=== Results 3.2 claim: 'retained all samples at every rarefaction depth tested' ===")
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
d = {s: sum(per[s].values()) for s in core}
lo = min(d.values()); who = [s for s in core if d[s] == lo][0]
print(f"   smallest core sample: {who} at {lo:,} reads")
for step in (57124, 58000, 60000, 65000):
    print(f"      depth {step:>7,}: {sum(1 for s in core if d[s] >= step)}/16 retained")
print("   Section 3.5 states the depth is capped at 57,124 because any greater depth loses 27-C.")

print("\n=== the yeast lineage: 71 ASVs overall against how many in the core 16 ===")
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
drop = {l.strip() for l in open(f"{O}/full_screen/nonfungal_asvs.txt") if l.strip()}
keep = yeast - drop
allsamp = [s for s in cols]
in_all = {a for s in allsamp for a in per[s] if a in keep}
in_core = {a for s in core for a in per[s] if a in keep}
print(f"   lineage ASVs retained overall: {len(keep)}")
print(f"   present in any of the 45 samples: {len(in_all)}")
print(f"   present in the 16 core samples : {len(in_core)}")
print(f"   therefore Bloemfontein-only    : {len(in_all - in_core)}")

print("\n=== Swanepoel on Jami et al. 2013, exact wording ===")
t = open(f"{C}/swanepoel2018_text.txt", encoding="utf-8").read()
for m in re.finditer(r"Jami", t):
    seg = " ".join(t[max(0, m.start()-330):m.start()+90].split())
    print(f"   ...{seg}...")
    break
