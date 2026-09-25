#!/usr/bin/env python
"""Item 23b and 23c. Factual description of where 25-C, 26-C and 27-C sit, and what they
contain. NO LABEL IS CHANGED by this script and none should be: these three samples are
inputs to the ordination being described, so using that ordination to reclassify them would
be circular. Distances are reported to the centroid of each group with the sample itself
EXCLUDED from its own group's centroid, otherwise a sample is pulled toward its own label."""
import csv, math, os, random
from collections import Counter
import zipfile
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"
DEPTH = 57124
TARGET = ["25-C", "26-C", "27-C"]

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
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
chr_a = [s for s in core if meta[s]["location"] == "Christiana" and meta[s]["condition_std"] == "asymptomatic"]
chr_s = [s for s in core if meta[s]["location"] == "Christiana" and meta[s]["condition_std"] == "symptomatic"]

rng = random.Random(42)
def rarefy(c, d):
    pool = []
    for a, v in c.items(): pool.extend([a]*v)
    return Counter(pool if len(pool) < d else rng.sample(pool, d))
rar = {s: rarefy(per[s], DEPTH) for s in core}
def bray(a, b):
    ks = set(a) | set(b); num = sum(abs(a.get(k,0)-b.get(k,0)) for k in ks)
    den = sum(a.values())+sum(b.values()); return num/den if den else 0.
def jac(a, b):
    A, Bs = set(a), set(b); return 1-len(A & Bs)/len(A | Bs) if (A | Bs) else 0.

print("=== 23d. group composition ===")
print(f"Christiana asymptomatic: {chr_a}  (n={len(chr_a)})")
print(f"Christiana symptomatic : {chr_s}  (n={len(chr_s)})")
print(f"The three samples in question ARE the entire Christiana asymptomatic group.")
print(f"Reclassifying them would leave Christiana with 0 asymptomatic samples and would remove")
print(f"the site's only within-site contrast, collapsing the core set from 16 to 13 usable for")
print(f"any Christiana condition comparison.\n")

print("=== 23b. mean dissimilarity to each Christiana group (self excluded) ===")
rows = []
for nm, fnc in (("Bray-Curtis", bray), ("Jaccard", jac)):
    print(f"\n-- {nm} --")
    print(f"{'sample':7} {'to asym group':>15} {'to sym group':>14} {'difference':>12}  closer to")
    for s in TARGET:
        oa = [x for x in chr_a if x != s]
        da = sum(fnc(rar[s], rar[x]) for x in oa)/len(oa)
        ds = sum(fnc(rar[s], rar[x]) for x in chr_s)/len(chr_s)
        print(f"{s:7} {da:>15.4f} {ds:>14.4f} {ds-da:>+12.4f}  {'asymptomatic' if da < ds else 'SYMPTOMATIC'}")
        rows.append([nm, s, f"{da:.4f}", f"{ds:.4f}", f"{ds-da:+.4f}",
                     "asymptomatic" if da < ds else "symptomatic"])
    # reference: how the symptomatic samples themselves score
    print(f"   for reference, the 6 symptomatic samples' own mean distance to the other")
    print(f"   symptomatic samples: ", end="")
    v = [sum(fnc(rar[s], rar[x]) for x in chr_s if x != s)/(len(chr_s)-1) for s in chr_s]
    print(f"{sum(v)/len(v):.4f} (range {min(v):.4f} to {max(v):.4f})")

print("\n=== 23c. Curvibasidium and yeast lineage ===")
tax10 = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax10[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
def genus(a):
    for part in tax10.get(a, "").split(";"):
        if part.strip().startswith("g__"): return part.strip()[3:]
    return ""
curv = {a for a in tax10 if genus(a) == "Curvibasidium" and pid.get(a, 0) >= 95}
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
def pct(s, S):
    t = sum(per[s].values()); return 100*sum(v for a, v in per[s].items() if a in S)/t if t else 0
print(f"{'sample':7} {'group':14} {'Curvibasidium %':>16} {'yeast lineage %':>16}")
for s in chr_a + chr_s:
    g = meta[s]["condition_std"] + (" *" if s in TARGET else "")
    print(f"{s:7} {g:14} {pct(s, curv):>16.2f} {pct(s, yeast):>16.2f}")
for nm, grp in (("asymptomatic (the 3)", chr_a), ("symptomatic (n=6)", chr_s)):
    c = [pct(s, curv) for s in grp]; y = [pct(s, yeast) for s in grp]
    print(f"   {nm:22} Curvibasidium mean {sum(c)/len(c):6.2f} (range {min(c):.2f}-{max(c):.2f})"
          f" | yeast mean {sum(y)/len(y):6.2f} (range {min(y):.2f}-{max(y):.2f})")

with open(f"{C}/item23_position.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["metric","sample","mean_dist_to_asymptomatic","mean_dist_to_symptomatic",
                "difference","closer_to"])
    w.writerows(rows)
print(f"\n-> {C}/item23_position.tsv")
print("NO LABELS CHANGED.")
