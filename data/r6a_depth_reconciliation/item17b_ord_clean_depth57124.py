#!/usr/bin/env python
"""H-0 pre-trained classifier availability; H-4 ordination (no PERMANOVA) and effort-corrected sharing."""
import csv, urllib.request, itertools, random
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"

print("=== H-0: pre-trained UNITE classifier availability ===")
for url in ["https://data.qiime2.org/classifiers/sklearn-1.4.2/unite/unite_ver10_dynamic_all_04.04.2024-Q2-2024.5.qza",
            "https://data.qiime2.org/2024.10/common/unite-ver10-dynamic-classifier.qza"]:
    try:
        req = urllib.request.Request(url, method="HEAD", headers={"User-Agent":"KMD/1.0"})
        with urllib.request.urlopen(req, timeout=30) as r:
            print(f"   AVAILABLE ({r.status}): {url}")
    except Exception as e:
        print(f"   not available: {url.split('/')[-1]}  [{type(e).__name__}]")
print("   QIIME 2 ships pre-trained classifiers for SILVA/Greengenes; UNITE classifiers are")
print("   conventionally user-trained, so the rank-threshold fallback stands as a documented limitation.")

# ---- data ----
lines = [l.rstrip("\n") for l in open(f"{O}/clean/feature-table-clean.tsv")
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
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]

DEPTH = 57124
rng = random.Random(42)
def rarefy(counts, depth):
    pool = []
    for a, v in counts.items(): pool.extend([a]*v)
    if len(pool) < depth: return Counter(pool)
    return Counter(rng.sample(pool, depth))
rar = {s: rarefy(per[s], DEPTH) for s in core}

def bray(a, b):
    ks = set(a) | set(b)
    num = sum(abs(a.get(k,0)-b.get(k,0)) for k in ks); den = sum(a.values())+sum(b.values())
    return num/den if den else 0.0
def jac(a, b):
    A, Bs = set(a), set(b)
    return 1 - len(A & Bs)/len(A | Bs) if (A|Bs) else 0.0

print(f"\n=== H-4: ORDINATION on the rarefied core 16 (depth {DEPTH:,}) ===")
print("PCoA coordinates only. NO PERMANOVA, NO p-values (every cell is n=1-2).")
import math
for name, fn in (("Bray-Curtis", bray), ("Jaccard", jac)):
    n = len(core)
    D = [[fn(rar[core[i]], rar[core[j]]) for j in range(n)] for i in range(n)]
    A = [[-0.5*D[i][j]**2 for j in range(n)] for i in range(n)]
    rm = [sum(r)/n for r in A]; cm = [sum(A[i][j] for i in range(n))/n for j in range(n)]
    gm = sum(rm)/n
    G = [[A[i][j]-rm[i]-cm[j]+gm for j in range(n)] for i in range(n)]
    # power iteration for first 2 axes
    def axis(G):
        v = [rng.random() for _ in range(n)]
        for _ in range(500):
            w = [sum(G[i][j]*v[j] for j in range(n)) for i in range(n)]
            nn = math.sqrt(sum(x*x for x in w)) or 1
            v = [x/nn for x in w]
        lam = sum(v[i]*sum(G[i][j]*v[j] for j in range(n)) for i in range(n))
        return lam, v
    l1, v1 = axis(G)
    G2 = [[G[i][j]-l1*v1[i]*v1[j] for j in range(n)] for i in range(n)]
    l2, v2 = axis(G2)
    tot = sum(G[i][i] for i in range(n))
    print(f"\n-- {name} PCoA (axis1 {100*l1/tot:.1f}%, axis2 {100*l2/tot:.1f}% of variation) --")
    pts = sorted(range(n), key=lambda i: v1[i]*math.sqrt(abs(l1)))
    for i in pts:
        s = core[i]
        print(f"   {s:6} {meta[s]['location']:11} {meta[s]['condition_std']:13} {meta[s]['tissue']:14} "
              f"ax1={v1[i]*math.sqrt(abs(l1)):+7.3f} ax2={v2[i]*math.sqrt(abs(l2)):+7.3f}")

# ---- effort-corrected ASV sharing ----
print("\n=== H-4: EFFORT-CORRECTED ASV SHARING (replaces the naive 81% vs 46-48%) ===")
print("Method: within each site, the larger condition group is subsampled WITHOUT replacement to the")
print("smaller group's sample count, 1,000 iterations; ASVs are counted as present if present in any")
print("sample of the subsample. Rarefied tables are used so per-sample read depth is also equal.")
for site in ("Christiana","Pretoria"):
    ss = [s for s in core if meta[s]["location"] == site]
    grp = defaultdict(list)
    for s in ss: grp[meta[s]["condition_std"]].append(s)
    a, b = "asymptomatic", "symptomatic"
    na, nb = len(grp[a]), len(grp[b])
    k = min(na, nb)
    tot_asv = len({x for s in ss for x in rar[s]})
    def naive(g): return len({x for s in grp[g] for x in rar[s]})
    print(f"\n--- {site}: {na} asymptomatic vs {nb} symptomatic; total ASVs detected = {tot_asv} ---")
    print(f"   NAIVE (unequal effort): asymptomatic {naive(a)} ({100*naive(a)/tot_asv:.0f}%), "
          f"symptomatic {naive(b)} ({100*naive(b)/tot_asv:.0f}%)")
    res = {}
    for g in (a, b):
        vals = []
        for _ in range(1000):
            pick = rng.sample(grp[g], k)
            vals.append(len({x for s in pick for x in rar[s]}))
        vals.sort()
        res[g] = (sum(vals)/len(vals), vals[24], vals[-25])
    print(f"   EFFORT-CORRECTED at n={k} per group, equal read depth:")
    for g in (a, b):
        m, lo, hi = res[g]
        print(f"      {g:13} mean richness = {m:6.1f} ASVs   95% range [{lo}, {hi}]")
    print(f"   difference (symptomatic - asymptomatic) = {res[b][0]-res[a][0]:+.1f} ASVs")
