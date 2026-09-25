#!/usr/bin/env python
"""Item 64. The low-abundance fungal fraction. Core 16 samples, fungi-only corrected table,
rarefied to the fixed depth of 57,124. No asymptotic estimators: DADA2 removes singletons by
design, so estimators that depend on singleton and doubleton frequencies are not interpretable
for this dataset.

An ASV is RARE IN A SAMPLE if it holds below 1 percent of that sample's reads. Because an ASV
can be rare in one sample and abundant in another, two views are reported and never conflated:
  per-sample  - the rare set within each sample
  always-rare - ASVs that never reach 1 percent in any core sample
"""
import csv, json, os, random, zipfile
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"
OUT = f"{C}/rare_fraction"; os.makedirs(OUT, exist_ok=True)
DEPTH, THRESH = 57124, 0.01

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
assert len(core) == 16
rng = random.Random(42)
def rarefy(c, d):
    pool = []
    for a, v in c.items(): pool.extend([a]*v)
    return Counter(rng.sample(pool, d))
rar = {s: rarefy(per[s], DEPTH) for s in core}
GROUPS = [("Christiana","asymptomatic"),("Christiana","symptomatic"),
          ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
def grp(s): return (meta[s]["location"], meta[s]["condition_std"])

# ---------- taxonomy ----------
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}
def ranks_of(t):
    out = {r: "" for _, r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part = part.strip()
        for pre, rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk == "species" and (b.endswith("_sp") or b.endswith("_sp.")))): v = ""
                out[rk] = v
    return out
final = {}
for r in csv.reader(open(f"{O}/final_taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): final[r[0]] = r[1]
def deepest(a):
    rk = ranks_of(final.get(a, "")); last = "none"
    for _, r in RANKS:
        if rk[r]: last = r
    return last
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
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

# ---------- 64a ----------
print("=== 64a. size of the rare fraction, threshold 1 percent within a sample ===")
rare_in, abund_in = {}, {}
for s in core:
    tot = sum(rar[s].values())
    rare_in[s] = {a for a, v in rar[s].items() if v/tot < THRESH}
    abund_in[s] = set(rar[s]) - rare_in[s]
print(f"{'sample':7} {'site':11} {'cond':13} {'ASVs':>6} {'rare':>6} {'%ASVs':>7} {'rare reads %':>13}")
rows = []
for s in sorted(core, key=lambda x: (grp(x), x)):
    tot = sum(rar[s].values()); n = len(rar[s]); nr = len(rare_in[s])
    rr = sum(v for a, v in rar[s].items() if a in rare_in[s])
    print(f"{s:7} {meta[s]['location']:11} {meta[s]['condition_std']:13} {n:>6} {nr:>6} "
          f"{100*nr/n:>6.1f}% {100*rr/tot:>12.2f}%")
    rows.append([s, meta[s]["location"], meta[s]["condition_std"], n, nr,
                 f"{100*nr/n:.1f}", f"{100*rr/tot:.2f}"])
allrare = set().union(*rare_in.values())
always = {a for a in allrare if all(rar[s].get(a, 0)/DEPTH < THRESH for s in core)}
everabund = allrare - always
tot_all = sum(sum(rar[s].values()) for s in core)
rr_all = sum(v for s in core for a, v in rar[s].items() if a in rare_in[s])
present = set().union(*[set(rar[s]) for s in core])
print(f"\noverall across the 16 core samples:")
print(f"   ASVs present            {len(present):>6}")
print(f"   rare in at least one    {len(allrare):>6} ({100*len(allrare)/len(present):.1f}%)")
print(f"   ALWAYS rare, never >=1% {len(always):>6} ({100*len(always)/len(present):.1f}%)")
print(f"   rare in one sample but abundant in another {len(everabund)}")
print(f"   reads in the rare fraction {rr_all:,} of {tot_all:,} ({100*rr_all/tot_all:.2f}%)")
print("\nby group:")
for g in GROUPS:
    ss = [s for s in core if grp(s) == g]
    u = set().union(*[rare_in[s] for s in ss])
    rr = sum(v for s in ss for a, v in rar[s].items() if a in rare_in[s])
    tt = sum(sum(rar[s].values()) for s in ss)
    print(f"   {g[0]:11} {g[1]:13} n={len(ss)}  distinct rare ASVs {len(u):>4}  "
          f"reads {100*rr/tt:5.2f}%")

# ---------- 64b ----------
print("\n=== 64b. taxonomic composition of the always-rare fraction ===")
dr = Counter(deepest(a) for a in always)
da = Counter(deepest(a) for a in (present - always))
print(f"{'deepest rank':10} {'always-rare':>12} {'  %':>6}   {'rest':>8} {'  %':>6}")
for r in ["none","kingdom","phylum","class","order","family","genus","species"]:
    print(f"   {r:10} {dr[r]:>9} {100*dr[r]/max(1,len(always)):>6.1f}%   "
          f"{da[r]:>6} {100*da[r]/max(1,len(present)-len(always)):>6.1f}%")
for rank in ("phylum","class","family"):
    c = Counter(ranks_of(final.get(a,""))[rank] or "(unresolved)" for a in always)
    print(f"\n   dominant {rank} among always-rare ASVs:")
    for k, v in c.most_common(6): print(f"      {k[:34]:34} {v:>4} ASVs")
def label(a):
    if a in yeast: return "Tremellomycetes yeast lineage"
    rk = ranks_of(final.get(a,"")); last = ("none","")
    for _, r in RANKS:
        if rk[r]: last = (r, rk[r])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"
print(f"\n   distinct taxa represented among always-rare ASVs: {len(set(label(a) for a in always))}")

# ---------- 64c ----------
print("\n=== 64c. assignment rate, always-rare vs the rest ===")
rest = present - always
for rank in ("phylum","class","order","family","genus"):
    ar = sum(1 for a in always if ranks_of(final.get(a,""))[rank])
    ab = sum(1 for a in rest if ranks_of(final.get(a,""))[rank])
    print(f"   {rank:8} always-rare {100*ar/len(always):5.1f}%   rest {100*ab/len(rest):5.1f}%   "
          f"difference {100*ar/len(always)-100*ab/len(rest):+5.1f} pp")
ur = sum(1 for a in always if deepest(a) in ("none",))
ub = sum(1 for a in rest if deepest(a) in ("none",))
print(f"\n   wholly unassigned: always-rare {ur} of {len(always)} ({100*ur/len(always):.1f}%), "
      f"rest {ub} of {len(rest)} ({100*ub/len(rest):.1f}%)")
print(f"   VERDICT: unassigned sequences are "
      f"{'MORE' if 100*ur/len(always) > 100*ub/len(rest) else 'NOT more'} common among rare ASVs")

# ---------- 64d ----------
print("\n=== 64d. rare fraction by condition and site ===")
print(f"{'group':28} {'n':>3} {'rare ASVs':>10} {'rare reads %':>13} {'distinct taxa':>14}")
for g in GROUPS:
    ss = [s for s in core if grp(s) == g]
    u = set().union(*[rare_in[s] for s in ss])
    rr = sum(v for s in ss for a, v in rar[s].items() if a in rare_in[s])
    tt = sum(sum(rar[s].values()) for s in ss)
    print(f"   {g[0]+' '+g[1]:25} {len(ss):>3} {len(u):>10} {100*rr/tt:>12.2f}% "
          f"{len(set(label(a) for a in u)):>14}")

# ---------- 64e ----------
print("\n=== 64e. prevalence of always-rare ASVs ===")
prev = Counter()
for a in always: prev[a] = sum(1 for s in core if rar[s].get(a, 0) > 0)
dist = Counter(prev.values())
print(f"   in 1 sample only      {dist[1]:>5} ({100*dist[1]/len(always):.1f}%)")
print(f"   in 2 to 4 samples     {sum(dist[i] for i in range(2,5)):>5}")
print(f"   in 5 to 8 samples     {sum(dist[i] for i in range(5,9)):>5}")
print(f"   in 9 to 15 samples    {sum(dist[i] for i in range(9,16)):>5}")
print(f"   in all 16 samples     {dist[16]:>5}")
wide = [a for a in always if prev[a] >= 12]
print(f"\n   always-rare ASVs present in 12 or more of the 16 samples: {len(wide)}")
for a in sorted(wide, key=lambda x: -prev[x])[:12]:
    print(f"      {a[:12]}  in {prev[a]:>2} samples  {label(a)[:48]}")

# ---------- 64f ----------
print("\n=== 64f. the yeast lineage and Curvibasidium in the rare fraction ===")
yr = [a for a in always if a in yeast]
ya = [a for a in present if a in yeast]
print(f"   yeast-lineage ASVs present in the core 16: {len(ya)}")
print(f"   of those, ALWAYS rare: {len(yr)}  ({100*len(yr)/max(1,len(ya)):.1f}%)")
print(f"   so the lineage has both abundant and rare members: {'yes' if yr and len(ya)>len(yr) else 'no'}")
cr = [a for a in present if a in curv]
crr = [a for a in cr if a in always]
insamp = sum(1 for s in core if any(a in rare_in[s] for a in cr))
print(f"   Curvibasidium ASVs in the core 16: {len(cr)}; always rare: {len(crr)}")
print(f"   samples where at least one Curvibasidium ASV falls in the rare fraction: {insamp} of 16")

with open(f"{OUT}/rare_per_sample.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample","site","condition","asvs_total","asvs_rare","pct_asvs_rare","pct_reads_rare"])
    w.writerows(rows)
json.dump({"always_rare": sorted(always), "depth": DEPTH, "threshold": THRESH},
          open(f"{OUT}/always_rare_asvs.json", "w"))
print(f"\n-> {OUT}/rare_per_sample.tsv and always_rare_asvs.json")
