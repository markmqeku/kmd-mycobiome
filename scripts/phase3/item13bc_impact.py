#!/usr/bin/env python
"""13b remove confirmed non-fungal ASVs and report before/after.
13c report the impact on every reported quantity. Nothing in FINAL_NUMBERS is written
by this script: it reports old and new side by side and leaves the decision open."""
import csv, os, zipfile, subprocess, re
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; S = f"{O}/full_screen"
C = f"{O}/clean"; os.makedirs(C, exist_ok=True)
DEPTH = 57124   # DEC-1: re-locked from 57,141 after non-fungal reads were removed from 27-C

drop = [l.strip() for l in open(f"{S}/nonfungal_asvs.txt") if l.strip()]
dropset = set(drop)
print(f"=== 13b REMOVAL: {len(drop)} confirmed non-fungal ASVs ===")

def load_table(p):
    lines = [l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hd = lines[0].lstrip("#").split("\t"); cols = [x for x in hd[1:] if x.strip()]
    per = {s: {} for s in cols}
    for l in lines[1:]:
        q = l.split("\t")
        for s, v in zip(cols, q[1:1+len(cols)]):
            v = int(float(v))
            if v: per[s][q[0]] = v
    return cols, per
cols, per = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
allasv = {a for s in cols for a in per[s]}
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]

before_n, before_reads = len(allasv), sum(sum(per[s].values()) for s in cols)
clean = {s: {a: v for a, v in per[s].items() if a not in dropset} for s in cols}
after_n = len({a for s in cols for a in clean[s]})
after_reads = sum(sum(clean[s].values()) for s in cols)
print(f"ASVs   before {before_n:>6,}   after {after_n:>6,}   removed {before_n-after_n}")
print(f"reads  before {before_reads:>9,}   after {after_reads:>9,}   removed {before_reads-after_reads:,}"
      f"  ({100*(before_reads-after_reads)/before_reads:.3f}% of the fungal table)")

print(f"\nper-sample frequencies (only samples that change):")
print(f"{'sample':7} {'site':12} {'before':>9} {'after':>9} {'lost':>7} {'%':>7} {'>=57,141 before':>16} {'after':>7}")
changed = []
for s in cols:
    b, a = sum(per[s].values()), sum(clean[s].values())
    if b != a:
        changed.append(s)
        print(f"{s:7} {meta.get(s,{}).get('location','?'):12} {b:>9,} {a:>9,} {b-a:>7,} "
              f"{100*(b-a)/b:>6.2f}% {'yes' if b>=DEPTH else 'no':>16} {'yes' if a>=DEPTH else 'no':>7}")
if not changed: print("   none")
kb = [s for s in core if sum(per[s].values()) >= DEPTH]
ka = [s for s in core if sum(clean[s].values()) >= DEPTH]
print(f"\ncore-16 retention at depth {DEPTH:,}: before {len(kb)}/16, after {len(ka)}/16 "
      f"-> {'UNCHANGED' if kb == ka else 'CHANGED: ' + str(set(kb) ^ set(ka))}")
core_hit = [a for a in drop if any(a in per[s] for s in core)]
print(f"removed ASVs occurring in any core-16 sample: {len(core_hit)}"
      + ("" if core_hit else "  -> Hill numbers, ordination and Curvibasidium cannot change"))
for a in core_hit:
    ss = [s for s in core if a in per[s]]
    print(f"   {a[:12]} in {len(ss)} core samples: {ss}, {sum(per[s][a] for s in ss):,} reads")

with open(f"{C}/feature-table-clean.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["#OTU ID"] + cols)
    for a in sorted({x for s in cols for x in clean[s]}):
        w.writerow([a] + [clean[s].get(a, 0) for s in cols])
print(f"\ncorrected table -> {C}/feature-table-clean.tsv")

# ---------------- 13c impact ----------------
print(f"\n\n=== 13c IMPACT ON REPORTED QUANTITIES ===")
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
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
YEAST = "Tremellomycetes yeast lineage (ITS2-unresolved)"
def label(a):
    if a in yeast: return YEAST
    rk = ranks_of(final.get(a, "")); last = ("none", "")
    for _, r in RANKS:
        if rk[r]: last = (r, rk[r])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"

groups = [("Christiana","asymptomatic"),("Christiana","symptomatic"),
          ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
print("\ncomposition, top labels by group, before vs after:")
for g in groups:
    ss = [s for s in core if meta[s]["location"] == g[0] and meta[s]["condition_std"] == g[1]]
    for tag, tab in (("before", per), ("after", clean)):
        tot = sum(sum(tab[s].values()) for s in ss)
        agg = Counter()
        for s in ss:
            for a, v in tab[s].items(): agg[label(a)] += v
        if tag == "before": bfr = (tot, agg)
        else:
            same = all(abs(100*bfr[1][k]/bfr[0] - 100*agg[k]/tot) < 0.005 for k in set(bfr[1]) | set(agg))
            print(f"   {g[0]:11} {g[1]:13} n={len(ss)}  reads {bfr[0]:>8,} -> {tot:>8,}  "
                  f"{'IDENTICAL to 2 dp' if same else 'CHANGED'}")
            if not same:
                for k in sorted(set(bfr[1]) | set(agg), key=lambda k: -agg[k])[:8]:
                    print(f"      {k[:40]:40} {100*bfr[1][k]/bfr[0]:6.2f}% -> {100*agg[k]/tot:6.2f}%")

print("\nCurvibasidium relative abundance per core sample, before vs after:")
# Use I-2's definition, not the final taxonomy: I-2 defined the group from the raw UNITE 10.0
# assignment at >=95% identity. Checking the final taxonomy instead silently returns zero ASVs
# and would make this test vacuous.
tax10 = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax10[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
_fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(_fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
cur = [a for a in allasv if ranks_of(tax10.get(a, ""))["genus"] == "Curvibasidium" and pid.get(a, 0) >= 95]
print(f"   Curvibasidium ASVs (UNITE 10.0 genus, >=95% identity): {len(cur)}; "
      f"any removed? {'YES' if set(cur) & dropset else 'no'}")
if not cur: raise SystemExit("ABORT: Curvibasidium set is empty, the test would be vacuous.")
diff = []
for s in core:
    b = 100*sum(per[s].get(a,0) for a in cur)/max(1,sum(per[s].values()))
    a_ = 100*sum(clean[s].get(a,0) for a in cur)/max(1,sum(clean[s].values()))
    if abs(b-a_) >= 0.005: diff.append((s, b, a_))
print(f"   samples whose value changes by >=0.01 pp: {len(diff)}")
for s, b, a_ in diff: print(f"      {s:7} {b:6.2f}% -> {a_:6.2f}%")

print("\nBloemfontein low-host inventory, before vs after:")
cols_a, per_a = load_table(f"{O}/table_exp/feature-table.tsv")
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
sel = []
for s in [x for x in cols_a if meta.get(x) and meta[x]["location"] == "Bloemfontein"]:
    tot = sum(per_a[s].values()); hs = sum(v for a, v in per_a[s].items() if a in host)
    if tot and 100*hs/tot <= 10: sel.append(s)
def inv(tab):
    """occ counts (sample, ASV) occurrences; dis counts DISTINCT ASVs. H-5 and FINAL_NUMBERS
    section 10 report the occurrence count but label it 'ASVs', which overstates every entry."""
    occ, dis = Counter(), defaultdict(set)
    for s in sel:
        for a in tab.get(s, {}):
            occ[label(a)] += 1; dis[label(a)].add(a)
    return occ, {k: len(v) for k, v in dis.items()}
(ob, db), (oa, da) = inv(per), inv(clean)
print(f"   samples: {sorted(sel)}")
print(f"   distinct ASVs present: {sum(db.values())} -> {sum(da.values())}")
print(f"   (sample, ASV) occurrences: {sum(ob.values())} -> {sum(oa.values())}")
print(f"   distinct taxa: {len(ob)} -> {len(oa)}   ({'unchanged' if len(ob)==len(oa) else 'CHANGED'})")
lost = set(ob) - set(oa)
if lost: print(f"   labels lost entirely: {sorted(lost)}")
print("\n   most frequent labels on the corrected table (distinct ASVs, then occurrences):")
for k, v in sorted(da.items(), key=lambda kv: -kv[1])[:10]:
    print(f"      {k[:50]:50} {v:>3} ASVs   ({oa[k]} occurrences)")
