#!/usr/bin/env python
"""G-4b fallback (sklearn OOM-infeasible): rank-specific identity cutoffs applied to the permissive
vsearch run, truncating each lineage to the deepest rank its best-hit identity supports.
Cutoffs (conservative ends of the ranges specified): species >=97, genus >=95, family >=90,
order >=85, class/phylum/kingdom >=80 (the vsearch floor).
G-4b cross-check: unassigned ASVs above 0.5% relative abundance vs NCBI nt.
G-4c comparison table.  G-5 composition on the fungi-only table."""
import csv, zipfile, os
from collections import defaultdict, Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUTOFF = {"species":97.0,"genus":95.0,"family":90.0,"order":85.0,"class":80.0,"phylum":80.0,"kingdom":80.0}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

def ranks_of(t):
    out = {r:"" for _,r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part = part.strip()
        for pre,rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk=="species" and (b.endswith("_sp") or b.endswith("_sp.")))): v = ""
                out[rk] = v
    return out

def load_tax(p):
    d = {}
    for r in csv.reader(open(p), delimiter="\t"):
        if r and r[0] not in ("Feature ID",) and not r[0].startswith("#"): d[r[0]] = r[1]
    return d

def best_pident(qza):
    z = zipfile.ZipFile(qza)
    f = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
    b = {}
    for line in z.read(f).decode("utf-8","replace").splitlines():
        p = line.split("\t")
        if len(p) < 3: continue
        try: pid = float(p[2])
        except ValueError: continue
        if p[0] not in b or pid > b[p[0]]: b[p[0]] = pid
    return b

tax80 = load_tax(f"{O}/taxonomy_unite10/exp/taxonomy.tsv")
tax97 = load_tax(f"{O}/taxonomy_unite10_id97/exp/taxonomy.tsv")
pid   = best_pident(f"{O}/taxonomy_unite10/search.qza")

# fungi-only ASV set
fungi = set()
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    if line.startswith(">"): fungi.add(line[1:].split()[0])

# ---- build the rank-truncated hybrid taxonomy ----
hybrid = {}
for a in fungi:
    rk = ranks_of(tax80.get(a,""))
    p = pid.get(a, 0.0)
    keep = {}
    for _, r in RANKS:
        if rk[r] and p >= CUTOFF[r]: keep[r] = rk[r]
        else: break                      # truncate at first unsupported rank
    hybrid[a] = keep

def deepest(d):
    last = (None,None)
    for _,r in RANKS:
        if d.get(r): last = (r, d[r])
    return last

print("=== G-4c: RANK RESOLUTION COMPARISON (fungi-only table, n=%d) ===" % len(fungi))
print(f"{'rank':9} {'vsearch 0.80':>14} {'vsearch 0.97':>14} {'hybrid (cutoffs)':>18}")
for _, r in RANKS:
    a80 = sum(1 for a in fungi if ranks_of(tax80.get(a,""))[r])
    a97 = sum(1 for a in fungi if ranks_of(tax97.get(a,""))[r])
    ah  = sum(1 for a in fungi if hybrid[a].get(r))
    print(f"{r:9} {a80:>6} ({100*a80/len(fungi):4.1f}%) {a97:>6} ({100*a97/len(fungi):4.1f}%) "
          f"{ah:>8} ({100*ah/len(fungi):4.1f}%)")
print("sklearn (UNITE 10.0, conf 0.7): NOT AVAILABLE - training OOM-killed (7.3 GB resident, 29 GB virtual)")

# ---- abundances on the fungi-only table ----
lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
per = {s: {} for s in cols}; tot = Counter()
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v; tot[q[0]] += v
grand = sum(tot.values())

# ---- G-4b: unassigned ASVs above 0.5% relative abundance, cross-checked vs NCBI ----
print("\n=== G-4b: high-abundance ASVs unassigned by UNITE, cross-checked against NCBI nt ===")
ncbi = {}
pth = f"{O}/nontarget_screen/blast_priority.tsv"
if os.path.exists(pth):
    for r in csv.reader(open(pth), delimiter="\t"):
        if len(r) < 8: continue
        if r[0] not in ncbi or float(r[2]) > ncbi[r[0]][0]: ncbi[r[0]] = (float(r[2]), float(r[4]), r[7])
big = [a for a in fungi if tot[a]/grand >= 0.005]
print(f"ASVs >=0.5% of fungal reads: {len(big)}")
print(f"{'ASV':10} {'%reads':>7} {'UNITE(hybrid) deepest':32} {'NCBI nt top hit':44}")
gain = 0
for a in sorted(big, key=lambda x: -tot[x]):
    r, v = deepest(hybrid[a])
    u = f"{v} ({r})" if r else "UNASSIGNED"
    n = ncbi.get(a)
    nn = f"{n[2][:40]} id={n[0]:.1f}" if n else "-"
    if not r and n: gain += 1
    print(f"{a[:9]:10} {100*tot[a]/grand:>6.2f}% {u:32} {nn:44}")
print(f"\nASVs unassigned by UNITE that gain a confident NCBI hit: {gain}")

# ---- G-5: composition per site x condition ----
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

print("\n\n=== G-5: TAXONOMIC COMPOSITION (fungi-only table, hybrid taxonomy) ===")
print("relative abundance of the deepest supported name; ALL earlier composition figures are obsolete\n")
groups = defaultdict(lambda: Counter())
gtot = Counter()
for s in cols:
    m = meta.get(s)
    if not m: continue
    key = (m["location"], m["condition_std"])
    for a, v in per[s].items():
        r, name = deepest(hybrid[a])
        lbl = f"{name} ({r})" if r else "unassigned"
        groups[key][lbl] += v; gtot[key] += v
for key in sorted(groups):
    print(f"--- {key[0]} / {key[1]}  (total {gtot[key]:,} reads, n={sum(1 for s in cols if meta.get(s) and (meta[s]['location'],meta[s]['condition_std'])==key)} samples) ---")
    for lbl, v in groups[key].most_common(10):
        print(f"     {lbl:44} {100*v/gtot[key]:5.1f}%")
    print()
