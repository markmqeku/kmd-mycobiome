#!/usr/bin/env python
"""G-3c: report where the yeast ASVs sit relative to Tremellomycetes references.
Placement is SHOWN from the tree and from identities; no genus name is assigned."""
import csv, re
from collections import Counter, defaultdict
D = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/yeast_placement"

# ---- identities ----
best = {}
for r in csv.reader(open(f"{D}/yeast_vs_tremello.tsv"), delimiter="\t"):
    if len(r) < 6: continue
    q, s, pid, qcov = r[0], r[1], float(r[2]), float(r[4])
    if q not in best or pid > best[q][1]:
        best[q] = (s, pid, qcov)
print("=== G-3c: nearest Tremellomycetes reference per yeast ASV ===")
gen = Counter(); reads_by_gen = Counter()
def readsof(q):
    m = re.search(r"reads(\d+)", q); return int(m.group(1)) if m else 0
for q, (s, pid, qcov) in best.items():
    g = s.split("_")[1] if s.startswith("REF_") else "?"
    gen[g] += 1; reads_by_gen[g] += readsof(q)
print(f"ASVs with a Tremellomycetes hit: {len(best)} of 78")
print(f"{'nearest reference genus':24} {'ASVs':>5} {'reads':>10}")
for g, n in gen.most_common():
    print(f"  {g:22} {n:>5} {reads_by_gen[g]:>10,}")
pid_all = [v[1] for v in best.values()]
if pid_all:
    print(f"\nbest-hit identity to any Tremellomycetes reference: "
          f"min={min(pid_all):.1f}  median={sorted(pid_all)[len(pid_all)//2]:.1f}  max={max(pid_all):.1f}")
print("\ntop 10 ASVs by reads:")
for q in sorted(best, key=readsof, reverse=True)[:10]:
    s, pid, qcov = best[q]
    print(f"  {q:34} -> {s[:34]:34} id={pid:5.1f} qcov={qcov:3.0f}")

# ---- tree topology: are the ASVs one clade or several? ----
try:
    tree = open(f"{D}/yeast_placement.tree").read()
except FileNotFoundError:
    print("\n(no tree file)"); raise SystemExit
tips = re.findall(r"[(,]([^(),:]+):", tree)
asv_tips = [t for t in tips if t.startswith("ASV_")]
print(f"\n=== tree: {len(tips)} tips ({len(asv_tips)} ASVs, {len(tips)-len(asv_tips)} references) ===")

# crude clade check: walk the newick, find maximal groups of consecutive ASV tips
order = [("ASV" if t.startswith("ASV_") else "REF") for t in tips]
runs = []
cur = None
for o in order:
    if o == "ASV":
        if cur is None: cur = 1
        else: cur += 1
    else:
        if cur: runs.append(cur); cur = None
if cur: runs.append(cur)
print(f"contiguous ASV blocks in tip order: {len(runs)} (sizes: {sorted(runs, reverse=True)[:10]})")
print("  a single large block means the ASVs form one clade; many small blocks means they are dispersed")
# which references sit adjacent to ASV blocks
adj = Counter()
for i, t in enumerate(tips):
    if t.startswith("ASV_"):
        for j in (i-1, i+1):
            if 0 <= j < len(tips) and tips[j].startswith("REF_"):
                adj[tips[j].split("_")[1]] += 1
print(f"reference genera adjacent to ASV tips: {dict(adj.most_common(8))}")
