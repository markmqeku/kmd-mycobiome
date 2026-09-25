#!/usr/bin/env python
"""T-3 DRY RUN: rarefaction retention at depth 5,000. Reports only. Computes no diversity statistic.
Grouping variables resolved by header text; fails loudly naming the header if absent."""
import csv, sys
from collections import Counter

DEPTH = 5000
META = sys.argv[1]
STATS = sys.argv[2:]           # DADA2 stats tsvs (non-chimeric == table frequency)

REQUIRED = ["location", "tissue", "condition_std", "age", "run_chemistry"]

freq = {}
for fn in STATS:
    with open(fn) as fh:
        rd = csv.reader(fh, delimiter="\t"); hdr = next(rd); next(rd)
        if "non-chimeric" not in hdr:
            sys.exit(f"FAIL: header 'non-chimeric' not found in {fn}. Headers present: {hdr}")
        i = hdr.index("non-chimeric")
        for row in rd:
            if row and row[0].strip():
                freq[row[0]] = int(float(row[i]))

with open(META) as fh:
    rd = csv.reader(fh, delimiter="\t"); hdr = next(rd); next(rd)
    missing = [c for c in REQUIRED if c not in hdr]
    if missing:
        sys.exit(f"FAIL: metadata header(s) not found: {missing}. Headers present: {hdr}")
    di = {c: i for i, c in enumerate(hdr)}
    meta = {r[0]: {c: r[di[c]] for c in hdr} for r in rd if r and r[0].strip()}

S = [s for s in freq if s in meta]
dropped = sorted([s for s in S if freq[s] < DEPTH], key=lambda s: freq[s])
kept = [s for s in S if freq[s] >= DEPTH]

print(f"=== RAREFACTION DRY RUN @ depth {DEPTH} ===")
print(f"samples with counts: {len(S)}   retained: {len(kept)}   DROPPED: {len(dropped)}")
print("dropped sample(s):", [(s, freq[s]) for s in dropped] or "none")
print(f"lowest retained: {min(((freq[s], s) for s in kept))}")

def block(title, keyfn, subset=None):
    print(f"\n--- {title} (retained n; * = n<3) ---")
    pool = subset if subset is not None else S
    before = Counter(keyfn(s) for s in pool)
    after = Counter(keyfn(s) for s in pool if freq[s] >= DEPTH)
    for k in sorted(before, key=lambda x: tuple(map(str, x if isinstance(x, tuple) else (x,)))):
        n = after.get(k, 0)
        flag = " *" if n < 3 else ""
        chg = "" if n == before[k] else f"   (was {before[k]})"
        print(f"  {str(k):48} {n}{flag}{chg}")

# Article 1: Bloemfontein, tissue x age
bloem = [s for s in S if meta[s]["location"] == "Bloemfontein"]
block("ARTICLE 1  Bloemfontein: tissue x age", lambda s: (meta[s]["tissue"], meta[s]["age"]), bloem)
print("   inferential core = Leaves + Twigs; Inflorescence + Seeds are descriptive-only (locked)")

# Article 2: within-site, matched tissue x condition_std
for site in ["Christiana", "Pretoria"]:
    sub = [s for s in S if meta[s]["location"] == site]
    block(f"ARTICLE 2  {site}: condition_std x tissue", lambda s: (meta[s]["condition_std"], meta[s]["tissue"]), sub)

block("run_chemistry (batch) x location", lambda s: (meta[s]["run_chemistry"], meta[s]["location"]))
print("\nNOTE: no diversity statistic computed. Awaiting 'go'.")
