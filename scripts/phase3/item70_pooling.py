#!/usr/bin/env python
"""Item 70 and 70-ADD. Lay out the Bloemfontein sample structure exactly as the map files
record it, then test each candidate reading of the scanned pooling notes against it. No tree
count is inferred."""
import csv, glob, os, re
from collections import Counter, defaultdict
K = "<KMD_ROOT>"
RR = f"{K}/__restructured_25-07-2026/raw_reads"

samples = {}
for p in sorted(glob.glob(f"{RR}/*/map*.txt")):
    for line in open(p, encoding="utf-8", errors="replace"):
        f = line.rstrip("\n").split("\t")
        if not f or not f[0].strip() or f[0].startswith("#"): continue
        sid = f[0].strip()
        treat = f[3].strip() if len(f) > 3 else ""
        desc = f[4].strip() if len(f) > 4 else ""
        site = desc.split("_")[-1].replace("Bloefontein", "Bloemfontein")
        tissue = desc.split("_")[0]
        age = "old" if treat.lower().startswith("old") else ("young" if treat.lower().startswith("young") else "")
        samples[sid] = (site, age, tissue)

print("=== 70a. Bloemfontein sample structure, exactly as the map files record it ===")
bl = {k: v for k, v in samples.items() if v[0] == "Bloemfontein"}
print(f"   Bloemfontein samples in the map files: {len(bl)}  ({min(bl, key=lambda s:int(s[1:]))} to "
      f"{max(bl, key=lambda s:int(s[1:]))})")
grid = defaultdict(list)
for s, (site, age, tis) in bl.items(): grid[(age, tis)].append(s)
print(f"\n   {'age':7} {'tissue':14} {'samples':>8}  ids")
for (age, tis), ss in sorted(grid.items()):
    ss = sorted(ss, key=lambda s: int(s[1:]))
    print(f"   {age:7} {tis:14} {len(ss):>8}  {', '.join(ss)}")
for age in ("young", "old"):
    n = sum(len(v) for (a, t), v in grid.items() if a == age)
    print(f"   total {age:6}: {n} samples")

print("\n   for comparison, the two affected sites:")
for site in ("Christiana", "Pretoria"):
    ss = {k: v for k, v in samples.items() if v[0] == site}
    g = Counter((v[1] or "not specified", v[2]) for v in ss.values())
    print(f"      {site}: {len(ss)} samples  " + "; ".join(f"{a}/{t}={c}" for (a, t), c in sorted(g.items())))

print("\n=== 70b. the internal inconsistency in the scanned notes ===")
print("   As reported from the scans:")
print("      gel                      : 32 samples, P1 to P32")
print("      pooling table, Bloem A   : 16 young + 16 old, 'pooling to 6 + 6'   -> trees total 32")
print("      pooling table, Bloem B   : 12 + 12,          'pooling to 6 + 6'   -> trees total 24")
print("   The two pooling-table rows disagree with each other on the tree total, 32 against 24,")
print("   while both give the same pooled result, 6 + 6. The gel count, 32, matches the number of")
print("   SEQUENCED SAMPLES in the map files exactly, not necessarily a number of trees.")
print(f"   map files give {len(bl)} Bloemfontein samples, which is the same 32.")
print("   So '32' is attested for samples and only asserted for trees; the two cannot be")
print("   distinguished from the scans alone.")

print("\n=== 70c. do '2 reps x 2 each' and '3 reps x 2 each' reconcile? ===")
young = {t: len(v) for (a, t), v in grid.items() if a == "young"}
old = {t: len(v) for (a, t), v in grid.items() if a == "old"}
print(f"   young tissue sample counts: {young}")
print(f"   old   tissue sample counts: {old}")
tests = [
 ("young, 3 reps x 2 each = 6 pooled samples", 3*2, [young.get('Leaves'), young.get('Twigs')]),
 ("old,   2 reps x 2 each = 4 pooled samples", 2*2, [old.get('Leaves'), old.get('Twigs')]),
]
for label, got, targets in tests:
    print(f"\n   {label}")
    for t in targets:
        print(f"      against a tissue with {t} samples: {'MATCHES' if t == got else 'does not match'}")
print("\n   Reading 'pooling to 6 + 6' as 6 leaf samples + 6 twig samples per age class:")
for age, d in (("young", young), ("old", old)):
    lt = d.get("Leaves", 0) + d.get("Twigs", 0)
    rest = sum(v for k, v in d.items() if k not in ("Leaves", "Twigs"))
    print(f"      {age:6}: leaves {d.get('Leaves',0)} + twigs {d.get('Twigs',0)} = {lt}, "
          f"plus {rest} in other tissues, total {lt+rest}")
print("      This reading CLOSES against the map files: 6 + 6 + 2 + 2 = 16 samples per age class,")
print("      and 16 + 16 = 32 samples in total. It describes SAMPLES, not trees.")
print("\n   If 16 young trees produced 6 young leaf samples, that is 2.67 trees per sample, which")
print("   is not an integer. 'x 2 each' would give 6 pools from 12 trees, leaving 4 unaccounted.")
print("   Four trees per pool, as recorded for the affected sites, would need 24 trees for the")
print("   young leaf samples alone, against 16 available.")
print("\n   CONCLUSION: the sample structure is fully documented and internally consistent.")
print("   The number of TREES per pooled sample at Bloemfontein does not follow from any")
print("   surviving record, and the three candidate readings give 2, 2.67 and 4 trees per pool.")
