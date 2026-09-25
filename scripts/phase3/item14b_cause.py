#!/usr/bin/env python
"""Item 14 continued. The two label sets differ mostly in formatting ('Fusarium' vs
'Fusarium (genus)'). Strip the rank suffix and compare NAMES, so the 192 -> 199 change
is attributed to real taxonomic differences rather than to cosmetics."""
import csv, zipfile, re
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
exec(open(f"{B}/phase3/item14_reconcile.py").read().split("# ---- B.")[0]
     .replace('print("=== A. SAMPLE SET ===")', 'pass  #'), globals())

RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
YEAST = "Tremellomycetes yeast lineage (ITS2-unresolved)"
def label_frozen(a):
    if a in yeast: return YEAST
    rk = ranks_of(tax10.get(a,""))
    if rk["genus"] and pid.get(a,0) >= 95: return rk["genus"]
    for _, r in reversed(RANKS[:5]):
        if rk[r] and pid.get(a,0) >= 85: return f"{rk[r]} ({r})"
    return "unassigned"
def label_final(a):
    if a in yeast: return YEAST
    rk = ranks_of(final.get(a,"")); last = ("none","")
    for _, r in RANKS:
        if rk[r]: last = (r, rk[r])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"
def norm(s): return re.sub(r"\s*\((kingdom|phylum|class|order|family|genus|species)\)$", "", s)

asvs = sorted({a for s, _ in sel for a in per_f.get(s, {})})
fro = Counter(norm(label_frozen(a)) for a in asvs)
fin = Counter(norm(label_final(a)) for a in asvs)
print(f"ASVs across the 6 low-host Bloemfontein samples: {len(asvs)}")
print(f"distinct taxon NAMES, frozen basis: {len(fro)}")
print(f"distinct taxon NAMES, final basis : {len(fin)}")
print(f"net change: {len(fin)-len(fro):+d}\n")

gone = sorted(set(fro) - set(fin)); new = sorted(set(fin) - set(fro))
print(f"names lost under the final basis ({len(gone)}):")
for k in gone: print(f"   - {k:44} was {fro[k]} ASVs")
print(f"\nnames gained under the final basis ({len(new)}):")
for k in new: print(f"   + {k:44} now {fin[k]} ASVs")

cat = Counter()
for a in asvs:
    x, y = norm(label_frozen(a)), norm(label_final(a))
    if x == y: cat["identical name"] += 1
    elif x == "unassigned": cat["was unassigned, now resolved"] += 1
    elif y == "unassigned": cat["was resolved, now unassigned"] += 1
    else: cat["resolved to a different name or rank"] += 1
print("\nper-ASV attribution:")
for k, v in cat.most_common(): print(f"   {k:38} {v:>4} ASVs")
print(f"\nyeast-lineage ASVs: {sum(1 for a in asvs if a in yeast)}, labelled identically under both bases.")
print("CONCLUSION: the 192 -> 199 change is NOT caused by the yeast relabelling.")
