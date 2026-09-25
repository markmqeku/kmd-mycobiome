#!/usr/bin/env python
"""KMD PROMPT 4, R4. Reads assigned to Didymellaceae and Mycosphaerellaceae (family) and to the genera Didymella
and Mycosphaerella, per site and condition, on the corrected fungal table (1,142 ASVs) with the final
rank-truncated taxonomy (species >= 97, genus >= 95, family >= 90, order >= 85 percent identity). Relative
abundance = reads in the group assigned to the taxon / all fungal reads of the group's libraries (pooled),
the basis used for composition in FINAL_NUMBERS section 7. Library-level shares are also given.
Bloemfontein: presence only (number of libraries with the taxon)."""
import csv, json
from collections import defaultdict
from pathlib import Path

O = Path(r"<KMD_ROOT>\__reanalysis_2026-06\phase3\outputs_27-07-2026")
tax = {r[0]: r[1] for r in csv.reader(open(O / "final_taxonomy.tsv", encoding="utf-8"), delimiter="\t") if len(r) > 1}
meta = {r["sample-id"]: r for r in csv.DictReader(open(O / "clean" / "metadata_v3_27-07-2026.tsv", encoding="utf-8"), delimiter="\t")
        if not r["sample-id"].startswith("#")}
rows = [l.rstrip("\n").split("\t") for l in open(O / "clean" / "feature-table-clean.tsv", encoding="utf-8") if not l.startswith("# ")]
cols = rows[0][1:]
tab = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
missing = [a for a in tab if a not in tax]
assert not missing, f"{len(missing)} ASVs without final taxonomy"

TAXA = {"Didymellaceae (family, all ASVs)": lambda t: "f__Didymellaceae" in t,
        "Didymellaceae, no deeper assignment": lambda t: t.rstrip(";").endswith("f__Didymellaceae"),
        "Didymella (genus)": lambda t: "g__Didymella" in t.split(";")[-1] and t.split(";")[-1] == "g__Didymella" or ";g__Didymella;" in t + ";",
        "Mycosphaerellaceae (family, all ASVs)": lambda t: "f__Mycosphaerellaceae" in t,
        "Mycosphaerellaceae, no deeper assignment": lambda t: t.rstrip(";").endswith("f__Mycosphaerellaceae"),
        "Mycosphaerella (genus)": lambda t: ";g__Mycosphaerella;" in t + ";"}
GROUPS = [("Christiana", "asymptomatic"), ("Christiana", "symptomatic"), ("Pretoria", "asymptomatic"), ("Pretoria", "symptomatic")]
out = {}
for (site, cond) in GROUPS:
    idx = [i for i, c in enumerate(cols) if meta[c]["location"] == site and meta[c]["condition_std"] == cond]
    total = sum(sum(v[i] for i in idx) for v in tab.values())
    for name, f in TAXA.items():
        asvs = [a for a in tab if f(tax[a])]
        reads = sum(sum(tab[a][i] for i in idx) for a in asvs)
        per_lib = []
        for i in idx:
            lt = sum(v[i] for v in tab.values())
            per_lib.append(100 * sum(tab[a][i] for a in asvs) / lt if lt else 0)
        out[(site, cond, name)] = dict(n=len(idx), group_reads=int(total), taxon_reads=int(reads), asvs_with_reads=sum(1 for a in asvs if sum(tab[a][i] for i in idx) > 0),
                                      pooled_pct=100 * reads / total, lib_min=min(per_lib), lib_max=max(per_lib))
b_idx = [i for i, c in enumerate(cols) if meta[c]["location"] == "Bloemfontein"]
bloem = {name: sum(1 for i in b_idx if sum(tab[a][i] for a in tab if f(tax[a])) > 0) for name, f in TAXA.items()}

print(f"{'group':28} {'taxon':42} {'n':>2} {'reads':>7} {'of':>8} {'pooled %':>9} {'library range %':>18} ASVs")
for (site, cond, name), v in out.items():
    print(f"{site+' '+cond:28} {name:42} {v['n']:>2} {v['taxon_reads']:>7,} {v['group_reads']:>8,} {v['pooled_pct']:>8.2f}  {v['lib_min']:>6.2f} to {v['lib_max']:<6.2f} {v['asvs_with_reads']}")
print(f"\nBloemfontein, presence in {len(b_idx)} libraries with fungal reads (P2 has none, 29 usable):", bloem)
for name, f in TAXA.items():
    if "genus" in name:
        print(f"  ASVs labelled {name}: {[(a[:8], tax[a].split(';')[-1]) for a in tab if f(tax[a])]}")
json.dump({"groups": {f"{k[0]}|{k[1]}|{k[2]}": v for k, v in out.items()}, "bloemfontein_presence": bloem},
          open(O.parent / "r4_didymella.json", "w", encoding="utf-8"), indent=1)
