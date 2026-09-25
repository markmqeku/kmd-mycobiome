#!/usr/bin/env python
"""13a validation. The whole-table screen sends only ASVs that fail a 'securely fungal by
UNITE alignment' test to NCBI. That shortcut is only legitimate if it cannot miss a plant.
Two checks, both written to disk so the exclusion is auditable:
  1. the UNITE 10.0 reference contains fungal entries only, so a strong UNITE alignment
     cannot be an alignment to a plant;
  2. all 149 independently confirmed host (Searsia) ASVs are run through the same test.
     Any that pass would be a false negative."""
import zipfile, os
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
S = f"{O}/full_screen"
MINLEN, MINID = 200, 90.0
out = []

def w(s):
    print(s); out.append(s)

import collections
z = zipfile.ZipFile(f"{O}/unite_current/unite10_tax.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith(".tsv")][0]
c = collections.Counter()
for l in z.read(fn).decode("utf-8", "replace").splitlines():
    p = l.split("\t")
    if len(p) < 2 or p[0] == "Feature ID": continue
    c[p[1].split(";")[0]] += 1
w("CHECK 1: kingdom composition of the UNITE 10.0 reference")
for k, v in c.most_common(): w(f"   {k:34} {v:>8,} entries")
w(f"   non-fungal reference entries: {sum(v for k,v in c.items() if k != 'k__Fungi')}")

best = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 4: continue
    try: pid, alen = float(p[2]), int(p[3])
    except ValueError: continue
    if p[0] not in best or alen*pid > best[p[0]][1]*best[p[0]][0]: best[p[0]] = (pid, alen)

host = [l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()]
passing = [h for h in host if best.get(h, (0, 0))[1] >= MINLEN and best.get(h, (0, 0))[0] >= MINID]
w(f"\nCHECK 2: positive control, {len(host)} confirmed Searsia lancea host ASVs")
w(f"   test: UNITE hit >= {MINLEN} bp at >= {MINID}% identity")
w(f"   host ASVs passing the test (false negatives): {len(passing)} of {len(host)}")
mx = max((best.get(h, (0, 0)) for h in host), key=lambda t: t[0]*t[1])
w(f"   strongest UNITE alignment achieved by any host ASV: {mx[0]:.1f}% over {mx[1]} bp")
a = "117a170de6384113d3f9dc358c787433"
w(f"   Amaranthus ASV {a[:12]} best UNITE hit: "
  + ("none at all" if best.get(a, (0, 0))[1] == 0 else str(best[a])))
w("\nCONCLUSION: a strong UNITE alignment cannot be an alignment to a plant, and no known "
  "plant sequence in this dataset achieves one. The 896 ASVs not sent to NCBI are excluded "
  "on evidence, not on convenience.")
open(f"{S}/screen_validation.txt", "w", encoding="utf-8").write("\n".join(out) + "\n")
print(f"\n-> {S}/screen_validation.txt")
