#!/usr/bin/env python
"""Item 29. Verify DRAFT_INTRODUCTION.md. READ ONLY; this script never writes to that file.
29b: no number from the present study may appear.
29c: every citation must exist in the consolidated reference list.
29d: zero em dashes."""
import csv, io, re, os
K = "<KMD_ROOT>"
C = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
intro = io.open(f"{K}/DRAFT_INTRODUCTION.md", encoding="utf-8").read()
fin = io.open(f"{K}/__reanalysis_2026-06/phase3/FINAL_NUMBERS_27-07-2026.md", encoding="utf-8").read()

print("=== 29b. numbers from the present study ===")
# every distinctive quantity this study reports
OWN = {"1,142":"corrected fungal ASVs","2,517,177":"corrected fungal reads","57,124":"rarefaction depth",
       "1,173":"pre-screen fungal ASVs","2,523,611":"pre-screen fungal reads","1,322":"total ASVs",
       "4,010,617":"total reads","1,487,006":"host reads","149":"host ASVs","199":"Bloemfontein taxa",
       "71":"yeast lineage ASVs","78":"yeast lineage, superseded","647,397":"yeast lineage reads",
       "37.1":"host percent","84.3":"Bloemfontein host percent","1,152":"intermediate ASV count",
       "6,434":"non-fungal reads removed","31":"non-fungal ASVs removed","896":"securely fungal",
       "277":"sent to NCBI","57,141":"superseded depth","93,085":"UNITE references"}
found = []
for v, what in OWN.items():
    for m in re.finditer(r"(?<![\d,.])" + re.escape(v) + r"(?![\d,.])", intro):
        ctx = " ".join(intro[max(0, m.start()-70):m.start()+70].split())
        found.append((v, what, ctx))
if found:
    for v, what, ctx in found: print(f"   FOUND '{v}' ({what}): ...{ctx}...")
else:
    print("   none. The Introduction contains no quantity from the present study.")

print("\n=== 29c. citations ===")
cit = set()
# Author forms: "X et al., 2012" | "X and Y, 2007" | "X, 2007" | "X et al. (2018)"
# The trailing year list is captured whole so "Crous et al., 2000, 2003" yields both years.
NAME = r"[A-Z][A-Za-zÀ-ſ'\-]+"
AUTH = rf"({NAME}(?:\s+and\s+{NAME})?(?:\s+et\s+al\.)?)"
for m in re.finditer(AUTH + r",?\s*\(?((?:\d{4}[a-z]?)(?:\s*,\s*\d{4}[a-z]?)*)\)?", intro):
    name = " ".join(m.group(1).split())
    if name.split()[0].lower() in ("figure", "table", "section", "province", "the", "in", "and"):
        continue
    for yr in re.findall(r"\d{4}[a-z]?", m.group(2)):
        cit.add((name, yr))
ref = {}
p = f"{C}/reference_list.tsv"
for r in csv.DictReader(open(p), delimiter="\t"):
    ref[r["citation_key"]] = r["status"]
# keep digits: stripping them collapses "Crous et al. 2000" and "Crous et al. 2003"
# onto one dict key and makes one of them look absent
def norm(s): return re.sub(r"[^a-z0-9]", "", s.lower())
keys = {norm(k): k for k in ref}
present, missing = [], []
for name, yr in sorted(cit, key=lambda x: x[0]):
    surname = name.split()[0]
    hit = next((keys[k] for k in keys if norm(surname) in k and yr in keys[k]), None)
    (present if hit else missing).append((f"{name} {yr}", hit))
print(f"   distinct citations detected: {len(cit)}")
print(f"   already in the consolidated reference list: {len(present)}")
for c, h in present: print(f"      ok      {c:34} -> {h} [{ref[h]}]")
print(f"   NOT in the list, must be added and verified: {len(missing)}")
for c, _ in missing: print(f"      MISSING {c}")

print("\n=== 29d. em dashes ===")
n = intro.count("—")
print(f"   em dashes: {n}  {'ok' if n == 0 else 'FAIL'}")
print(f"   en dashes: {intro.count(chr(0x2013))} (not prohibited, reported for information)")

with open(f"{C}/item29_missing_references.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["citation_as_used", "status"])
    for c, _ in missing: w.writerow([c, "NOT VERIFIED, not in consolidated list"])
print(f"\n-> {C}/item29_missing_references.tsv")
