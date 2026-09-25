#!/usr/bin/env python
"""Item 20. Final verification. Checks the manuscript prose against FINAL_NUMBERS, the
corrected GenBank file, and the supplementary tables. Reports mismatches; changes nothing."""
import csv, re, io, os, zipfile
B = "<KMD_ROOT>/__reanalysis_2026-06"
K = "<KMD_ROOT>"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"
fails = []

def rd(p): return io.open(p, encoding="utf-8", errors="replace").read()
meth, res = rd(f"{K}/DRAFT_METHODS.md"), rd(f"{K}/DRAFT_RESULTS.md")
fin = rd(f"{B}/phase3/FINAL_NUMBERS_27-07-2026.md")

print("=== 20a. prose against FINAL_NUMBERS ===")
# values that must appear, and values that must NOT survive anywhere
MUST = [("1,142", "corrected fungal ASV count"), ("2,517,177", "corrected fungal read total"),
        ("57,124", "re-locked rarefaction depth"), ("199", "Bloemfontein taxon count"),
        ("71", "yeast lineage size")]
STALE = [("57,141", "superseded depth"), ("192 taxa", "superseded inventory count"),
         ("542 of 1,173", "superseded guild denominator"),
         ("assigned none of its 78", "superseded yeast lineage size")]
for doc, name in ((meth, "DRAFT_METHODS"), (res, "DRAFT_RESULTS")):
    for v, what in MUST:
        if v not in doc: print(f"   NOTE  {name}: '{v}' ({what}) not present")
    for v, what in STALE:
        n = doc.count(v)
        if n:
            # a superseded value is allowed only where it is explicitly labelled as superseded
            ctx = [m.start() for m in re.finditer(re.escape(v), doc)]
            bad = [i for i in ctx if not re.search(
                r"supersed|before that correction|previously|then included|first delimited|was 27-C",
                doc[max(0, i-360):i+360], re.I)]
            if bad:
                fails.append(f"{name}: unlabelled superseded value '{v}' ({what}) x{len(bad)}")
                print(f"   FAIL  {name}: '{v}' appears unlabelled {len(bad)}x")
            else:
                print(f"   ok    {name}: '{v}' appears {n}x, all labelled as superseded")
for v, what in [("1,142", "ASVs"), ("2,517,177", "reads"), ("57,124", "depth"), ("199", "inventory")]:
    print(f"   {'ok  ' if v in fin else 'FAIL'}  FINAL_NUMBERS contains {v} ({what})")
    if v not in fin: fails.append(f"FINAL_NUMBERS missing {v}")

print("\n=== 20b. corrected GenBank file ===")
gb = f"{C}/genbank_submission/Tremellomycetes_ASVs.fasta"
hdrs = [l for l in rd(gb).splitlines() if l.startswith(">")]
drop = {l.strip() for l in open(f"{O}/full_screen/nonfungal_asvs.txt") if l.strip()}
bad = [h for h in hdrs if (m := re.search(r"md5=([0-9a-f]+)", h)) and m.group(1) in drop]
print(f"   sequences: {len(hdrs)}   (expected 71)")
print(f"   non-fungal sequences remaining: {len(bad)}   (expected 0)")
if len(hdrs) != 71: fails.append(f"GenBank has {len(hdrs)} sequences, expected 71")
if bad: fails.append(f"GenBank still contains {len(bad)} non-fungal sequences")
nseq = sum(1 for l in rd(gb).splitlines() if l and not l.startswith(">"))
print(f"   sequence lines: {nseq} (one per record: {'ok' if nseq == len(hdrs) else 'CHECK'})")

print("\n=== 20c. supplementary tables on the corrected table ===")
tab = [l for l in rd(f"{C}/feature-table-clean.tsv").splitlines() if l.strip()]
n_final = len(tab) - 1
print(f"   corrected feature table: {n_final} ASVs")
checks = [("Table S1 catalogue", f"{C}/asv_catalogue/ASV_catalogue.tsv", None),
          ("Table S2 inventory",  f"{C}/TableS2_bloemfontein_inventory.tsv", None),
          ("Table S3 thresholds", f"{C}/TableS3_taxonomy_thresholds.tsv", str(n_final)),
          ("Table S8 Hill",       f"{C}/TableS8_hill_numbers.tsv", None),
          ("guild table",         f"{C}/guild_assignments_per_ASV.tsv", None)]
for name, p, needle in checks:
    if not os.path.exists(p):
        print(f"   FAIL  {name}: missing"); fails.append(f"{name} missing"); continue
    rows = [l for l in rd(p).splitlines() if l.strip() and not l.startswith("#")]
    extra = ""
    if name == "Table S1 catalogue" or name == "guild table":
        n = len(rows) - 1
        extra = f", {n} data rows {'ok' if n == n_final else 'MISMATCH'}"
        if n != n_final: fails.append(f"{name} has {n} rows, table has {n_final}")
    if needle and needle not in rd(p):
        extra += f", denominator {needle} NOT stated"; fails.append(f"{name} denominator")
    elif needle: extra += f", denominator n={needle} stated"
    print(f"   ok    {name}: {len(rows)} lines{extra}")
sup = [l for l in rd(f"{C}/TableS8_hill_numbers.tsv").splitlines() if "57124" in l]
print(f"   Table S8 rows at depth 57,124: {len(sup)} (expected 48)")
if len(sup) != 48: fails.append("Table S8 depth rows")

print("\n=== SUMMARY ===")
print("ALL CHECKS PASSED" if not fails else f"{len(fails)} PROBLEM(S):")
for f in fails: print(f"   - {f}")
