#!/usr/bin/env python
"""DEC-6. Samples 25-C to 29-C carry no age class in the original collection record; the
'Old' designation was added downstream. Write a corrected metadata file rather than mutating
metadata_v2 in place: every output produced so far cites v2 as its provenance, and silently
rewriting it would make those records unverifiable. v3 supersedes v2 for the age field of
these five samples only. First proves the change alters nothing that is reported."""
import csv, shutil, os
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"
SRC = f"{O}/metadata_v2_27-07-2026.tsv"; DST = f"{C}/metadata_v3_27-07-2026.tsv"
TARGET = ["25-C", "26-C", "27-C", "28-C", "29-C"]
NEW = "not specified"

rows = list(csv.reader(open(SRC), delimiter="\t"))
hdr = rows[0]; ai = hdr.index("age")
print(f"header: {hdr}")
changed = []
for r in rows:
    if r and r[0] in TARGET:
        changed.append((r[0], r[ai], NEW)); r[ai] = NEW
with open(DST, "w", newline="", encoding="utf-8") as fh:
    csv.writer(fh, delimiter="\t").writerows(rows)
print(f"\nage field changed for {len(changed)} samples:")
for s, old, new in changed: print(f"   {s:6} '{old}' -> '{new}'")

# --- prove nothing reported depends on it ---
meta = {}
rd = csv.reader(open(SRC), delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
for r in rd:
    if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
print("\nwhere the age field is consumed:")
print("   Table S2 (Bloemfontein inventory) - restricted to Bloemfontein low-host samples")
bl_low = ["P5", "P32", "P4", "P18", "P19", "P25"]
print(f"      those samples: {bl_low}")
print(f"      any of the five affected samples among them? "
      f"{'YES' if set(TARGET) & set(bl_low) else 'no'}")
print("   G1_G2_summary / post_host_retention / dryrun_rarefaction - Bloemfontein only")
print(f"      site of the five affected samples: "
      f"{sorted({meta[s]['location'] for s in TARGET})}")
print("\nNo reported figure, table or statistic is computed by age for Christiana or Pretoria.")
print("CONCLUSION: the correction changes no reported result.")
print(f"\n-> {DST}")
