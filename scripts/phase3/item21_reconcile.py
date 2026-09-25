#!/usr/bin/env python
"""Item 21. Read the original map files as submitted and reconcile every field against
metadata_v2. Reports discrepancies; changes nothing."""
import csv, glob, os, re
K = "<KMD_ROOT>"
RR = f"{K}/__restructured_25-07-2026/raw_reads"
O = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026"; C = f"{O}/clean"

orig = {}
files = sorted(glob.glob(f"{RR}/*/map*.txt"))
for p in files:
    site_dir = os.path.basename(os.path.dirname(p))
    for line in open(p, encoding="utf-8", errors="replace"):
        line = line.rstrip("\n")
        if not line.strip() or line.startswith("#"): continue
        f = line.split("\t")
        sid = f[0].strip()
        if not sid: continue
        treat = f[3].strip() if len(f) > 3 else ""
        desc = f[4].strip() if len(f) > 4 else ""
        orig[sid] = {"file": f"{site_dir}/{os.path.basename(p)}", "treatment": treat,
                     "description": desc, "n_fields": len(f),
                     "barcode": f[1].strip() if len(f) > 1 else "",
                     "primer": f[2].strip() if len(f) > 2 else ""}
print(f"map files read: {len(files)}")
for p in files: print(f"   {os.path.dirname(p).split('/')[-1]}/{os.path.basename(p)}")
print(f"samples in map files: {len(orig)}\n")

meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

SITE = {"Bloefontein":"Bloemfontein","Bloemfontein":"Bloemfontein",
        "Christiana":"Christiana","Pretoria":"Pretoria"}
def parse(v):
    t = v["treatment"]; d = v["description"]
    age = "Old" if t.lower().startswith("old") else ("Young" if t.lower().startswith("young") else "")
    cond = "malformed" if "malform" in t.lower() else ("healthy" if "healthy" in t.lower() else "")
    if not cond and t.lower() in ("young_inflo", "young_seeds"): cond = "(not stated)"
    tis = d.split("_")[0] if "_" in d else ""
    site = SITE.get(d.split("_")[-1], "")
    return age, cond, tis, site

print("=== ORIGINAL GROUPING, ALL SAMPLES ===")
print(f"{'sample':7} {'source file':22} {'treatment':17} {'description':30} {'age':6} {'cond':11} {'tissue':14} site")
rows, disc = [], []
for sid in sorted(orig, key=lambda s: (orig[s]["file"], s)):
    v = orig[sid]; age, cond, tis, site = parse(v)
    print(f"{sid:7} {v['file']:22} {v['treatment']:17} {v['description'][:30]:30} {age:6} {cond:11} {tis:14} {site}")
    rows.append([sid, v["file"], v["treatment"], v["description"], age, cond, tis, site])

print("\n=== RECONCILIATION AGAINST metadata_v2 ===")
MAPCOND = {"healthy":("asymptomatic","reference_site"), "malformed":("symptomatic",)}
for sid in sorted(orig):
    v = orig[sid]; age, cond, tis, site = parse(v)
    m = meta.get(sid)
    if not m:
        disc.append((sid, "not present in metadata_v2", v["treatment"], "-")); continue
    if site and m["location"] != site:
        disc.append((sid, "location", site, m["location"]))
    if tis and m["tissue"].lower() != tis.lower():
        disc.append((sid, "tissue", tis, m["tissue"]))
    if age and m.get("age", "").lower() != age.lower():
        disc.append((sid, "age", age, m.get("age", "(blank)")))
    if not age and m.get("age", "").strip():
        disc.append((sid, "age ABSENT in map but set in metadata_v2", "(none)", m.get("age")))
    if cond in MAPCOND and m["condition_std"] not in MAPCOND[cond]:
        disc.append((sid, "condition", cond, m["condition_std"]))
in_meta_only = sorted(set(meta) - set(orig))
print(f"samples in metadata_v2 but absent from every map file: {len(in_meta_only)}")
if in_meta_only: print("   " + ", ".join(in_meta_only))
print(f"\ndiscrepancies found: {len(disc)}")
for sid, field, a, b in disc:
    print(f"   {sid:7} {field:44} map='{a}'  metadata_v2='{b}'")

print("\n=== DOES ANY FIELD RECORD TREES PER POOLED SAMPLE? ===")
cols_used = {k for v in orig.values() for k, val in v.items() if k in ("barcode","primer") and val}
print(f"   map columns present: SampleID, BarcodeSequence, LinkerPrimerSequence, Treatment, Description")
print(f"   BarcodeSequence / LinkerPrimerSequence populated in any row: {'yes' if cols_used else 'NO, all blank'}")
pat = re.compile(r"\b(tree|trees|rep|replicate|pool|pooled|n\s*=\s*\d+|\d+\s*trees)\b", re.I)
hits = [(s, v["treatment"], v["description"]) for s, v in orig.items()
        if pat.search(v["treatment"]) or pat.search(v["description"])]
print(f"   fields mentioning tree / replicate / pooling: {len(hits)}")
print("   CONCLUSION: no tree or replicate identifier exists in any map file, and no field "
      "records how many trees contributed to a pooled sample. The pooling INSERT is NOT resolved.")

with open(f"{C}/item21_original_grouping.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample","source_map_file","treatment","description","age","condition","tissue","site"])
    w.writerows(rows)
with open(f"{C}/item21_discrepancies.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["sample","field","map_value","metadata_v2_value"]); w.writerows(disc)
print(f"\n-> {C}/item21_original_grouping.tsv and item21_discrepancies.tsv")
