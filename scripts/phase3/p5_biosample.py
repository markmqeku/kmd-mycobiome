#!/usr/bin/env python
"""KMD PROMPT 5, D1: apply the verified deposition terms to the BioSample batch file and the SRA metadata for the 45
usable libraries of this study. Terms were verified in OLS and NCBI Taxonomy by p5_verify_terms.py. Re-runnable from the
pre-Prompt-5 copies in _retired."""
import csv, json, shutil
from pathlib import Path

P3 = Path(__file__).resolve().parent
ROOT = P3.parents[1]
MD = ROOT / "NCBI_SUBMISSION" / "SRA_metadata"
READS = ROOT / "NCBI_SUBMISSION" / "SRA_reads"
BS, SRA = MD / "BioSample_MIMARKS_survey_plant-associated_6.0.tsv", MD / "SRA_metadata.tsv"
V = json.load(open(P3 / "p5_verified_terms.json", encoding="utf-8"))
assert all(v["ok"] for v in V["ols"].values())
assert V["taxonomy"]["Searsia lancea"][0]["taxid"] == "298678"
for f in (BS, SRA):
    b = ROOT / "_retired" / (f.stem + "_pre-prompt5_backup.tsv")
    if not b.exists(): shutil.copy2(f, b)
S12 = {r["library"]: r for r in csv.DictReader(open(P3 / "outputs_27-07-2026" / "clean" / "TableS12_library_metadata.tsv",
                                                      encoding="utf-8"), delimiter="\t")}


def term(cid):
    return f"{V['ols'][cid]['label']} [{cid}]"


SITE = {"Bloemfontein": {"broad": term("ENVO:01000177"), "local": term("ENVO:00000467"), "pool": "single tree"},
        "Christiana": {"broad": term("ENVO:01000178"), "local": term("ENVO:01000447"),
                       "pool": "tissue from four trees of the same site and condition pooled per library"},
        "Pretoria": {"broad": term("ENVO:01000178"), "local": term("ENVO:00000078"),
                     "pool": "tissue from four trees of the same site and condition pooled per library"}}
STRUC = {"Leaves": term("PO:0025034"), "Twigs": term("PO:0025073"), "Inflorescence": term("PO:0009049"),
         "Seeds": term("PO:0009010"), "Malformation": term("PO:0009006")}
TISSUE_WORD = {"Leaves": "leaf", "Twigs": "twig", "Inflorescence": "inflorescence", "Seeds": "seed",
               "Malformation": "malformed shoot"}
DISEASE = {"symptomatic": "karee malformation disease, symptomatic tissue",
           "asymptomatic": "karee malformation disease, asymptomatic tissue at affected site",
           "reference_site": "reference site, no recorded karee malformation disease"}
OTHER = {"CR", "PR", "S1", "S2", "S3", "S4"}

rows = list(csv.DictReader(open(ROOT / "_retired" / "BioSample_MIMARKS_survey_plant-associated_6.0_pre-prompt5_backup.tsv",
                                encoding="utf-8"), delimiter="\t"))
fields = list(rows[0].keys())
for extra in ("host_sex", "samp_pooling"):
    if extra not in fields: fields.append(extra)
empty = []
for r in rows:
    lib = r["*sample_name"]
    assert lib not in OTHER and lib in S12, lib
    site = S12[lib]["site"]; cond = r["condition"]; tis = r["tissue"]
    r["*env_broad_scale"] = SITE[site]["broad"]
    r["*env_local_scale"] = SITE[site]["local"]
    r["*env_medium"] = term("ENVO:01001121")
    r["plant_struc"] = STRUC[tis]
    r["host_disease"] = DISEASE[cond]
    r["host_taxid"] = "298678"
    sex = S12[lib]["tree sex"]
    r["host_sex"] = sex if sex in ("male", "female") else ""
    r["samp_pooling"] = SITE[site]["pool"]
    cond_word = {"symptomatic": "symptomatic tissue", "asymptomatic": "asymptomatic tissue",
                 "reference_site": "reference site"}[cond]
    unit = "single-tree library" if site == "Bloemfontein" else "library pooling four trees"
    r["sample_title"] = f"ITS2 {unit} {lib}: Searsia lancea {TISSUE_WORD[tis]}, {site}, {cond_word}"
    for k in fields:
        if k.startswith("*") and not r.get(k, "").strip(): empty.append((lib, k))
    assert "MARK_TO_CHOOSE" not in "\t".join(r.get(k, "") for k in fields), lib
with open(BS, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)

# SRA metadata: library wording and pooling; filenames against the 90 files
srows = list(csv.DictReader(open(ROOT / "_retired" / "SRA_metadata_pre-prompt5_backup.tsv", encoding="utf-8"), delimiter="\t"))
sfields = list(srows[0].keys())
for r in srows:
    lib = r["sample_name"]; site = S12[lib]["site"]
    pool = ("tissue from four trees of the same site and condition pooled per library" if site != "Bloemfontein"
            else "single tree per library")
    if "pooled" not in r["design_description"] and "single tree" not in r["design_description"]:
        r["design_description"] = f"{pool}; " + r["design_description"]
    r["title"] = r["title"].replace("ITS2 amplicon of", "ITS2 amplicon library of")
with open(SRA, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=sfields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(srows)
named = sorted([r["filename"] for r in srows] + [r["filename2"] for r in srows])
on_disk = sorted(p.name for p in READS.glob("*.fastq.gz"))
print(f"BioSample rows {len(rows)}; empty mandatory fields {empty or 'none'}")
print(f"SRA rows {len(srows)}; filenames named {len(named)}, files on disk {len(on_disk)}, identical sets {named == on_disk}")
print("sites:", {s: sum(1 for r in rows if S12[r['*sample_name']]['site'] == s) for s in SITE},
      "| host_sex filled:", sum(1 for r in rows if r["host_sex"]))
