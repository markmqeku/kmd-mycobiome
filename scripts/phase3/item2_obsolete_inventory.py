#!/usr/bin/env python
"""Item 2: inventory every value in the two superseded drafts now contradicted by FINAL_NUMBERS.
Reads the ACCEPTED-changes view (w:t only; w:delText excluded), so this is what each draft would say
if its tracked changes were accepted. Report only - the drafts are not modified."""
import zipfile, re, sys
from xml.etree import ElementTree as ET
NS = '{http://schemas.openxmlformats.org/wordprocessingml/2006/main}'

DRAFTS = [
    ("A1", r"<KMD_ROOT>/__work_A1_2026-06-17/Article1_SUBMISSION_DRAFT_TRACKED.docx"),
    ("A2", r"<KMD_ROOT>/__work_A2_2026-06-17/Article2_Combined_TRACKED.docx"),
]

def paras(path):
    xml = zipfile.ZipFile(path).read('word/document.xml').decode('utf-8', 'replace')
    root = ET.fromstring(xml)
    out = []
    for p in root.iter(NS + 'p'):
        out.append("".join(t.text or "" for t in p.iter(NS + 't')))
    return out

# pattern -> (label, superseding value per FINAL_NUMBERS)
NUM = [
    (r"\b258\b",            "258 ASVs (Bloemfontein, old pipeline)",      "1,173 fungal ASVs total; Bloemfontein excluded from quantitative analysis"),
    (r"\b316\b",            "316 ASVs (Christiana, old pipeline)",        "superseded; single merged fungi-only table = 1,173 ASVs"),
    (r"\b465\b",            "465 ASVs (Pretoria, old pipeline)",          "superseded; 1,173 fungal ASVs"),
    (r"\b491\b",            "491 ASVs (Pretoria, unsupported)",           "superseded; 1,173 fungal ASVs"),
    (r"477,?786",           "477,786 reads (mislabelled 'filtered')",     "4,010,617 total reads; 2,523,611 fungal after host removal"),
    (r"1,?027,?325",        "1,027,325 raw reads (Bloemfontein)",         "4,010,617 total across 45 samples"),
    (r"1,?545,?777",        "1,545,777 (Pretoria filtered)",              "superseded"),
    (r"745,?249",           "745,249 (Christiana filtered)",              "superseded"),
    (r"\b114\b",            "rarefaction depth 114",                     "depth 57,141 (deepest retaining all 16 core samples)"),
    (r"\b81\s?%",           "81% ASV sharing",                            "effort-corrected: Christiana 143.0 vs 181.1; Pretoria 210.0 vs 165.6"),
    (r"\b82\s?%",           "82% ASV sharing",                            "as above; naive contrast was an effort artifact"),
    (r"\b83\s?%",           "83% ASV sharing",                            "as above"),
    (r"\b43\s?%",           "43% healthy ASV share",                      "as above"),
    (r"\b46\s?%",           "46% healthy ASV share",                      "as above"),
    (r"\b48\s?%",           "48% healthy ASV share",                      "as above"),
    (r"\b259\b",            "259 ASVs malformed (Christiana)",            "superseded by effort-corrected sharing"),
    (r"\b257\b",            "257 ASVs malformed (Christiana, recount)",   "superseded"),
    (r"\b153\b",            "153 ASVs healthy (Christiana)",              "superseded"),
    (r"\b151\b",            "151 ASVs healthy (Christiana, recount)",     "superseded"),
    (r"\b406\b",            "406 ASVs malformed (Pretoria)",              "superseded"),
    (r"\b378\b",            "378 ASVs malformed (Pretoria, recount)",     "superseded"),
    (r"\b213\b",            "213 ASVs healthy (Pretoria)",                "superseded"),
    (r"\b354\b|\b358\b|\b338\b", "mean read length (354/358/338 bp)",     "not carried forward; per-run truncation 228/200 and 280/230"),
]
PVAL = re.compile(r"p\s*[=<>]\s*0?\.?\d+(?:e-?\d+)?|p\s*=\s*\d+e-\d+", re.I)
FLAGS = [
    (re.compile(r"\bhealthy\b", re.I),        "'healthy' used for tissue from symptomatic trees", "condition_std: asymptomatic / symptomatic / reference_site"),
    (re.compile(r"metagenomic", re.I),        "'metagenomic' misnomer",                            "amplicon metabarcoding"),
    (re.compile(r"emergen\w*", re.I),         "'emergence' applied to a taxon",                    "Curvibasidium is a resident that is enriched; direction consistent, magnitude sample-driven"),
    (re.compile(r"Didymella", re.I),          "Didymella abundance claim",                         "NOT REPRODUCIBLE: 0.00-0.01% in all four groups; resolves only to Didymellaceae"),
    (re.compile(r"Mycosphaerella", re.I),     "Mycosphaerella abundance claim",                    "NOT REPRODUCIBLE: 0.00% in all four groups at every threshold"),
    (re.compile(r"(elevat|increas|greater|higher|dramatic)\w*\s+(?:\w+\s+){0,3}divers", re.I),
                                              "claim that diversity is elevated in malformed tissue", "direction inconsistent between sites; Pretoria reverses; no test performed"),
    (re.compile(r"dysbiosis", re.I),          "dysbiosis framing",                                 "not supported; causality unestablished"),
    (re.compile(r"UNITE", re.I),              "UNITE citation",                                    "UNITE release 10.0 (own DOI); Nilsson et al. 2019 as description paper; v8.0 never used"),
]

def clip(s, n=150):
    s = re.sub(r"\s+", " ", s).strip()
    return s[:n] + ("..." if len(s) > n else "")

total = 0
for tag, path in DRAFTS:
    print("=" * 100)
    print(f"### {tag}: {path.split('/')[-1]}")
    print("=" * 100)
    ps = paras(path)
    for i, para in enumerate(ps, 1):
        if not para.strip(): continue
        hits = []
        for pat, label, sup in NUM:
            if re.search(pat, para): hits.append(("VALUE", label, sup))
        for pv in set(PVAL.findall(para)):
            hits.append(("P-VALUE", f"'{pv.strip()}'", "no inferential testing performed; all p-values withdrawn"))
        for rx, label, sup in FLAGS:
            if rx.search(para): hits.append(("FRAMING", label, sup))
        if hits:
            print(f"\n[{tag} para {i}] {clip(para, 220)}")
            for kind, label, sup in hits:
                total += 1
                print(f"    - {kind:8} {label}")
                print(f"      -> superseded by: {sup}")
print("\n" + "=" * 100)
print(f"TOTAL flagged occurrences: {total}")
print("Report only. Neither draft was modified.")
