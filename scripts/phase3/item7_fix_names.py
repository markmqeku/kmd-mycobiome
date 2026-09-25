#!/usr/bin/env python
"""Repair the nearest-reference names in placement_results.tsv.
The earlier label builder dropped the genus. Re-fetch the authoritative title per accession
from NCBI and rewrite the name column. Identities, alignments and supports are unchanged."""
import csv, re, json, urllib.request, urllib.parse, time, os
P = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/placement"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
rows = list(csv.reader(open(f"{P}/placement_results.tsv"), delimiter="\t"))
hdr, body = rows[0], rows[1:]
accs = sorted({m.group(0) for r in body for m in [re.search(r"[A-Z]{1,2}\d{5,6}\.\d", r[2])] if m})
print(f"accessions to resolve: {len(accs)}")
titles = {}
for i in range(0, len(accs), 60):
    chunk = accs[i:i+60]
    try:
        d = json.loads(urllib.request.urlopen(
            urllib.request.Request(f"{E}/esummary.fcgi?db=nuccore&id={urllib.parse.quote(','.join(chunk))}&retmode=json",
                                   headers={"User-Agent":"KMD/1.0"}), timeout=60).read().decode())
        for k, v in d.get("result", {}).items():
            if k == "uids": continue
            titles[v.get("accessionversion","")] = v.get("title","")
    except Exception as e:
        print("  esummary failed:", type(e).__name__)
    time.sleep(0.4)

def clean(title, acc):
    if not title: return acc
    t = re.split(r"\b(isolate|strain|culture|voucher|clone|internal transcribed|18S|small subunit|28S)\b",
                 title, maxsplit=1)[0].strip(" ,;")
    parts = t.split()
    name = " ".join(parts[:2]) if len(parts) >= 2 else t
    return f"{name} ({acc})"

fixed = 0
for r in body:
    m = re.search(r"[A-Z]{1,2}\d{5,6}\.\d", r[2])
    if not m: continue
    acc = m.group(0)
    new = clean(titles.get(acc, ""), acc)
    if new and new != r[2]:
        r[2] = new; fixed += 1
    # sisters column
    def fix_s(s):
        out = []
        for part in s.split(";"):
            mm = re.search(r"[A-Z]{1,2}\d{5,6}\.\d", part)
            out.append(clean(titles.get(mm.group(0), ""), mm.group(0)) if mm else part.strip())
        return "; ".join(out)
    if len(r) >= 10 and r[9] and r[9] != "none resolved":
        r[9] = fix_s(r[9])
        r[8] = re.sub(r"sister to .*$", f"sister to {r[9]}", r[8])
with open(f"{P}/placement_results.tsv","w",newline="",encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)
print(f"repaired {fixed} nearest-reference names")
print(f"\n{'neighbourhood':20} {'ASV':22} {'nearest reference':38} {'id%':>6} {'aln':>5} {'supp':>5}")
for r in body:
    print(f"{r[0]:20} {r[1][:22]:22} {r[2][:38]:38} {r[3]:>6} {r[4]:>5} {r[7] or '-':>5}")
