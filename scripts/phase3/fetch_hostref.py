#!/usr/bin/env python
"""Fetch host-family (Anacardiaceae / Searsia / Rhus) ITS references from NCBI for the non-target screen."""
import urllib.parse, urllib.request, time, os, sys

D = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen/hostref"
os.makedirs(D, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "KMD-reanalysis/1.0"})
    with urllib.request.urlopen(req, timeout=60) as r:
        return r.read().decode("utf-8", "replace")

QUERIES = [
    ("Searsia",       'Searsia[Organism] AND (internal transcribed spacer[All Fields] OR ITS[All Fields])'),
    ("Rhus",          'Rhus[Organism] AND (internal transcribed spacer[All Fields] OR ITS[All Fields])'),
    ("Anacardiaceae", 'Anacardiaceae[Organism] AND internal transcribed spacer[All Fields]'),
]
out = os.path.join(D, "host_ref.fasta")
open(out, "w").close()
total = 0
for label, q in QUERIES:
    url = f"{E}/esearch.fcgi?db=nuccore&retmax=150&term=" + urllib.parse.quote(q)
    try:
        xml = get(url)
    except Exception as e:
        print(f"  {label}: esearch FAILED {e}"); continue
    ids = [x.split("<")[0] for x in xml.split("<Id>")[1:]]
    print(f"  {label}: {len(ids)} ids")
    if not ids: continue
    for i in range(0, len(ids), 100):
        chunk = ",".join(ids[i:i+100])
        try:
            fa = get(f"{E}/efetch.fcgi?db=nuccore&id={chunk}&rettype=fasta&retmode=text")
            with open(out, "a") as fh: fh.write(fa)
            total += fa.count(">")
        except Exception as e:
            print(f"    efetch failed: {e}")
        time.sleep(0.5)
print(f"host reference sequences written: {total} -> {out}")
