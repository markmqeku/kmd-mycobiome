#!/usr/bin/env python
"""KMD PROMPT 5, D0c: Figure 8 tip labels. The reference tips in the neighbourhood trees carry names built from the
GenBank title with the genus dropped ("REF_tinkyukukuisolate_PZ779985.1"). This fetches the organism name of every
reference accession in the eight trees from NCBI (E-utilities esummary, db=nuccore) and writes
placement/reference_organism_names.tsv, which make_fig8_fig5.py uses for the tip labels. Nothing else is changed."""
import csv, json, re, time, urllib.parse, urllib.request
from pathlib import Path
try:
    import truststore; truststore.inject_into_ssl()
except ImportError:
    pass

P = Path(__file__).resolve().parent / "outputs_27-07-2026" / "placement"
accs = set()
for t in P.glob("*_tree.contree"):
    accs |= set(re.findall(r"REF_[^_,:()]+_([A-Z]{1,2}_?\d{5,9}\.\d)", t.read_text()))
accs = sorted(accs)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=nuccore&retmode=json&id="
names = {}
for i in range(0, len(accs), 100):
    chunk = accs[i:i + 100]
    d = json.loads(urllib.request.urlopen(urllib.request.Request(E + ",".join(chunk), headers={"User-Agent": "KMD-verify/1.0"}),
                                          timeout=60).read())["result"]
    for uid in d.get("uids", []):
        r = d[uid]
        names[r.get("accessionversion")] = r.get("organism", "")
    time.sleep(0.4)
with open(P / "reference_organism_names.tsv", "w", encoding="utf-8", newline="") as fh:
    w = csv.writer(fh, delimiter="\t", lineterminator="\n"); w.writerow(["accession", "organism"])
    for a in accs: w.writerow([a, names.get(a, "")])
print(f"{len(accs)} reference accessions; organism found for {sum(1 for a in accs if names.get(a))}; missing "
      f"{[a for a in accs if not names.get(a)]}")
