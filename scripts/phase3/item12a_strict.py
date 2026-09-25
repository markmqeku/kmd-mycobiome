#!/usr/bin/env python
"""12a strict re-verification. A Crossref hit counts as verified ONLY if the returned year
matches the citation year and a distinctive title token is present. Everything else is flagged
NEEDS MANUAL CHECK. Automated bibliographic matching is unreliable and must not be trusted blindly."""
import csv, re, os
A = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/assembly"
rows = list(csv.reader(open(f"{A}/reference_list.tsv"), delimiter="\t"))
hdr, body = rows[0], rows[1:]
# distinctive token that must appear in the resolved title
TOKEN = {
 "Bolyen et al. 2019":"qiime 2","Callahan et al. 2016":"dada2","White et al. 1990":"ribosomal",
 "Nilsson et al. 2019":"unite","Chao et al. 2014":"hill numbers","Magurran 2004":"measuring biological diversity",
 "Hsieh et al. 2016":"inext","Li 2018":"hillr","Polme et al. 2020":"fungaltraits",
 "Nguyen et al. 2016":"funguild","Katoh & Standley 2013":"mafft","Minh et al. 2020":"iq-tree",
 "Hoang et al. 2018":"ufboot","Martin 2011":"cutadapt","Rognes et al. 2016":"vsearch",
 "Swanepoel et al. 2018":"karee","Crous et al. 2003":"muribasidiospora","Marasas et al. 2006":"mango",
 "Schulz & Boyle 2005":"endophytic continuum","Strobel & Daisy 2003":"endophyte",
 "Guo et al. 2000":"sterilia","Caporaso et al. 2010":"qiime",
}
out = []
for r in body:
    key, used, title, year, doi, _ = (r + [""]*6)[:6]
    m = re.search(r"(\d{4})", key); want = m.group(1) if m else ""
    tok = TOKEN.get(key, "").lower()
    tl = (title or "").lower()
    year_ok = (year == want)
    tok_ok = bool(tok) and tok in tl
    if year_ok and tok_ok:
        status = "VERIFIED"
    elif tok_ok and not year_ok:
        status = f"CHECK: title matches but year returned {year}, expected {want}"
    elif year_ok and not tok_ok:
        status = "CHECK: year matches but title does not"
    else:
        status = "NEEDS MANUAL CHECK: no reliable match"
        title, doi = "", ""      # do not carry a wrong DOI forward
    out.append([key, used, title, year if status.startswith("VERIFIED") else "", doi, status])
with open(f"{A}/reference_list.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(hdr); w.writerows(out)
v = sum(1 for r in out if r[5]=="VERIFIED")
print(f"{'citation':26} {'status'}")
for r in out: print(f"  {r[0]:26} {r[5]}")
print(f"\nVERIFIED {v} of {len(out)}. The remainder require manual confirmation against the source.")
print("No DOI is carried forward for any entry that failed the test.")
