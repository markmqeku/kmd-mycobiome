#!/usr/bin/env python
"""Item 72. Verify the flagged references against Crossref, targeting the PUBLISHED version and
rejecting review records, preprints and commentaries. A match must satisfy year, a distinctive
title token, AND a record-type guard. Anything short of that is reported unverified with a
statement of what must be checked by hand. Nothing is guessed."""
import csv, json, re, time, urllib.parse, urllib.request
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
UA = "KMD-reference-check (mailto:<email>)"

# key: (search string, expected year, distinctive title token, expected container substring or None)
TARGETS = {
 # GROUP A
 "Chao et al. 2014": ("Rarefaction and extrapolation with Hill numbers a framework for sampling and estimation in species diversity studies", "2014", "hill numbers", "ecological monographs"),
 "Crous et al. 2003": ("Muribasidiospora indica causing a prominent leaf spot disease on Rhus lancea in South Africa", "2003", "muribasidiospora", None),
 "Hoang et al. 2018": ("UFBoot2 Improving the Ultrafast Bootstrap Approximation", "2018", "ufboot", "molecular biology and evolution"),
 "Magurran 2004": ("Measuring Biological Diversity", "2004", "measuring biological diversity", None),
 "Marasas et al. 2006": ("Mango malformation disease and the associated Fusarium species", "2006", "malformation", None),
 "Minh et al. 2020": ("IQ-TREE 2 New Models and Efficient Methods for Phylogenetic Inference in the Genomic Era", "2020", "iq-tree 2", "molecular biology and evolution"),
 "Rognes et al. 2016": ("VSEARCH a versatile open source tool for metagenomics", "2016", "vsearch", "peerj"),
 # GROUP B
 "Bolyen et al. 2019": ("Reproducible interactive scalable and extensible microbiome data science using QIIME 2", "2019", "qiime 2", "nature biotechnology"),
 "Põlme et al. 2020": ("FungalTraits a user-friendly traits database of fungi and fungus-like stramenopiles", "2020", "fungaltraits", "fungal diversity"),
 # GROUP C
 "Coates-Palgrave 2002": ("Keith Coates-Palgrave Trees of Southern Africa", "2002", "trees of southern africa", None),
 "Crous et al. 2000": ("Phytopathogenic fungi from South Africa", "2000", "phytopathogenic fungi", None),
 "Jami et al. 2013": ("Botryosphaeriaceae species overlap on four unrelated native South African hosts", "2013", "botryosphaeriaceae", None),
 "Koekemoer et al. 2013": ("Guide to plant families of southern Africa Strelitzia", "2013", "plant families", None),
 "Krishnan et al. 2009": ("Mango Mangifera indica malformation an unsolved mystery", "2009", "malformation", None),
 "Lange et al. 2012": ("Effects of different soil ameliorants on karee trees Searsia lancea", "2012", "ameliorant", None),
 "Moffett 2007": ("Name changes in the Old World Rhus and recognition of Searsia Anacardiaceae", "2007", "searsia", "bothalia"),
 "Van Wyk & Gericke 2007": ("People's Plants A Guide to Useful Plants of Southern Africa", "2007", "people's plants", None),
}
BAD_TYPE = re.compile(r"peer-review|posted-content|component|grant|dataset", re.I)
BAD_CONTAINER = re.compile(r"faculty opinions|post-publication|peer review|preprint|biorxiv|cabi compendium", re.I)

def query(q, rows=5):
    url = ("https://api.crossref.org/works?rows=%d&select=title,issued,DOI,container-title,"
           "author,type,ISBN,publisher,volume,page&query.bibliographic=" % rows) + urllib.parse.quote(q)
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for _ in range(3):
        try:
            return json.loads(urllib.request.urlopen(req, timeout=50).read().decode()
                              )["message"]["items"]
        except Exception:
            time.sleep(4)
    return []

out = []
for key, (q, want_year, token, want_container) in TARGETS.items():
    items = query(q)
    chosen, why = None, "no candidate satisfied year, title and record type"
    for it in items:
        title = (it.get("title") or [""])[0]
        cont = (it.get("container-title") or [""])[0]
        typ = it.get("type", "")
        dp = it.get("issued", {}).get("date-parts", [[None]])
        yr = str(dp[0][0]) if dp and dp[0] and dp[0][0] else ""
        if BAD_TYPE.search(typ) or BAD_CONTAINER.search(cont): continue
        if token not in title.lower(): continue
        if yr != want_year: continue
        if want_container and want_container not in cont.lower(): continue
        chosen = it; why = ""
        break
    if chosen:
        cont = (chosen.get("container-title") or [""])[0]
        dp = chosen.get("issued", {}).get("date-parts", [[None]])
        out.append([key, "VERIFIED", (chosen.get("title") or [""])[0], str(dp[0][0]),
                    cont, (chosen.get("author") or [{}])[0].get("family", ""),
                    chosen.get("DOI", ""),
                    f"type={chosen.get('type','')}; vol={chosen.get('volume','')}; "
                    f"pages={chosen.get('page','')}; isbn={','.join(chosen.get('ISBN') or []) or 'n/a'}; "
                    f"publisher={chosen.get('publisher','')}"])
    else:
        best = ""
        if items:
            t0 = (items[0].get("title") or [""])[0]
            c0 = (items[0].get("container-title") or [""])[0]
            best = f" Closest Crossref result: '{t0[:70]}' in '{c0[:50]}', type {items[0].get('type','')}"
        out.append([key, "UNVERIFIED", "", "", "", "", "", why + "." + best])
    time.sleep(1)

with open(f"{C}/reference_list.tsv") as fh:
    rows = list(csv.reader(fh, delimiter="\t"))
hdr, body = rows[0], rows[1:]
res = {r[0]: r for r in out}
for r in body:
    if r[0] in res and res[r[0]][1] == "VERIFIED":
        r[1:] = res[r[0]][1:]
    elif r[0] in res:
        r[7] = res[r[0]][7]
with open(f"{C}/reference_list.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)

nv = sum(1 for r in out if r[1] == "VERIFIED")
print(f"attempted {len(out)}; newly VERIFIED {nv}; still unverified {len(out)-nv}\n")
for r in out:
    print(f"{r[1]:11} {r[0]}")
    if r[1] == "VERIFIED":
        print(f"            {r[2][:88]}")
        print(f"            {r[4]}  {r[3]}  doi:{r[6]}")
        print(f"            {r[7]}")
    else:
        print(f"            {r[7][:200]}")
tot = sum(1 for r in body if r[1] == "VERIFIED")
print(f"\nreference list overall: {len(body)} entries, {tot} verified, {len(body)-tot} not")
