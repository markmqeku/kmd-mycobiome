#!/usr/bin/env python
"""20e. Re-run the reference check and report each flagged entry individually: the citation
as used, what Crossref actually returned, and the specific reason the match was rejected.
A returned DOI is not evidence of a correct match; that assumption is what produced the
earlier false 'all verified' result."""
import csv, json, re, time, urllib.parse, urllib.request
A = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/assembly"
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
UA = "KMD-reference-check (mailto:<email>)"
QUERY = {
 "Bolyen et al. 2019": "Reproducible interactive scalable and extensible microbiome data science using QIIME 2",
 "Callahan et al. 2016": "DADA2 High-resolution sample inference from Illumina amplicon data",
 "White et al. 1990": "Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics",
 "Nilsson et al. 2019": "The UNITE database for molecular identification of fungi",
 "Chao et al. 2014": "Rarefaction and extrapolation with Hill numbers a framework for sampling and estimation in species diversity studies",
 "Magurran 2004": "Measuring Biological Diversity",
 "Hsieh et al. 2016": "iNEXT an R package for rarefaction and extrapolation of species diversity Hill numbers",
 "Li 2018": "hillR taxonomic functional and phylogenetic diversity and similarity through Hill Numbers",
 "Polme et al. 2020": "FungalTraits a user-friendly traits database of fungi and fungus-like stramenopiles",
 "Nguyen et al. 2016": "FUNGuild an open annotation tool for parsing fungal community datasets by ecological guild",
 "Katoh & Standley 2013": "MAFFT multiple sequence alignment software version 7 improvements in performance and usability",
 "Minh et al. 2020": "IQ-TREE 2 new models and efficient methods for phylogenetic inference in the genomic era",
 "Hoang et al. 2018": "UFBoot2 improving the ultrafast bootstrap approximation",
 "Martin 2011": "Cutadapt removes adapter sequences from high-throughput sequencing reads",
 "Rognes et al. 2016": "VSEARCH a versatile open source tool for metagenomics",
 "Swanepoel et al. 2018": "karee Searsia lancea malformation",
 "Crous et al. 2003": "Muribasidiospora",
 "Marasas et al. 2006": "mango malformation Fusarium",
 "Schulz & Boyle 2005": "The endophytic continuum",
 "Strobel & Daisy 2003": "Bioprospecting for microbial endophytes and their natural products",
 "Guo et al. 2000": "Identification of endophytic fungi from Livistona chinensis based on morphology and rDNA sequences",
 "Caporaso et al. 2010": "QIIME allows analysis of high-throughput community sequencing data",
}
TOKEN = {
 "Bolyen et al. 2019":"qiime 2","Callahan et al. 2016":"dada2","White et al. 1990":"ribosomal",
 "Nilsson et al. 2019":"unite","Chao et al. 2014":"hill numbers","Magurran 2004":"measuring biological diversity",
 "Hsieh et al. 2016":"inext","Li 2018":"hillr","Polme et al. 2020":"fungaltraits",
 "Nguyen et al. 2016":"funguild","Katoh & Standley 2013":"mafft","Minh et al. 2020":"iq-tree",
 "Hoang et al. 2018":"ufboot","Martin 2011":"cutadapt","Rognes et al. 2016":"vsearch",
 "Swanepoel et al. 2018":"karee","Crous et al. 2003":"muribasidiospora","Marasas et al. 2006":"mango",
 "Schulz & Boyle 2005":"endophytic continuum","Strobel & Daisy 2003":"endophyte",
 # 'sterilia' was the wrong token here: it belongs to a different Guo paper. Crossref was
 # returning the correct New Phytologist record all along and my test rejected it.
 "Guo et al. 2000":"livistona","Caporaso et al. 2010":"qiime",
}
def crossref(q):
    url = ("https://api.crossref.org/works?rows=1&select=title,issued,DOI,container-title,author"
           "&query.bibliographic=" + urllib.parse.quote(q))
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for _ in range(3):
        try:
            d = json.loads(urllib.request.urlopen(req, timeout=45).read().decode())
            it = d.get("message", {}).get("items", [])
            if not it: return None
            r = it[0]
            yr = ""
            dp = r.get("issued", {}).get("date-parts", [[None]])
            if dp and dp[0] and dp[0][0]: yr = str(dp[0][0])
            au = r.get("author", [{}])
            first = (au[0].get("family", "") if au else "")
            return {"title": (r.get("title") or [""])[0], "year": yr, "doi": r.get("DOI", ""),
                    "journal": (r.get("container-title") or [""])[0], "first_author": first}
        except Exception:
            time.sleep(5)
    return None

rows, flagged = [], []
for key, q in QUERY.items():
    want = (re.search(r"(\d{4})", key) or [None, ""])[1] if re.search(r"(\d{4})", key) else ""
    want = re.search(r"(\d{4})", key).group(1)
    r = crossref(q)
    tok = TOKEN.get(key, "").lower()
    if r is None:
        status, why = "NEEDS MANUAL CHECK", "Crossref returned no record for this query"
    else:
        yok, tok_ok = (r["year"] == want), (tok and tok in r["title"].lower())
        if yok and tok_ok: status, why = "VERIFIED", ""
        elif tok_ok and not yok:
            status, why = "FLAGGED", f"title matches but Crossref year is {r['year']}, citation says {want}"
        elif yok and not tok_ok:
            status, why = "FLAGGED", f"year matches but returned title does not contain '{tok}'"
        else:
            status, why = "FLAGGED", f"neither year nor title matched (returned {r['year']}, '{r['title'][:60]}')"
    rows.append([key, status, (r or {}).get("title", ""), (r or {}).get("year", ""),
                 (r or {}).get("journal", ""), (r or {}).get("first_author", ""),
                 (r or {}).get("doi", "") if status == "VERIFIED" else "", why])
    if status != "VERIFIED": flagged.append(rows[-1])
    time.sleep(1)

with open(f"{C}/reference_list.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["citation_key","status","crossref_title","crossref_year","crossref_journal",
                "crossref_first_author","doi_if_verified","reason_if_flagged"])
    w.writerows(rows)
v = sum(1 for r in rows if r[1] == "VERIFIED")
print(f"VERIFIED {v} of {len(rows)}. No DOI is carried for any entry that failed.\n")
print(f"=== {len(flagged)} FLAGGED REFERENCES, individually ===\n")
for r in flagged:
    print(f"{r[0]}")
    print(f"   as used in the manuscript : {QUERY[r[0]]}")
    print(f"   Crossref returned         : {r[2][:88] or '(no record)'}")
    if r[2]:
        print(f"                               {r[5]}, {r[3]}, {r[4][:60]}")
    print(f"   why it failed             : {r[7]}\n")
print(f"-> {C}/reference_list.tsv")
