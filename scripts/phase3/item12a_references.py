#!/usr/bin/env python
"""12a: consolidated reference list from DRAFT_METHODS, DRAFT_RESULTS and the 6f harvest,
verified against Crossref where a DOI can be resolved. Nothing invented; unverifiable entries flagged.
Also rewrites SRA_metadata.tsv excluding the four write-off samples whose raw files are unusable."""
import csv, json, os, urllib.parse, urllib.request, time
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; A = f"{O}/assembly"; os.makedirs(A, exist_ok=True)

# (citation key, search string, where it is used)
REFS = [
 ("Bolyen et al. 2019","Reproducible interactive scalable and extensible microbiome data science using QIIME 2","Methods 2.3"),
 ("Callahan et al. 2016","DADA2 High-resolution sample inference from Illumina amplicon data","Methods 2.3"),
 ("White et al. 1990","Amplification and direct sequencing of fungal ribosomal RNA genes for phylogenetics","Methods 2.2"),
 ("Nilsson et al. 2019","The UNITE database for molecular identification of fungi handling dark taxa","Methods 2.5"),
 ("Chao et al. 2014","Rarefaction and extrapolation with Hill numbers a framework for sampling and estimation","Methods 2.7"),
 ("Magurran 2004","Measuring Biological Diversity","Methods 2.7"),
 ("Hsieh et al. 2016","iNEXT an R package for rarefaction and extrapolation of species diversity","Methods 2.7"),
 ("Li 2018","hillR taxonomic functional and phylogenetic diversity and similarity through Hill numbers","Methods 2.7"),
 ("Polme et al. 2020","FungalTraits a user-friendly traits database of fungi and fungus-like stramenopiles","Results 3.9"),
 ("Nguyen et al. 2016","FUNGuild an open annotation tool for parsing fungal community datasets by ecological guild","Results 3.9"),
 ("Katoh & Standley 2013","MAFFT multiple sequence alignment software version 7 improvements in performance and usability","Methods 2.6"),
 ("Minh et al. 2020","IQ-TREE 2 new models and efficient methods for phylogenetic inference in the genomic era","Methods 2.6 / item 7"),
 ("Hoang et al. 2018","UFBoot2 improving the ultrafast bootstrap approximation","Methods 2.6 / item 7"),
 ("Martin 2011","Cutadapt removes adapter sequences from high-throughput sequencing reads","Methods 2.3"),
 ("Rognes et al. 2016","VSEARCH a versatile open source tool for metagenomics","Methods 2.5"),
 ("Swanepoel et al. 2018","First report of karee malformation disease and investigations into its possible cause","Introduction, background"),
 ("Crous et al. 2003","Muribasidiospora indica causing a prominent leaf spot disease on Rhus lancea in South Africa","Introduction"),
 ("Marasas et al. 2006","Mango malformation disease and the associated Fusarium species","Introduction"),
 ("Schulz & Boyle 2005","The endophytic continuum","Introduction"),
 ("Strobel & Daisy 2003","Bioprospecting for microbial endophytes and their natural products","Introduction"),
 ("Guo et al. 2000","Identification of non-sporulating endophytic fungi mycelia sterilia by ribosomal DNA sequences","Introduction"),
 ("Caporaso et al. 2010","QIIME allows analysis of high-throughput community sequencing data","superseded methods, cite only if QIIME1 history is discussed"),
]
def crossref(q):
    try:
        u = "https://api.crossref.org/works?rows=1&query.bibliographic=" + urllib.parse.quote(q)
        d = json.loads(urllib.request.urlopen(
            urllib.request.Request(u, headers={"User-Agent":"KMD-reanalysis/1.0 (mailto:<email>)"}),
            timeout=45).read().decode())
        it = d["message"]["items"]
        if not it: return None
        i = it[0]
        yr = (i.get("issued",{}).get("date-parts",[[None]])[0] or [None])[0]
        return {"title": (i.get("title") or [""])[0], "year": yr, "doi": i.get("DOI",""),
                "container": (i.get("container-title") or [""])[0], "type": i.get("type","")}
    except Exception as e:
        return {"error": type(e).__name__}

rows = []
for key, q, used in REFS:
    r = crossref(q); time.sleep(0.35)
    if not r or "error" in r:
        rows.append([key, used, "", "", "", "NOT VERIFIED" + (f" ({r['error']})" if r and 'error' in r else " (no Crossref match)")])
        print(f"  {key:26} NOT VERIFIED")
        continue
    ok = "verified" if r["doi"] else "no DOI returned"
    rows.append([key, used, r["title"][:110], str(r["year"] or ""), r["doi"], ok])
    print(f"  {key:26} {r['year']}  {r['doi'] or 'no DOI'}")
with open(f"{A}/reference_list.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["citation_key","used_in","resolved_title","year","doi","verification"])
    w.writerows(rows)
print(f"\nwrote {A}/reference_list.tsv  ({len(rows)} entries; "
      f"{sum(1 for r in rows if r[5]=='verified')} verified)")

# ---- SRA: exclude the four write-offs ----
WRITEOFF = {"P6","P15","P24","YT-P"}
p = f"{A}/SRA_metadata.tsv"
if os.path.exists(p):
    rows2 = list(csv.reader(open(p, encoding="utf-8"), delimiter="\t"))
    hdr, body = rows2[0], [r for r in rows2[1:] if r and r[0] not in WRITEOFF]
    with open(p,"w",newline="",encoding="utf-8") as fh:
        w=csv.writer(fh,delimiter="\t"); w.writerow(hdr); w.writerows(body)
    print(f"SRA_metadata.tsv rewritten: {len(body)} samples "
          f"(excluded {sorted(WRITEOFF)}, whose raw files are truncated or unpaired)")
