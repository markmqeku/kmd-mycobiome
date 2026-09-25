#!/usr/bin/env python
"""KMD PROMPT 5, D2: the GenBank file for all 71 lineage ASVs, rRNA-ITS Eukaryote workflow, following the GenBank
help desk guidance (case CAS-1702031-X8R7H0): unique identifiers in "clone", never "isolate".

Two tests per ASV: (1) a blastn hit to the Tremellomycetes reference set (yeast_placement/yeast_vs_tremello.tsv, the
former 68); (2) inside the clade holding all Tremellomycetes references in the rerun tree (fig3_71_placement.tsv, the 70).
Both: "uncultured Tremellomycetes" (taxid 543743). One only: "uncultured fungus" (taxid 175245). Taxids verified in NCBI
Taxonomy by p5_verify_terms.py. Clone names KMD_ITS2_ASV001 to 071 by total read count (fungi-only table), ties by hash.
collection_date: the ISO 8601 date, or interval, spanned by the libraries in which the ASV occurs."""
import csv, hashlib, json, shutil
from pathlib import Path

P3 = Path(__file__).resolve().parent
ROOT = P3.parents[1]
O = P3 / "outputs_27-07-2026"
GB = ROOT / "NCBI_SUBMISSION" / "GenBank"
V = json.load(open(P3 / "p5_verified_terms.json", encoding="utf-8"))["taxonomy"]
TAX = {"uncultured Tremellomycetes": V["uncultured Tremellomycetes"][0]["taxid"],
       "uncultured fungus": V["uncultured fungus"][0]["taxid"]}
assert TAX == {"uncultured Tremellomycetes": "543743", "uncultured fungus": "175245"}


def rd(p):
    s, k = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): k = l[1:].split()[0]; s[k] = ""
        elif k: s[k] += l.strip()
    return s


seqs = [v.upper() for v in rd(O / "clean" / "genbank_submission" / "Tremellomycetes_ASVs.fasta").values()]
assert len(seqs) == 71
md5 = {hashlib.md5(s.encode()).hexdigest(): s for s in seqs}
# feature table: reads and occurrence per library
rows = [r for r in csv.reader(open(O / "clean" / "feature-table-clean.tsv"), delimiter="\t") if r]
head = next(r for r in rows if r[0].startswith("#OTU"))
tab = {r[0]: [float(x) for x in r[1:]] for r in rows if not r[0].startswith("#")}
libs = head[1:]
meta = {r["library"]: r for r in csv.DictReader(open(O / "clean" / "TableS12_library_metadata.tsv", encoding="utf-8"), delimiter="\t")}
DATE = {"Bloemfontein": "2016-04-26", "Christiana": "2017-05-10", "Pretoria": "2017-05-18"}
# test 1: blastn hit, via exact sequence to the placement input names
yin = rd(O / "yeast_placement" / "yeast_asvs.fasta")
name_of = {v.upper(): k for k, v in yin.items()}
hits = {l.split("\t")[0] for l in open(O / "yeast_placement" / "yeast_vs_tremello.tsv") if l.strip()}
# test 2: rerun tree
tree = {r["asv"]: r["placed_within_Tremellomycetes"] == "yes"
        for r in csv.DictReader(open(O / "yeast_placement_71" / "fig3_71_placement.tsv"), delimiter="\t")}

recs = []
for h, s in md5.items():
    assert h in tab, h
    counts = tab[h]
    occ = [libs[i] for i, c in enumerate(counts) if c > 0]
    sites = sorted({meta[l]["site"] for l in occ})
    dates = sorted(DATE[x] for x in sites)
    t1 = name_of[s] in hits
    t2 = tree[h[:8]]
    recs.append({"hash": h, "seq": s, "reads": int(sum(counts)), "n_lib": len(occ), "sites": sites,
                 "date": dates[0] if dates[0] == dates[-1] else f"{dates[0]}/{dates[-1]}", "blast": t1, "tree": t2})
recs.sort(key=lambda r: (-r["reads"], r["hash"]))
NOTE = ("amplicon sequence variant (DADA2) of the ITS2 region, primers ITS3 and ITS4, from merged Illumina MiSeq read pairs "
        "from two runs (2x250 bp and 2x301 bp)")
out, cmap = [], []
for i, r in enumerate(recs, 1):
    clone = f"KMD_ITS2_ASV{i:03d}"
    both = r["blast"] and r["tree"]
    org = "uncultured Tremellomycetes" if both else "uncultured fungus"
    extra = "" if both else ("; one of 71 ASVs of an unresolved lineage provisionally assigned to Tremellomycetes; "
                             + ("inside the Tremellomycetes reference clade of the placement tree but without a blastn hit "
                                "to the Tremellomycetes reference set" if r["tree"] else
                                "with a blastn hit to the Tremellomycetes reference set but outside the reference clade of "
                                "the placement tree"))
    r.update(clone=clone, organism=org)
    out.append({"Sequence_ID": clone, "clone": clone, "organism": org,
                "isolation_source": "surface-sterilised internal tissue of Searsia lancea", "host": "Searsia lancea",
                "geo_loc_name": "South Africa", "collection_date": r["date"], "environmental_sample": "TRUE",
                "BioProject": "", "note": NOTE + extra})
    cmap.append({"clone": clone, "ASV_ID_md5": r["hash"], "total_reads": r["reads"], "libraries": r["n_lib"],
                 "sites": ", ".join(r["sites"]), "blastn_hit_Tremellomycetes_refs": "yes" if r["blast"] else "no",
                 "inside_Tremellomycetes_clade_rerun_tree": "yes" if r["tree"] else "no",
                 "organism": org, "taxid": TAX[org]})

# retire the Prompt 3 GenBank files, write the new ones
old = ROOT / "_retired" / "NCBI_GenBank_pre-prompt5"
old.mkdir(exist_ok=True)
for f in ("Tremellomycetes_68_ASVs.fasta", "source_table.tsv", "excluded_3_ASVs.tsv"):
    if (GB / f).exists(): shutil.move(str(GB / f), str(old / f))
with open(GB / "KMD_ITS2_71_ASVs.fasta", "w", newline="\n") as fh:
    for r in recs: fh.write(f">{r['clone']}\n{r['seq']}\n")
for name, data in (("source_modifiers.tsv", out), ("CLONE_MAP.tsv", cmap)):
    with open(GB / name, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0].keys()), delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(data)

# Supplementary Table S1: add the clone name to the ASV catalogue (tsv and xlsx)
cat = O / "asv_catalogue_final" / "ASV_catalogue.tsv"
bak = ROOT / "_retired" / "ASV_catalogue_pre-prompt5_backup.tsv"
if not bak.exists(): shutil.copy2(cat, bak)
crow = list(csv.DictReader(open(bak, encoding="utf-8"), delimiter="\t"))
by = {r["hash"]: r["clone"] for r in recs}
cf = list(crow[0].keys())
cf.insert(1, "genbank_clone")
for r in crow: r["genbank_clone"] = by.get(r["ASV_ID_md5"], "")
assert sum(1 for r in crow if r["genbank_clone"]) == 71
with open(cat, "w", encoding="utf-8", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=cf, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(crow)
try:
    import openpyxl
    wb = openpyxl.Workbook(); ws = wb.active; ws.title = "ASV_catalogue"; ws.append(cf)
    for r in crow: ws.append([int(r[c]) if c in ("length_bp", "total_reads", "prevalence_n_samples") and r[c].isdigit() else r[c] for c in cf])
    xb = ROOT / "_retired" / "ASV_catalogue_pre-prompt5_backup.xlsx"
    if not xb.exists(): shutil.copy2(cat.with_suffix(".xlsx"), xb)
    wb.save(cat.with_suffix(".xlsx")); xl = "xlsx rewritten"
except ImportError:
    xl = "openpyxl missing: xlsx NOT updated"
n_both = sum(1 for r in recs if r["blast"] and r["tree"])
print(f"71 ASVs: both tests {n_both}; tree only {[r['hash'][:8] for r in recs if r['tree'] and not r['blast']]}; "
      f"BLAST only {[r['hash'][:8] for r in recs if r['blast'] and not r['tree']]}; neither "
      f"{[r['hash'][:8] for r in recs if not r['blast'] and not r['tree']]}")
print("uncultured fungus clones:", [(r["clone"], r["hash"][:8], r["reads"]) for r in recs if r["organism"] == "uncultured fungus"])
print("date values:", sorted({r['date'] for r in recs}), "|", xl)
