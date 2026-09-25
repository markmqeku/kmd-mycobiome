#!/usr/bin/env python
"""Item 72 continued. Two corrections and one better source.

1. Magurran 2004 was 'verified' to a Choice Reviews Online record, which is a BOOK REVIEW, not
   the book. Same failure class as the Faculty Opinions and PeerJ peer-review records. Demoted,
   and 'choice reviews' added to the record-type guard.
2. QIIME2 artifacts carry a citations.bib written by the plugins themselves. That is a
   provenance-grade source for the software references, better than a bibliographic search.
   Harvest it rather than guessing.
"""
import csv, os, re, zipfile
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
C = f"{O}/clean"
ARTIFACTS = ["unite_current/unite10_seqs.qza", "unite_current/unite10_tax.qza",
             "taxonomy_unite10/search.qza", "table_exp/../table.qza"]
ARTIFACTS = [a for a in ARTIFACTS if os.path.exists(os.path.join(O, a))]
for root, _, files in os.walk(O):
    for f in files:
        if f.endswith(".qza"):
            p = os.path.relpath(os.path.join(root, f), O)
            if p not in ARTIFACTS: ARTIFACTS.append(p)

entries = {}
for a in ARTIFACTS:
    try: z = zipfile.ZipFile(os.path.join(O, a))
    except Exception: continue
    for n in z.namelist():
        if not n.endswith(".bib"): continue
        txt = z.read(n).decode("utf-8", "replace")
        for blk in txt.split("@")[1:]:
            key = blk.split("{", 1)[1].split(",", 1)[0].strip() if "{" in blk else ""
            def fld(name):
                m = re.search(name + r"\s*=\s*\{(.+?)\}\s*,?\s*\n", blk, re.S)
                return " ".join(m.group(1).split()) if m else ""
            e = {"title": fld("title"), "doi": fld("doi"), "journal": fld("journal"),
                 "year": fld("year"), "volume": fld("volume"), "pages": fld("pages"),
                 "author": fld("author"), "booktitle": fld("booktitle"),
                 "publisher": fld("publisher")}
            if e["title"] and e["doi"]: entries.setdefault(e["doi"], e)
print(f"artifacts scanned: {len(ARTIFACTS)};  distinct DOI-bearing citations found: {len(entries)}\n")

WANT = {
 "Bolyen et al. 2019":   ("qiime 2", "2019"),
 "Callahan et al. 2016": ("dada2", "2016"),
 "Rognes et al. 2016":   ("vsearch", "2016"),
 "Martin 2011":          ("cutadapt", "2011"),
 "Nilsson et al. 2019":  ("unite", "2019"),
}
matched = {}
for key, (tok, yr) in WANT.items():
    for doi, e in entries.items():
        if tok in e["title"].lower() and (e["year"] == yr or not e["year"]):
            matched[key] = (doi, e); break
print("=== software references recoverable from artifact provenance ===")
for key in WANT:
    if key in matched:
        doi, e = matched[key]
        j = e["journal"] or e["booktitle"] or e["publisher"]
        print(f"   {key:24} doi:{doi}")
        print(f"      {e['title'][:92]}")
        print(f"      {j}  {e['year']}  vol {e['volume']}  pp {e['pages']}")
    else:
        print(f"   {key:24} not present in any artifact bib")

# ---- apply: demote Magurran, promote anything provenance confirms ----
p = f"{C}/reference_list.tsv"
rows = list(csv.reader(open(p), delimiter="\t")); hdr, body = rows[0], rows[1:]
changed = []
for r in body:
    if r[0] == "Magurran 2004" and r[1] == "VERIFIED":
        r[1] = "UNVERIFIED"; r[6] = ""
        r[7] = ("Crossref returned a Choice Reviews Online record, which is a BOOK REVIEW, not the "
                "book. Magurran (2004) Measuring Biological Diversity is a Blackwell monograph; "
                "confirm edition, publisher and ISBN by hand. No DOI carried.")
        changed.append(("demoted", r[0]))
    if r[0] in matched:
        doi, e = matched[r[0]]
        if r[1] != "VERIFIED" or not r[6]:
            r[1] = "VERIFIED"; r[2] = e["title"]; r[3] = e["year"]
            r[4] = (e["journal"] or e["booktitle"] or e["publisher"])
            r[6] = doi
            r[7] = ("Verified from the QIIME2 artifact citation bundle written by the plugin itself "
                    f"(vol {e['volume']}, pp {e['pages']})")
            changed.append(("verified from provenance", r[0]))
with open(p, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)
print("\n=== changes applied ===")
for what, k in changed: print(f"   {what:26} {k}")
v = sum(1 for r in body if r[1] == "VERIFIED")
print(f"\nreference list: {len(body)} entries, {v} verified, {len(body)-v} not verified")
