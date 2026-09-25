#!/usr/bin/env python
"""Item 8: verify current accepted status for genera appearing in the FINAL composition
(95% genus threshold). Scoped: names from the superseded UNITE 8.2 comparison are excluded."""
import csv, zipfile, re, json, urllib.request, urllib.parse, time
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

def genus_of(t, p, thr=95.0):
    if not t or t.lower().startswith("unassigned") or p < thr: return ""
    for part in t.split(";"):
        part = part.strip()
        if part.startswith("g__"):
            v = part[3:].strip()
            if v.lower() in PH or v.lower().endswith("_incertae_sedis"): return ""
            return v
    return ""

tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v

lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
reads = Counter()
for l in lines[1:]:
    q = l.split("\t")
    reads[q[0]] = sum(int(float(x)) for x in q[1:1+len(cols)])
grand = sum(reads.values())
g_reads = Counter()
for a in tax:
    g = genus_of(tax[a], pid.get(a, 0.0))
    if g: g_reads[g] += reads.get(a, 0)
genera = [g for g, v in g_reads.items() if v > 0]
print(f"Genera in the FINAL composition at the 95% threshold: {len(genera)}")
print("(scope: genera with >0 reads in the fungi-only table; UNITE 8.2-only names excluded)\n")

def check(name):
    out = {}
    try:
        u = ("http://www.indexfungorum.org/ixfwebservice/fungus.asmx/NameSearch?SearchText="
             + urllib.parse.quote(name) + "&AnywhereInText=false&MaxNumber=6")
        x = urllib.request.urlopen(urllib.request.Request(u, headers={"User-Agent":"KMD/1.0"}), timeout=45).read().decode("utf-8","replace")
        recs = x.split("<IndexFungorum>")
        for rec in recs[1:]:
            nm = re.search(r"<NAME_x0020_OF_x0020_FUNGUS>(.*?)</", rec)
            au = re.search(r"<AUTHORS>(.*?)</", rec)
            cu = re.search(r"<CURRENT_x0020_NAME>(.*?)</", rec)
            rk = re.search(r"<RANK>(.*?)</", rec)
            if nm and nm.group(1).strip().lower() == name.lower():
                out = {"name": nm.group(1), "authors": au.group(1) if au else "",
                       "current": cu.group(1) if cu else "", "rank": rk.group(1) if rk else ""}
                break
    except Exception as e:
        out = {"error": f"{type(e).__name__}"}
    return out

print(f"{'genus':22} {'reads%':>7}  {'IF authors':34} {'IF current name':26} status")
for g in sorted(genera, key=lambda x: -g_reads[x]):
    r = check(g); time.sleep(0.3)
    if "error" in r:
        print(f"{g:22} {100*g_reads[g]/grand:6.2f}%  LOOKUP FAILED ({r['error']})")
        continue
    if not r:
        print(f"{g:22} {100*g_reads[g]/grand:6.2f}%  {'NOT FOUND in Index Fungorum':34} {'':26} FLAG")
        continue
    cur = r.get("current","").strip()
    status = "accepted (current = itself)" if (not cur or cur.lower() == g.lower()) else f"SYNONYM -> {cur}"
    print(f"{g:22} {100*g_reads[g]/grand:6.2f}%  {r.get('authors','')[:34]:34} {cur[:26]:26} {status}")
