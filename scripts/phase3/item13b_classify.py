#!/usr/bin/env python
"""13a analysis. Classify every screened ASV by the kingdom of its best-bitscore NCBI hit.
Hits are ranked by BITSCORE, never by identity alone: a short high-identity hit must not
outrank a long one. Lineages come from the NCBI taxonomy service, not from guessing at
the description string."""
import csv, json, os, time, urllib.request, urllib.parse
from collections import defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; S = f"{O}/full_screen"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
MINID, MINLEN = 90.0, 150            # evidence floor for calling an ASV non-fungal

hits = defaultdict(list)
srcs = [p for p in (f"{S}/screen_hits.tsv", f"{S}/screen_hits_b.tsv") if os.path.exists(p)]
for src in srcs:
    for r in csv.reader(open(src), delimiter="\t"):
        if len(r) < 11: continue
        q, sacc, pid, alen, qcov, ev, bits, taxids, names, sk, title = r[:11]
        try: hits[q].append((float(bits), float(pid), int(alen), float(qcov), sacc,
                             taxids.split(";")[0], names, title))
        except ValueError: continue
print("hit files read: " + ", ".join(os.path.basename(p) for p in srcs))

# Coverage check: an ASV with no hit line is NOT evidence of anything. Report the gap
# rather than letting a silently unscreened sequence look like a clean result.
sent = {l[1:].split()[0] for l in open(f"{S}/screen_set.fasta") if l.startswith(">")}
nohit = sent - set(hits)
print(f"screen set: {len(sent)} ASVs | returned hits: {len(set(hits) & sent)} | no hit: {len(nohit)}")
for q in hits: hits[q].sort(reverse=True)
print(f"ASVs with at least one NCBI hit: {len(hits)}")

# ---- lineages for every taxid seen in a top hit ----
taxids = sorted({h[5] for q in hits for h in hits[q] if h[5].isdigit()})
lin = {}
cache = f"{S}/taxid_lineage.json"
if os.path.exists(cache): lin = json.load(open(cache))
todo = [t for t in taxids if t not in lin]
print(f"taxids needing lineage: {len(todo)} of {len(taxids)}")
for i in range(0, len(todo), 150):
    chunk = todo[i:i+150]
    for attempt in (1, 2, 3):
        try:
            x = urllib.request.urlopen(
                f"{E}/efetch.fcgi?db=taxonomy&id={urllib.parse.quote(','.join(chunk))}&retmode=xml",
                timeout=120).read().decode("utf-8", "replace")
            break
        except Exception as e:
            print(f"   efetch attempt {attempt} failed: {type(e).__name__}"); time.sleep(10); x = ""
    # Parse as XML, not by regex: <TaxId> also occurs nested inside <LineageEx>, and a
    # greedy scan silently keys the wrong record (it returns 131567 'cellular organisms'
    # for every taxon after the first).
    import xml.etree.ElementTree as ET
    try:
        root = ET.fromstring(x)
    except ET.ParseError as e:
        print(f"   XML parse failed for this chunk: {e}"); continue
    for t in root.findall("Taxon"):                    # top-level records only
        tid = t.findtext("TaxId")
        if tid:
            lin[tid] = {"name": t.findtext("ScientificName") or "",
                        "lineage": t.findtext("Lineage") or ""}
    time.sleep(1)
json.dump(lin, open(cache, "w"))
print(f"lineages resolved: {len(lin)}")

NONFUNGAL = ("Viridiplantae","Metazoa","Bacteria","Archaea","Amoebozoa",
             "Rhodophyta","Stramenopiles","Alveolata","Rhizaria","Viruses")
def kingdom(taxid):
    d = lin.get(taxid)
    if not d: return "unknown"
    L = d["lineage"] + "; " + d["name"]
    if "Fungi" in L: return "Fungi"
    for k in NONFUNGAL:
        if k in L: return k
    # 'uncultured eukaryote', 'environmental samples' and similar are NOT evidence of a
    # non-fungal sequence: most environmental ITS2 clones are fungal. Never remove on these.
    if "Eukaryota" in L: return "unplaced eukaryote"
    return "unclassified"

idxrows = {r["asv_id"]: r for r in csv.DictReader(open(f"{S}/asv_index.tsv"), delimiter="\t")}
out, nonfungal = [], []
for q, hl in hits.items():
    bits, pid, alen, qcov, sacc, tid, names, title = hl[0]
    kd = kingdom(tid)
    ir = idxrows.get(q, {})
    strong = (pid >= MINID and alen >= MINLEN)
    # Adversarial guard: a single mislabelled record must not condemn an ASV. If any fungal
    # hit scores within 5% of the best non-fungal hit, the evidence is contested, not decisive.
    contest = [h for h in hl if kingdom(h[5]) == "Fungi" and h[0] >= 0.95*bits]
    if kd in NONFUNGAL and strong and not contest: call = "NON-FUNGAL"
    elif kd in NONFUNGAL and strong and contest:   call = "contested"
    elif kd == "Fungi":                            call = "fungal"
    else:                                          call = "inconclusive"
    row = [q, call, kd, f"{pid:.1f}", alen, f"{qcov:.0f}", f"{bits:.0f}", sacc,
           lin.get(tid, {}).get("name", names), title[:70],
           ir.get("reads",""), ir.get("n_samples",""), ir.get("samples","")]
    out.append(row)
    if call == "NON-FUNGAL": nonfungal.append(row)

for q in sorted(nohit):
    ir = idxrows.get(q, {})
    out.append([q, "no NCBI hit", "n/a", "", "", "", "", "", "", "",
                ir.get("reads",""), ir.get("n_samples",""), ir.get("samples","")])

hdr = ["asv_id","call","kingdom","identity","aln_len","qcov","bitscore","accession",
       "organism","hit_title","reads","n_samples","samples"]
with open(f"{S}/screen_calls.tsv","w",newline="",encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr)
    w.writerows(sorted(out, key=lambda r: (r[1] != "NON-FUNGAL", -int(r[10] or 0))))
with open(f"{S}/nonfungal_asvs.txt","w") as fh:
    for r in sorted(nonfungal, key=lambda r: -int(r[10] or 0)): fh.write(r[0] + "\n")

from collections import Counter
print("\ncall summary:", dict(Counter(r[1] for r in out)))
print("kingdom of best hit:", dict(Counter(r[2] for r in out)))
print(f"\n=== NON-FUNGAL ASVs ({len(nonfungal)}) ===")
print(f"{'asv':14} {'kingdom':16} {'%id':>6} {'len':>5} {'cov':>4} {'reads':>8} {'n':>3} organism")
for r in sorted(nonfungal, key=lambda r: -int(r[10] or 0)):
    print(f"{r[0][:12]:14} {r[2]:16} {r[3]:>6} {r[4]:>5} {r[5]:>4} {int(r[10]):>8,} {r[11]:>3} {r[8][:34]}")
    print(f"{'':14} samples: {r[12]}")
print(f"\ntotal reads in non-fungal ASVs: {sum(int(r[10] or 0) for r in nonfungal):,}")
print(f"calls -> {S}/screen_calls.tsv ; ids -> {S}/nonfungal_asvs.txt")
