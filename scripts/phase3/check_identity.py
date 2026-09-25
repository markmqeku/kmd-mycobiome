import zipfile, csv, io, statistics as st
B = "<KMD_ROOT>/__reanalysis_2026-06"
z = zipfile.ZipFile(f"{B}/phase3/outputs_27-07-2026/taxonomy_unite10/search.qza")
f = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
txt = z.read(f).decode("utf-8", "replace")

# BLAST6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
best = {}
for line in txt.splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    q, pid = p[0], p[2]
    try: pid = float(pid)
    except ValueError: continue
    if q not in best or pid > best[q]: best[q] = pid

cat = list(csv.DictReader(open(f"{B}/phase3/outputs_27-07-2026/asv_catalogue/ASV_catalogue.tsv"), delimiter="\t"))
by_id = {r["ASV_ID_md5"]: r for r in cat}

def summarise(label, ids):
    v = [best[i] for i in ids if i in best]
    if not v:
        print(f"{label}: no hits recorded"); return
    print(f"{label}: n={len(v)}  best-hit %identity  min={min(v):.1f}  median={st.median(v):.1f}  max={max(v):.1f}"
          f"  |  >=97%: {sum(1 for x in v if x>=97)}  <90%: {sum(1 for x in v if x<90)}")

print("=== BEST-HIT PERCENT IDENTITY (from vsearch search results) ===")
sp  = [r["ASV_ID_md5"] for r in cat if r["species"]]
gen = [r["ASV_ID_md5"] for r in cat if r["genus"] and not r["species"]]
hom = [r["ASV_ID_md5"] for r in cat if "Homophron" in r["display_label"]]
summarise("ASVs given a SPECIES name ", sp)
summarise("ASVs given a GENUS only   ", gen)
summarise("Homophron spadiceum ASVs  ", hom)
print()
print("Homophron ASVs with best-hit identity >=97%:",
      sum(1 for i in hom if best.get(i, 0) >= 97), "of", len(hom))
print("Species-level ASVs with best-hit identity <97%:",
      sum(1 for i in sp if best.get(i, 0) < 97), "of", len(sp))
top = sorted(hom, key=lambda i: -int(by_id[i]["total_reads"]))[:6]
print("\ntop Homophron ASVs (reads, best-hit identity):")
for i in top:
    print(f"   {by_id[i]['display_label']:34} reads={int(by_id[i]['total_reads']):>8,}  best_pident={best.get(i,float('nan')):.1f}")
