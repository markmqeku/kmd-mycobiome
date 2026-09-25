#!/usr/bin/env python
"""13a. Define the NCBI screen set for the WHOLE 1,173-ASV fungal table on evidence,
not on a convenience cap.

For every ASV, take the best UNITE 10.0 alignment (identity AND alignment length).
An ASV is 'securely fungal by alignment' only if it aligns to a fungal reference over
>= MINLEN bp at >= MINID identity: a plant, algal or protist ITS2 cannot do that.
Every ASV that fails that test goes to NCBI nt, whatever its assigned rank.
The pass/fail counts are printed so the exclusion is auditable rather than silent."""
import csv, zipfile, os
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; S = f"{O}/full_screen"
MINLEN, MINID = 200, 90.0

best = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 4: continue
    try: pid, alen = float(p[2]), int(p[3])
    except ValueError: continue
    sc = (pid, alen)
    if p[0] not in best or alen*pid > best[p[0]][1]*best[p[0]][0]:
        best[p[0]] = sc

idx = list(csv.DictReader(open(f"{S}/asv_index.tsv"), delimiter="\t"))
seqs, sid, buf = {}, None, []
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line)
if sid: seqs[sid] = "".join(buf)

need, secure = [], []
for r in idx:
    a = r["asv_id"]
    pid, alen = best.get(a, (0.0, 0))
    (secure if (alen >= MINLEN and pid >= MINID) else need).append((a, pid, alen, int(r["reads"]), r["tier"]))
need.sort(key=lambda x: -x[3])

print(f"criterion for 'securely fungal by alignment': UNITE hit >= {MINLEN} bp at >= {MINID}% identity")
print(f"  securely fungal, NOT sent to NCBI: {len(secure):>5} ASVs, {sum(x[3] for x in secure):>9,} reads")
print(f"  sent to NCBI nt:                   {len(need):>5} ASVs, {sum(x[3] for x in need):>9,} reads")
from collections import Counter
c = Counter(x[4] for x in need)
print("  screen set by resolution tier: " + ", ".join(f"tier{t}={c[t]}" for t in sorted(c)))
print(f"\n  weakest 'secure' ASV: {min(secure, key=lambda x: x[1]*x[2])[1:3]} (identity, aln_len)")

with open(f"{S}/screen_set.fasta", "w") as fh:
    for a, pid, alen, rd, t in need: fh.write(f">{a}\n{seqs[a]}\n")
with open(f"{S}/screen_basis.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["asv_id","sent_to_ncbi","unite_best_identity","unite_best_aln_len","reads","tier"])
    for a, pid, alen, rd, t in need:   w.writerow([a,"yes",f"{pid:.1f}",alen,rd,t])
    for a, pid, alen, rd, t in secure: w.writerow([a,"no",f"{pid:.1f}",alen,rd,t])
print(f"\nscreen set -> {S}/screen_set.fasta ({len(need)} sequences)")
print(f"basis      -> {S}/screen_basis.tsv (all 1,173 ASVs, with the reason for each)")
