#!/usr/bin/env python
"""H-2a corrected: select nearest reference by BITSCORE (alignment quality), not by raw %identity.
Selecting by %identity alone lets a 31 bp perfect match outrank a 350 bp 97% match - which is what
distorted the earlier G-3 summary. Rerun and report identity WITH alignment length and coverage."""
import csv, re, subprocess, os
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
D = f"{O}/yeast_placement"

# re-blast with bitscore and more targets so the ranking is meaningful
out = f"{D}/yeast_vs_tremello_full.tsv"
if not os.path.exists(out):
    subprocess.run(["blastn", "-db", f"{D}/tremellodb", "-query", f"{D}/yeast_asvs.fasta",
                    "-outfmt", "6 qseqid sseqid pident length qcovhsp evalue bitscore",
                    "-max_target_seqs", "20", "-evalue", "1e-5", "-num_threads", "4", "-out", out], check=True)

rows = [r for r in csv.reader(open(out), delimiter="\t") if len(r) >= 7]
best = {}
for r in rows:
    q, s, pid, alen, qcov, ev, bits = r[0], r[1], float(r[2]), int(r[3]), float(r[4]), r[5], float(r[6])
    if q not in best or bits > best[q][5]:
        best[q] = (s, pid, alen, qcov, ev, bits)

def reads(q):
    m = re.search(r"reads(\d+)", q); return int(m.group(1)) if m else 0
seqs, sid, buf = {}, None, []
for line in open(f"{D}/yeast_asvs.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:], []
    else: buf.append(line)
if sid: seqs[sid] = "".join(buf)

print("=== H-2a CORRECTED: nearest Tremellomycetes reference selected by BITSCORE ===")
print(f"{'ASV':30} {'reference':32} {'%id':>6} {'aln':>5} {'qcovhsp':>8} {'ASVlen':>7} {'bits':>7}")
top = sorted(best, key=reads, reverse=True)[:10]
for q in top:
    s, pid, alen, qcov, ev, bits = best[q]
    print(f"{q[:30]:30} {s[:32]:32} {pid:>6.1f} {alen:>5} {qcov:>7.0f}% {len(seqs.get(q,'')):>7} {bits:>7.0f}")

gen = Counter(); rd = Counter()
for q, v in best.items():
    g = v[0].split("_")[1] if v[0].startswith("REF_") else "?"
    gen[g] += 1; rd[g] += reads(q)
print(f"\nnearest reference genus (by bitscore), all {len(best)} ASVs with a hit:")
for g, n in gen.most_common():
    print(f"   {g:20} {n:>4} ASVs   {rd[g]:>10,} reads")

al = [v[2] for v in best.values()]; pi = [v[1] for v in best.values()]
print(f"\nalignment length to nearest reference: min={min(al)} median={sorted(al)[len(al)//2]} max={max(al)}")
print(f"identity to nearest reference        : min={min(pi):.1f} median={sorted(pi)[len(pi)//2]:.1f} max={max(pi):.1f}")
short = sum(1 for v in best.values() if v[2] < 100)
print(f"ASVs whose best alignment is <100 bp : {short} of {len(best)}  "
      f"(these carry essentially no genus-level information)")

# compare against the NCBI nt evidence, which had full-length alignments
print("\n=== for contrast, the NCBI nt hits (from the earlier remote BLAST) ===")
try:
    n = {}
    for r in csv.reader(open(f"{O}/nontarget_screen/blast_priority.tsv"), delimiter="\t"):
        if len(r) >= 8:
            q = r[0]; pid = float(r[2]); alen = int(r[3]); qcov = float(r[4])
            if q not in n or pid > n[q][0]: n[q] = (pid, alen, qcov, r[7])
    shown = 0
    for q, v in sorted(n.items(), key=lambda x: -x[1][1]):
        if "Cryptococcus" in v[3] and shown < 6:
            print(f"   {q[:9]} id={v[0]:.1f} aln={v[1]} qcov={v[2]:.0f}%  {v[3][:48]}")
            shown += 1
except FileNotFoundError:
    print("   (no remote BLAST file)")
