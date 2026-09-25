#!/usr/bin/env python
"""F-4 decisive test: is host sequence PRESENT IN THE RAW READS of Christiana/Pretoria?
If raw reads contain host but the ASV table does not -> pipeline/filtering effect.
If raw reads genuinely lack host -> upstream (wet-lab or sample) effect.
Observation only; no cause asserted."""
import gzip, glob, os
B = "<KMD_ROOT>/__reanalysis_2026-06"
RAW = None
for d in sorted(glob.glob("<KMD_ROOT>/__restructured_*"), reverse=True):
    if os.path.isdir(os.path.join(d, "raw_reads")):
        RAW = os.path.join(d, "raw_reads"); break

# most abundant host ASV sequence
host = {l.strip() for l in open(f"{B}/phase3/outputs_27-07-2026/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
seqs, sid, buf = {}, None, []
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line.strip())
if sid: seqs[sid] = "".join(buf)

# rank host ASVs by abundance
lines = [l.rstrip("\n") for l in open(f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); ns = len([h for h in hdr[1:] if h.strip()])
tot = {}
for l in lines[1:]:
    p = l.split("\t"); tot[p[0]] = sum(int(float(x)) for x in p[1:1+ns])
top_host = sorted(host, key=lambda a: -tot.get(a, 0))[:3]

# use internal 30-mers as probes (avoid primer region, tolerate ends)
probes = []
for a in top_host:
    s = seqs[a]
    probes.append((a[:8], s[60:90], s[150:180]))
print("host probes (30-mers from the 3 most abundant host ASVs):")
for t in probes: print(f"   {t[0]}  p1={t[1]}  p2={t[2]}")

SAMPLES = [("Bloemfontein", "P16"), ("Bloemfontein", "P26"), ("Bloemfontein", "P5"),
           ("Christiana", "25-C"), ("Christiana", "OL-C"), ("Christiana", "27-C"),
           ("Pretoria", "28-C"), ("Pretoria", "OL-P"), ("Pretoria", "29-C")]

def openf(p): return gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, errors="replace")
def find(sid, rn):
    for ext in (".fastq", ".fastq.gz"):
        g = glob.glob(f"{RAW}/{sid}_S*_L001_{rn}_001{ext}")
        if g: return g[0]
    return None

N = 200000  # reads scanned per sample
print(f"\nscanning first {N:,} R1 reads per sample for host probes\n")
print(f"{'site':14} {'sample':7} {'reads':>9} {'probe hits':>11} {'% of scanned':>13}")
for site, sid in SAMPLES:
    p = find(sid, "R1")
    if not p: print(f"{site:14} {sid:7} FILE NOT FOUND"); continue
    hits = 0; n = 0
    with openf(p) as fh:
        for i, line in enumerate(fh):
            if i >= N * 4: break
            if i % 4 == 1:
                n += 1
                s = line.strip()
                if any(p1 in s or p2 in s for _, p1, p2 in probes): hits += 1
    print(f"{site:14} {sid:7} {n:>9,} {hits:>11,} {100*hits/n if n else 0:>12.2f}%")
