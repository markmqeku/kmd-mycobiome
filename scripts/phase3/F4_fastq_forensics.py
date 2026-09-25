#!/usr/bin/env python
"""F-4 / F-2 / F-3: raw FASTQ forensics per site. Observations only; no cause asserted."""
import gzip, os, glob, statistics as st
from collections import Counter

RAW = None
for d in sorted(glob.glob("<KMD_ROOT>/__restructured_*"), reverse=True):
    if os.path.isdir(os.path.join(d, "raw_reads")):
        RAW = os.path.join(d, "raw_reads"); break
print(f"resolved raw_reads: {RAW}\n")

ITS3 = "GCATCGATGAAGAACGCAGC"
ITS4 = "TCCTCCGCTTATTGATATGC"
# common Illumina adapter fragments
ADAPTERS = {"TruSeq_R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            "TruSeq_R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
            "Nextera":   "CTGTCTCTTATACACATCT"}

SAMPLES = [("Bloemfontein", "P16"), ("Bloemfontein", "P26"), ("Bloemfontein", "P5"),
           ("Christiana", "25-C"), ("Christiana", "OL-C"),
           ("Pretoria", "28-C"), ("Pretoria", "OL-P")]

def openf(p):
    return gzip.open(p, "rt", errors="replace") if p.endswith(".gz") else open(p, errors="replace")

def find_file(sid, rn):
    for ext in (".fastq", ".fastq.gz"):
        g = glob.glob(f"{RAW}/{sid}_S*_L001_{rn}_001{ext}")
        if g: return g[0]
    return None

N = 4000
for site, sid in SAMPLES:
    p = find_file(sid, "R1"); p2 = find_file(sid, "R2")
    if not p:
        print(f"{site} {sid}: R1 NOT FOUND"); continue
    heads, lens, quals, prim_at0, prim_any, adapters = [], [], [], 0, 0, Counter()
    with openf(p) as fh:
        for i, line in enumerate(fh):
            if i >= N * 4: break
            m = i % 4
            if m == 0: heads.append(line.strip())
            elif m == 1:
                s = line.strip(); lens.append(len(s))
                if s.startswith(ITS3): prim_at0 += 1
                if ITS3 in s[:60]: prim_any += 1
                for an, aseq in ADAPTERS.items():
                    if aseq[:15] in s: adapters[an] += 1
            elif m == 3:
                q = line.strip()
                if q: quals.append(st.mean(ord(c) - 33 for c in q))
    n = len(lens)
    lc = Counter(lens)
    hdr0 = heads[0] if heads else ""
    inst = hdr0.split(":")[0:4] if ":" in hdr0 else [hdr0]
    print(f"--- {site:13} {sid:6} ({os.path.basename(p)}) ---")
    print(f"    header       : {hdr0[:72]}")
    print(f"    instrument/run/flowcell/lane : {':'.join(inst)}")
    print(f"    reads sampled: {n}   distinct lengths: {len(lc)}")
    print(f"    top lengths  : {lc.most_common(4)}")
    print(f"    max len      : {max(lens)}   uniform-at-max: {100*lc[max(lens)]/n:.1f}%")
    print(f"    ITS3 at pos0 : {100*prim_at0/n:5.1f}%     ITS3 within first 60 bp: {100*prim_any/n:5.1f}%")
    print(f"    adapter frags: {dict(adapters) if adapters else 'none detected'}")
    print(f"    mean read Q  : {st.mean(quals):.1f}")
    print()

# F-2: ITS4 on R2, per site
print("=== F-2: ITS4 at start of R2, per site ===")
for site, sid in SAMPLES:
    p2 = find_file(sid, "R2")
    if not p2: continue
    tot = hit = 0
    with openf(p2) as fh:
        for i, line in enumerate(fh):
            if i >= N * 4: break
            if i % 4 == 1:
                tot += 1
                if line.strip().startswith(ITS4): hit += 1
    print(f"   {site:13} {sid:6}: {100*hit/tot:5.1f}%")

# F-3: length distribution, host vs fungal ASVs
print("\n=== F-3: ASV sequence length, host vs fungal ===")
B = "<KMD_ROOT>/__reanalysis_2026-06"
host = {l.strip() for l in open(f"{B}/phase3/outputs_27-07-2026/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
seqs, sid_, buf = {}, None, []
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid_: seqs[sid_] = "".join(buf)
        sid_, buf = line[1:].split()[0], []
    else: buf.append(line.strip())
if sid_: seqs[sid_] = "".join(buf)
hl = [len(s) for a, s in seqs.items() if a in host]
fl = [len(s) for a, s in seqs.items() if a not in host]
for lab, v in (("HOST (n=%d)" % len(hl), hl), ("FUNGAL/other (n=%d)" % len(fl), fl)):
    print(f"   {lab:22} min={min(v)} q25={int(st.quantiles(v,n=4)[0])} median={int(st.median(v))} "
          f"q75={int(st.quantiles(v,n=4)[2])} max={max(v)}")
