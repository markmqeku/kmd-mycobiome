#!/usr/bin/env python
"""H-1 depth for the core 16; H-2 tighten the yeast evidence."""
import csv, urllib.request, urllib.parse, zipfile, re
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"

meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
per = {s: {} for s in cols}
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
depth = {s: sum(per[s].values()) for s in cols}
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]

print("=== H-1: DEPTH for the core 16 (Christiana 9 + Pretoria 7) ===")
print(f"core samples found: {len(core)}")
srt = sorted(((depth[s], s) for s in core))
print("per-sample fungal depth (ascending):")
for d, s in srt:
    print(f"   {s:6} {meta[s]['location']:11} {meta[s]['condition_std']:13} {meta[s]['tissue']:14} {d:>9,}")
lowest, lowest_s = srt[0]
print(f"\nlowest core sample: {lowest_s} = {lowest:,} reads")
print("\nretention in steps beyond 10,000:")
step = 10000
while True:
    kept = sum(1 for s in core if depth[s] >= step)
    print(f"   depth {step:>7,}: {kept}/16" + ("   <-- first loss" if kept < 16 else ""))
    if kept < 16: break
    step += 5000
    if step > 300000: break
print(f"\nHIGHEST DEPTH RETAINING ALL 16 = {lowest:,} (the minimum core sample, {lowest_s})")
print(f"  practical primary depth (rounded down to 1,000): {1000*(lowest//1000):,}")

# ---------- H-2 ----------
print("\n\n=== H-2a: alignment extent behind the '100% identity' Filobasidium hits ===")
D = f"{O}/yeast_placement"
rows = [r for r in csv.reader(open(f"{D}/yeast_vs_tremello.tsv"), delimiter="\t") if len(r) >= 6]
best = {}
for r in rows:
    q, s, pid, alen, qcov = r[0], r[1], float(r[2]), int(r[3]), float(r[4])
    if q not in best or pid > best[q][1]:
        best[q] = (s, pid, alen, qcov)
def reads(q):
    m = re.search(r"reads(\d+)", q); return int(m.group(1)) if m else 0
top = sorted((q for q in best), key=reads, reverse=True)[:10]
print(f"{'ASV':34} {'reference':30} {'%id':>6} {'aln_len':>8} {'qcov':>6} {'ASV_len':>8}")
seqs = {}
sid = None; buf = []
for line in open(f"{D}/yeast_asvs.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:], []
    else: buf.append(line)
if sid: seqs[sid] = "".join(buf)
for q in top:
    s, pid, alen, qcov = best[q]
    print(f"{q:34} {s[:30]:30} {pid:>6.1f} {alen:>8} {qcov:>5.0f}% {len(seqs.get(q,'')):>8}")

print("\n=== H-2b: the ASVs with NO Tremellomycetes hit ===")
allq = set(seqs)
nohit = sorted(allq - set(best), key=reads, reverse=True)
print(f"count: {len(nohit)}")
ncbi = {}
try:
    for r in csv.reader(open(f"{O}/nontarget_screen/blast_priority.tsv"), delimiter="\t"):
        if len(r) >= 8 and (r[0] not in ncbi or float(r[2]) > ncbi[r[0]][0]):
            ncbi[r[0]] = (float(r[2]), float(r[4]), r[7])
except FileNotFoundError: pass
for q in nohit:
    md5 = q.split("_")[1]
    hit = next((v for k, v in ncbi.items() if k.startswith(md5)), None)
    print(f"   {q:34} len={len(seqs[q]):>4} reads={reads(q):>7,}  "
          f"NCBI: {hit[2][:46]+' id='+str(hit[0]) if hit else 'no confident hit'}")

print("\n=== H-2c: do the top references postdate UNITE 10.0? ===")
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
accs = sorted({best[q][0].split("_")[-1] for q in top})
print(f"accessions: {accs}")
try:
    ids = ",".join(accs)
    x = urllib.request.urlopen(f"{E}/esummary.fcgi?db=nuccore&id={urllib.parse.quote(ids)}&retmode=json", timeout=60).read().decode()
    import json
    d = json.loads(x)
    for k, v in d.get("result", {}).items():
        if k == "uids": continue
        print(f"   {v.get('accessionversion','?'):16} createdate={v.get('createdate','?')} "
              f"updatedate={v.get('updatedate','?')}  {v.get('title','')[:52]}")
except Exception as e:
    print(f"   esummary failed: {type(e).__name__}: {e}")
print("   UNITE 10.0 release date: 2024-04-04 (version 10.0, released 4 April 2024)")
