#!/usr/bin/env python
"""KMD PROMPT 4, R0c. Size-selection hypothesis. Read-weighted ASV length of the 149 host ASVs against the
71 ASVs of the unresolved Tremellomycetes lineage (and all fungal ASVs), from the denoised sequences and
the full 45-sample table. ASV length is the ITS2 amplicon after primer removal; library fragments on the
Bioanalyzer are longer by the adapters and indices, equally for host and fungal amplicons."""
import re, statistics as st
from collections import Counter
from pathlib import Path

RA = Path(r"<KMD_ROOT>\__reanalysis_2026-06")
O = RA / "phase3" / "outputs_27-07-2026"


def fasta(p):
    d, k = {}, None
    for l in open(p, encoding="utf-8"):
        l = l.strip()
        if l.startswith(">"): k = l[1:].split()[0]; d[k] = ""
        elif k: d[k] += l
    return d


seqs = fasta(RA / "phase2" / "merged2" / "rep_exp" / "dna-sequences.fasta")
host = {l.split()[0] for l in open(O / "nontarget_screen" / "host_asv_ids_ALL.txt", encoding="utf-8") if l.strip()}
lineage = set(re.findall(r"md5=([0-9a-f]{32})", open(O / "clean" / "genbank_submission" / "Tremellomycetes_ASVs.fasta", encoding="utf-8").read()))
fungal = {l[1:].split()[0] for l in open(O / "clean" / "rep_fungi_clean.fasta", encoding="utf-8") if l.startswith(">")}
rows = [l.rstrip("\n").split("\t") for l in open(O / "table_exp" / "feature-table.tsv", encoding="utf-8") if not l.startswith("# ")]
cols = rows[0][1:]
tab = {r[0]: [float(x) for x in r[1:]] for r in rows[1:]}
chr_cols = [i for i, c in enumerate(cols) if c.endswith("-C") and c not in ("28-C", "29-C", "30-C")]
pta_cols = [i for i, c in enumerate(cols) if c.endswith("-P") or c in ("28-C", "29-C")]
assert len(host) == 149 and len(lineage) == 71 and all(a in seqs for a in host | lineage)


def dist(ids, colsel=None, label=""):
    w = Counter()
    for a in ids:
        v = tab.get(a)
        if not v: continue
        n = sum(v) if colsel is None else sum(v[i] for i in colsel)
        if n: w[len(seqs[a])] += n
    tot = sum(w.values())
    if not tot:
        print(f"{label:44} no reads"); return w
    exp = sorted((L for L, n in w.items() for _ in range(1)), key=lambda L: -w[L])
    cum, med = 0, None
    for L in sorted(w):
        cum += w[L]
        if med is None and cum >= tot / 2: med = L
    mode = max(w, key=w.get)
    within = sum(n for L, n in w.items() if abs(L - 372) <= 10) / tot
    lo = min(w); hi = max(w)
    print(f"{label:44} reads {int(tot):>10,}  mode {mode} bp ({w[mode]/tot:5.1%})  median {med}  range {lo} to {hi}  within 372 +/- 10 bp: {within:6.1%}")
    return w


print("ASV length = ITS2 amplicon after primer removal\n")
wh = dist(host, None, "host ASVs (149), all samples")
wl = dist(lineage, None, "lineage ASVs (71), all samples")
dist(lineage, chr_cols, "lineage ASVs, Christiana (9 libraries)")
dist(lineage, pta_cols, "lineage ASVs, Pretoria (7 libraries)")
dist(fungal - lineage, None, "all other fungal ASVs (1,071), all samples")
dist(fungal - lineage, chr_cols, "all other fungal ASVs, Christiana")
top = sorted(((sum(tab[a]), len(seqs[a]), a[:8]) for a in lineage if a in tab), reverse=True)[:10]
print("\nten most abundant lineage ASVs (reads, length, id):")
for n, L, a in top: print(f"   {int(n):>8,}  {L} bp  {a}")
hl = Counter(len(seqs[a]) for a in host)
print("\nhost ASV lengths (unweighted count of ASVs):", dict(sorted(hl.items())))
