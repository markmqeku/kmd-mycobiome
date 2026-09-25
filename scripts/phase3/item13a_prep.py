#!/usr/bin/env python
"""13a prep. Stratify all 1,173 fungal ASVs by taxonomic resolution and abundance,
then write tiered FASTA query files for a whole-table NCBI nt screen.
Tier 1 is the highest-risk set (unassigned / kingdom-only); tiers run in order so the
riskiest sequences are screened first even if the remote BLAST is slow."""
import csv, os
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
S = f"{O}/full_screen"; os.makedirs(S, exist_ok=True)
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
ORDER = [r for _, r in RANKS]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

def deepest(t):
    if not t or t.lower().startswith("unassigned"): return "none"
    last = "none"
    for part in t.split(";"):
        part = part.strip()
        for pre, rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk == "species" and (b.endswith("_sp") or b.endswith("_sp.")))):
                    continue
                last = rk
    return last

tax = {}
for r in csv.reader(open(f"{O}/final_taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"):
        tax[r[0]] = r[1]

# read counts from the fungi-only table
lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
reads, insamp = {}, {}
for l in lines[1:]:
    q = l.split("\t")
    vals = [int(float(v)) for v in q[1:1+len(cols)]]
    reads[q[0]] = sum(vals)
    insamp[q[0]] = [s for s, v in zip(cols, vals) if v]

seqs, sid, buf = {}, None, []
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line)
if sid: seqs[sid] = "".join(buf)

TIER = {"none":1,"kingdom":1,"phylum":2,"class":2,"order":3,"family":3,"genus":4,"species":4}
rows = []
for a in seqs:
    d = deepest(tax.get(a, ""))
    rows.append((TIER[d], -reads.get(a, 0), a, d, reads.get(a, 0), len(insamp.get(a, []))))
rows.sort()

c = Counter(r[3] for r in rows)
print("deepest resolved rank across the 1,173 fungal ASVs:")
for r in ["none"] + ORDER:
    if c[r]: print(f"   {r:9} {c[r]:>5} ASVs   {sum(x[4] for x in rows if x[3]==r):>9,} reads")
print()
for t in (1, 2, 3, 4):
    sub = [r for r in rows if r[0] == t]
    with open(f"{S}/tier{t}.fasta", "w") as fh:
        for _, _, a, d, rd, ns in sub:
            fh.write(f">{a}\n{seqs[a]}\n")
    print(f"tier {t}: {len(sub):>5} ASVs, {sum(x[4] for x in sub):>9,} reads -> tier{t}.fasta")

with open(f"{S}/asv_index.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["asv_id","tier","deepest_rank","final_taxon","reads","n_samples","samples","length"])
    for t, _, a, d, rd, ns in rows:
        w.writerow([a, t, d, tax.get(a, ""), rd, ns, ";".join(insamp.get(a, [])), len(seqs[a])])
print(f"\nindex -> {S}/asv_index.tsv")
