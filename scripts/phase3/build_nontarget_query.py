#!/usr/bin/env python
"""T-4a: assemble the non-target screen query set.
   (i)  the 78 ASVs called Homophron under the permissive vsearch 0.80 run
   (ii) every ASV unassigned at kingdom level under the final id97 run
Writes a FASTA sorted by descending total reads, and reports read coverage."""
import csv, os
B = "<KMD_ROOT>/__reanalysis_2026-06"
OUT = f"{B}/phase3/outputs_27-07-2026/nontarget_screen"
os.makedirs(OUT, exist_ok=True)

def load_tax(p):
    d = {}
    for r in csv.reader(open(p), delimiter="\t"):
        if r and r[0] not in ("Feature ID",) and not r[0].startswith("#"):
            d[r[0]] = r[1]
    return d

tax80 = load_tax(f"{B}/phase3/outputs_27-07-2026/taxonomy_unite10/exp/taxonomy.tsv")
tax97 = load_tax(f"{B}/phase3/outputs_27-07-2026/taxonomy_unite10_id97/exp/taxonomy.tsv")

def has_kingdom(t):
    if not t or t.lower().startswith("unassigned"):
        return False
    for part in t.split(";"):
        part = part.strip()
        if part.startswith("k__"):
            v = part[3:].strip().lower()
            return v not in ("", "unidentified", "unclassified")
    return False

homophron = {a for a, t in tax80.items() if "Homophron" in t}
no_kingdom = {a for a in tax97 if not has_kingdom(tax97[a])}
query = homophron | no_kingdom

# abundances
tabp = f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv"
lines = [l.rstrip("\n") for l in open(tabp) if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); samples = [h for h in hdr[1:] if h.strip()]
counts, total = {}, {}
for l in lines[1:]:
    p = l.split("\t")
    c = [int(float(x)) for x in p[1:1+len(samples)]]
    counts[p[0]] = c; total[p[0]] = sum(c)
grand = sum(total.values())

# sequences
seqs, sid, buf = {}, None, []
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line.strip())
if sid: seqs[sid] = "".join(buf)

q = sorted(query, key=lambda a: -total.get(a, 0))
qreads = sum(total.get(a, 0) for a in q)
print(f"Homophron(0.80) ASVs        : {len(homophron)}")
print(f"kingdom-unassigned(id97)    : {len(no_kingdom)}")
print(f"UNION query set             : {len(q)} ASVs")
print(f"reads in query set          : {qreads:,} / {grand:,} = {100*qreads/grand:.1f}% of dataset")

cum = 0
for i, a in enumerate(q, 1):
    cum += total.get(a, 0)
    if cum >= 0.95 * qreads:
        print(f"top {i} ASVs cover 95% of query-set reads")
        break

with open(f"{OUT}/nontarget_query_all.fasta", "w") as fh:
    for a in q:
        fh.write(f">{a}\n{seqs.get(a,'')}\n")
# priority subset: enough to cover 99% of query reads, capped for remote BLAST throughput
cum, pri = 0, []
for a in q:
    pri.append(a); cum += total.get(a, 0)
    if cum >= 0.99 * qreads and len(pri) >= 50: break
pri = pri[:250]
with open(f"{OUT}/nontarget_query_priority.fasta", "w") as fh:
    for a in pri:
        fh.write(f">{a}\n{seqs.get(a,'')}\n")
print(f"priority FASTA: {len(pri)} ASVs = {100*sum(total.get(a,0) for a in pri)/grand:.1f}% of dataset reads")
with open(f"{OUT}/query_abundance.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["ASV_ID_md5","total_reads","in_homophron78","kingdom_unassigned_id97"])
    for a in q: w.writerow([a, total.get(a,0), int(a in homophron), int(a in no_kingdom)])
print("wrote:", OUT)
