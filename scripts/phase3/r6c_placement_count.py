#!/usr/bin/env python
"""KMD PROMPT 4, R6c. Which of the 71 lineage ASVs place within Tremellomycetes in the rerun tree?

The tree (IQ-TREE consensus, fig3_71_tree.contree) is rooted on the Curvibasidium references
(Microbotryomycetes, outside Tremellomycetes). An ASV places within Tremellomycetes if it lies inside the
smallest clade containing all 48 Tremellomycetes references. Run with the qiime2 environment's Python
(Biopython). Writes a TSV of every ASV's placement and prints the count."""
import csv, re
from Bio import Phylo

D = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/yeast_placement_71"
t = Phylo.read(f"{D}/fig3_71_tree.contree", "newick")
tips = {c.name: c for c in t.get_terminals()}
curv = [n for n in tips if n.startswith("REF_Curvibasidium")]
trem = [n for n in tips if n.startswith("REF_") and not n.startswith("REF_Curvibasidium")]
asvs = [n for n in tips if n.startswith("ASV_")]
print(f"tips: {len(tips)}  ASVs {len(asvs)}  Tremellomycetes refs {len(trem)}  Curvibasidium refs {len(curv)}")
t.root_with_outgroup(*curv)
out_mrca = t.common_ancestor(*curv)
print("Curvibasidium references monophyletic after rooting:",
      set(c.name for c in out_mrca.get_terminals()) == set(curv))
tm = t.common_ancestor(*trem)
inside = {c.name for c in tm.get_terminals()}
extra = sorted(n for n in inside if not n.startswith(("REF_", "ASV_", "UNNAMED_")))
print("Tremellomycetes clade: bootstrap", tm.confidence, "| tips", len(inside),
      "| contains Curvibasidium refs:", any(n in inside for n in curv))
placed = [a for a in asvs if a in inside]
outside = [a for a in asvs if a not in inside]
print(f"ASVs placed within Tremellomycetes: {len(placed)} of {len(asvs)}")
print("outside:", outside)
print("MK019123.1 inside:", any(n.startswith("UNNAMED_") for n in inside))
# where does each outside ASV sit? report its sister group
for a in outside:
    node = tips[a]
    path = t.get_path(node)
    parent = path[-2] if len(path) > 1 else t.root
    sis = [c.name for c in parent.get_terminals() if c.name != a]
    print(f"   {a}: sister group of {len(sis)} tips, e.g. {sis[:4]}; support {parent.confidence}")
# the main lineage clade: smallest clade containing the ASVs that sit inside Tremellomycetes with MK019123.1
with open(f"{D}/fig3_71_placement.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["asv", "reads", "placed_within_Tremellomycetes"])
    for a in asvs:
        m = re.match(r"ASV_([0-9a-f]+)_reads(\d+)", a)
        w.writerow([m.group(1), m.group(2), "yes" if a in inside else "no"])
aln = [l for l in open(f"{D}/fig3_71_aln.fasta") if l.startswith(">")]
seq = "".join(open(f"{D}/fig3_71_aln.fasta").read().split(">")[1].split("\n")[1:])
print("alignment: sequences", len(aln), "| columns", len(seq))
