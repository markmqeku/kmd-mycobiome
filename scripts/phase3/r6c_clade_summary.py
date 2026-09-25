#!/usr/bin/env python
"""KMD PROMPT 4, R6c. Summary values from the rerun Figure 3 tree (fig3_71_tree.contree, rooted on the four
Curvibasidium references) and from the original 78-ASV tree under the same rooting, for FINAL_NUMBERS.
Read counts are taken from the tip labels (ASV_<md5>_reads<n>). Run with the qiime2 environment's Python."""
import json, re
from Bio import Phylo

O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"


def rd(p):
    s, k = {}, None
    for l in open(p):
        l = l.rstrip("\n")
        if l.startswith(">"): k = l[1:].split()[0]; s[k] = ""
        elif k: s[k] += l.strip()
    return s


def reads(n):
    return int(re.search(r"_reads(\d+)", n).group(1))


def rooted(path):
    t = Phylo.read(path, "newick")
    names = [c.name for c in t.get_terminals()]
    curv = [n for n in names if n.startswith("REF_Curvibasidium")]
    t.root_with_outgroup(t.common_ancestor(*curv))
    trem = [n for n in names if n.startswith("REF_") and n not in curv]
    tm = t.common_ancestor(*trem)
    return t, names, curv, tm, {c.name for c in tm.get_terminals()}


out = {}
t, names, curv, tm, inside = rooted(f"{O}/yeast_placement_71/fig3_71_tree.contree")
asvs = [n for n in names if n.startswith("ASV_")]
tot = sum(reads(a) for a in asvs)
un = [n for n in names if n.startswith("UNNAMED_")][0]
node = [c for c in t.get_terminals() if c.name == un][0]
path = t.get_path(node)
lin = None
for anc in reversed(path[:-1]):                      # largest ancestor of MK019123.1 holding no reference
    if any(c.name.startswith("REF_") for c in anc.get_terminals()): break
    lin = anc
lin_asvs = [c.name for c in lin.get_terminals() if c.name.startswith("ASV_")]
aln = rd(f"{O}/yeast_placement_71/fig3_71_aln.fasta")
out["rerun"] = {
    "tips": len(names), "asvs": len(asvs), "alignment_columns": len(next(iter(aln.values()))),
    "tremellomycetes_clade_support": tm.confidence, "asvs_inside": sum(a in inside for a in asvs),
    "asvs_outside": [a for a in asvs if a not in inside],
    "curvibasidium_inside": sum(c in inside for c in curv),
    "mk019123_clade_asvs": len(lin_asvs), "mk019123_clade_support": lin.confidence,
    "mk019123_clade_reads": sum(reads(a) for a in lin_asvs), "lineage_reads_in_tree": tot,
    "mk019123_clade_read_pct": round(100 * sum(reads(a) for a in lin_asvs) / tot, 2),
    "most_abundant_in_mk_clade": max(asvs, key=reads) in lin_asvs,
}
# the original 78-ASV tree, same rooting and criterion
t0, names0, curv0, tm0, inside0 = rooted(f"{O}/yeast_placement/fig3_tree.contree")
kept = {v.upper() for v in rd(f"{O}/clean/genbank_submission/Tremellomycetes_ASVs.fasta").values()}
inp = rd(f"{O}/yeast_placement/fig3_input.fasta")
a0 = [n for n in names0 if n.startswith("ASV_")]
alg = [a for a in a0 if inp[a].upper() not in kept]
out["original_78"] = {"tremellomycetes_clade_support": tm0.confidence,
                      "retained_inside": sum(a in inside0 for a in a0 if a not in alg), "retained": len(a0) - len(alg),
                      "green_algal_inside": sum(a in inside0 for a in alg), "green_algal": len(alg)}
json.dump(out, open(f"{O}/yeast_placement_71/r6c_clade_summary.json", "w"), indent=1)
print(json.dumps(out, indent=1))
