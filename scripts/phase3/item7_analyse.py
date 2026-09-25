#!/usr/bin/env python
"""Item 7 c-f: per-ASV nearest reference (bitscore-ranked), placement from the neighbourhood
tree, and node support. Reports PLACEMENT, never identification."""
import csv, os, re, glob
from collections import defaultdict
from Bio import Phylo
P = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/placement"
NEIGH = ["Didymellaceae","Microbotryomycetes","Mortierellaceae","Dothideomycetes",
         "Coniochaetales","Dothideaceae","Saccotheciaceae","Sclerotiniaceae"]
rows = []
for nb in NEIGH:
    bl = f"{P}/{nb}_blast.tsv"; tr = f"{P}/{nb}_tree.contree"
    if not os.path.exists(bl): print(f"{nb}: no blast"); continue
    best = {}
    for r in csv.reader(open(bl), delimiter="\t"):
        if len(r) < 7: continue
        q, s, pid, alen, qcov, ev, bits = r[0], r[1], float(r[2]), int(r[3]), float(r[4]), r[5], float(r[6])
        if q not in best or bits > best[q][5]: best[q] = (s, pid, alen, qcov, ev, bits)
    tree = None; support = {}
    if os.path.exists(tr):
        try:
            tree = Phylo.read(tr, "newick"); tree.root_at_midpoint()
            for tip in tree.get_terminals():
                if tip.name and tip.name.startswith("ASV_"):
                    path = tree.get_path(tip)
                    conf = None
                    for cl in reversed(path[:-1]):
                        if cl.confidence is not None:
                            try: conf = float(cl.confidence); break
                            except (TypeError, ValueError): pass
                    sisters = []
                    if len(path) >= 2:
                        parent = path[-2]
                        sisters = [t.name for t in parent.get_terminals() if t.name != tip.name]
                    support[tip.name] = (conf, sisters[:3])
        except Exception as e:
            print(f"{nb}: tree read failed {type(e).__name__}")
    print(f"\n=== {nb} ===")
    for q in sorted(best, key=lambda x: -float(re.search(r"rel([\d.]+)", x).group(1)) if re.search(r"rel([\d.]+)", x) else 0):
        s, pid, alen, qcov, ev, bits = best[q]
        conf, sis = support.get(q, (None, []))
        refname = s.replace("REF_", "").replace("_", " ")
        sis_s = "; ".join(x.replace("REF_", "").replace("_", " ")[:34] for x in sis) or "none resolved"
        # placement wording rules
        if conf is not None and conf >= 70 and sis:
            place = f"placed within {nb}, sister to {sis_s}"
        elif sis:
            place = f"placed within {nb}, position weakly supported"
        else:
            place = f"placed within {nb}, no resolved sister group"
        print(f"  {q:26} nearest={refname[:40]:40} id={pid:5.1f}% aln={alen:4} "
              f"support={'-' if conf is None else int(conf)}")
        rows.append([nb, q, refname, f"{pid:.1f}", alen, f"{qcov:.0f}", ev,
                     "" if conf is None else f"{conf:.0f}", place, sis_s])
with open(f"{P}/placement_results.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["neighbourhood","ASV","nearest_reference","identity_pct","alignment_len",
                "query_cov_pct","evalue","node_support_UFBoot","placement","sister_taxa"])
    w.writerows(rows)
print(f"\nwrote {P}/placement_results.tsv  ({len(rows)} ASVs)")
hi = sum(1 for r in rows if r[7] and float(r[7]) >= 70)
g95 = sum(1 for r in rows if float(r[3]) >= 95)
print(f"  ASVs with node support >=70 UFBoot : {hi}/{len(rows)}")
print(f"  ASVs whose nearest reference is >=95% identical: {g95}/{len(rows)}")
print("  NOTE: no genus or species is assigned from topology. Placement wording only.")
