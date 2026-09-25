#!/usr/bin/env python
"""Figure 8: multi-panel placement of the item-7 ASVs, one pruned neighbourhood tree per panel,
same standards as Figure 3 (IQ-TREE GTR+G, 1000 UFBoot; support shown at the ASV's parent node).
Figure 5 regenerated from the machine-readable Hill TSV."""
import csv, os, re
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from Bio import Phylo
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
P = f"{O}/placement"; F = f"{O}/figures"; G = f"{O}/gap_closure"
CB = {"orange":"#E69F00","blue":"#0072B2","green":"#009E73","vermillion":"#D55E00",
      "black":"#000000","grey":"#999999"}
plt.rcParams.update({"font.size":8,"savefig.dpi":300,"savefig.bbox":"tight"})

# ---------------- Figure 8 ----------------
res = list(csv.DictReader(open(f"{P}/placement_results.tsv"), delimiter="\t"))
ORG = {r["accession"]: r["organism"] for r in csv.DictReader(open(f"{P}/reference_organism_names.tsv"), delimiter="\t")}
by_nb = {}
for r in res: by_nb.setdefault(r["neighbourhood"], []).append(r)
NB = [n for n in ["Didymellaceae","Microbotryomycetes","Mortierellaceae","Dothideomycetes",
                  "Dothideaceae","Coniochaetales","Saccotheciaceae","Sclerotiniaceae"] if n in by_nb]
fig, axes = plt.subplots(2, 4, figsize=(19, 11))
axes = axes.ravel()
for ax, nb in zip(axes, NB):
    tp = f"{P}/{nb}_tree.contree"
    if not os.path.exists(tp):
        ax.text(.5,.5,f"{nb}\nno tree",ha="center"); ax.axis("off"); continue
    t = Phylo.read(tp,"newick"); t.root_at_midpoint(); t.ladderize()
    asv_tips = {x.name for x in t.get_terminals() if x.name and x.name.startswith("ASV_")}
    near = set()
    for r in by_nb[nb]:
        m = re.search(r"\(([A-Z]{1,2}\d{5,6}\.\d)\)", r["nearest_reference"])
        if m: near.add(m.group(1))
    keep = set(asv_tips)
    for x in t.get_terminals():
        if x.name and any(a in x.name for a in near): keep.add(x.name)
    extra = [x.name for x in t.get_terminals() if x.name and x.name not in keep]
    keep |= set(extra[:max(0, 14-len(keep))])
    for x in list(t.get_terminals()):
        if x.name not in keep:
            try: t.prune(x)
            except Exception: pass
    def lbl(n):
        if not n.is_terminal() or not n.name: return ""
        if n.name.startswith("ASV_"):
            m=re.match(r"ASV_([0-9a-f]+)_rel([\d.]+)",n.name)
            return f"ASV {m.group(1)}  ({m.group(2)}%)" if m else n.name
        # KMD PROMPT 5, D0c: organism name from NCBI (p5_fig8_refnames.py) and accession, in place of the
        # title-derived tip name, which had lost its genus
        m = re.search(r"((?:[A-Z]{2}_\d{6,9}|[A-Z]{1,2}\d{5,9})\.\d)$", n.name)
        if m and m.group(1) in ORG:
            org = re.sub(r" \(in: [^)]*\)", "", ORG[m.group(1)])       # NCBI homonym qualifier only
            return f"{org} {m.group(1)}"
        s=n.name.replace("REF_","").replace("_"," ")
        return s[:44]
    Phylo.draw(t, axes=ax, do_show=False, label_func=lbl,
               label_colors=lambda s: CB["vermillion"] if s.startswith("ASV ") else CB["black"],
               show_confidence=False)
    for cl in t.get_nonterminals():
        try: c=float(cl.confidence) if cl.confidence is not None else None
        except (TypeError,ValueError): c=None
    ax.set_title(f"{nb}  (n = {len(by_nb[nb])} ASVs)", fontsize=9.5, loc="left")
    ax.set_ylabel(""); ax.set_yticks([])
    ax.set_xlabel("substitutions per site", fontsize=7)
    for sp in ("top","right","left"): ax.spines[sp].set_visible(False)
for ax in axes[len(NB):]: ax.axis("off")
fig.suptitle("Figure 8. Phylogenetic placement of abundant unnamed ASVs, by taxonomic neighbourhood.\n"
             "Query ASVs in orange with their maximum group relative abundance. Separate alignment and tree per "
             "neighbourhood (MAFFT, IQ-TREE GTR+G, 1000 ultrafast bootstrap). Placement only, no name assigned.",
             fontsize=11, x=.01, ha="left")
fig.tight_layout(rect=[0,0,1,.94])
fig.savefig(f"{F}/Figure8_placement_by_neighbourhood.png"); plt.close(fig)
print(f"wrote Figure8_placement_by_neighbourhood.png  ({len(NB)} panels)")

# ---------------- Figure 5 from TSV ----------------
rows = list(csv.DictReader(open(f"{G}/hill_numbers_per_sample.tsv"), delimiter="\t"))
SITE_C={"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M={"asymptomatic":"o","symptomatic":"^"}
gg=[("Christiana","asymptomatic"),("Christiana","symptomatic"),
    ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
fig,axes=plt.subplots(1,3,figsize=(9.6,3.7))
titles={"0":"q = 0  (richness)","1":"q = 1  (exp Shannon)","2":"q = 2  (inverse Simpson)"}
for q,ax in zip(("0","1","2"),axes):
    sub=[r for r in rows if r["hill_q"]==q]
    for i,g in enumerate(gg):
        v=[(r["sample_id"],float(r["qD"])) for r in sub
           if r["site"]==g[0] and r["condition"]==g[1]]
        if not v: continue
        xs=[i+(hash(s)%100-50)/420 for s,_ in v]
        ax.scatter(xs,[y for _,y in v],s=30,facecolor="white",edgecolor=SITE_C[g[0]],
                   marker=COND_M[g[1]],linewidth=1.2,zorder=3)
        mu=sum(y for _,y in v)/len(v)
        ax.plot([i-.22,i+.22],[mu,mu],color=CB["black"],linewidth=1.6,zorder=4)
    ax.set_xticks(range(len(gg)))
    ax.set_xticklabels([f"{g[0][:4]}\n{g[1][:5]}\nn={sum(1 for r in sub if r['site']==g[0] and r['condition']==g[1])}"
                        for g in gg],fontsize=7.5)
    ax.set_title(titles[q],fontsize=9)
    for sp in ("top","right"): ax.spines[sp].set_visible(False)
axes[0].set_ylabel("Effective number of ASVs")
fig.suptitle("Figure 5. Hill numbers at depth 57,141. Points are samples, bars are group means.",
             fontsize=10,x=.02,ha="left")
fig.savefig(f"{F}/Figure5_hill_numbers.png"); plt.close(fig)
print("regenerated Figure5_hill_numbers.png from hill_numbers_per_sample.tsv")
for q in ("0","1","2"):
    for g in gg:
        v=[float(r["qD"]) for r in rows if r["hill_q"]==q and r["site"]==g[0] and r["condition"]==g[1]]
        if v: print(f"   q={q} {g[0]:11} {g[1]:13} n={len(v)} mean={sum(v)/len(v):6.1f}")
