#!/usr/bin/env python
"""KMD PROMPT 4, R6c. Figure 3 redrawn from the rerun placement tree (71 lineage ASVs, the green-algal ASVs
excluded; r6c_fig3_rerun.sh). IQ-TREE GTR+G consensus tree, 1000 ultrafast bootstrap replicates, rooted on
the four Curvibasidium references (Microbotryomycetes), which fall outside every Tremellomycetes reference
after rooting. Study ASVs, named references with their accessions, and the unnamed match MK019123.1 are
drawn in different colours; ultrafast bootstrap values of 70 or more are printed at their nodes; a scale
bar gives substitutions per site. No title is drawn: the caption is in the manuscript. Run with the qiime2
environment's Python (Biopython, matplotlib)."""
import os, re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from Bio import Phylo

O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
D = f"{O}/yeast_placement_71"
F = f"{O}/figures"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})

t = Phylo.read(f"{D}/fig3_71_tree.contree", "newick")
names = [c.name for c in t.get_terminals()]
curv = [n for n in names if n.startswith("REF_Curvibasidium")]
t.root_with_outgroup(t.common_ancestor(*curv))
t.ladderize(reverse=True)
trem = [n for n in names if n.startswith("REF_") and n not in curv]
trem_clade = t.common_ancestor(*trem)
inside = {c.name for c in trem_clade.get_terminals()}
asvs = [n for n in names if n.startswith("ASV_")]
n_in = sum(a in inside for a in asvs)
assert not any(c in inside for c in curv)

VERM, BLUE, BLACK, GREY = "#D55E00", "#0072B2", "#000000", "#777777"

# rectangular layout: x = distance from the root, y = tip order
depth = t.depths()
if not max(depth.values()):
    depth = t.depths(unit_branch_lengths=True)
tips = t.get_terminals()
ypos = {c: i for i, c in enumerate(tips)}
for cl in t.find_clades(order="postorder"):
    if cl not in ypos:
        ypos[cl] = (ypos[cl.clades[0]] + ypos[cl.clades[-1]]) / 2

H = 9.2; W = 6.85                                 # 174 x 234 mm, full page width
fig, ax = plt.subplots(figsize=(W, H))
LW = 0.5
for cl in t.find_clades(order="preorder"):
    x1 = depth[cl]
    for ch in cl.clades:
        ax.plot([x1, x1], [ypos[cl], ypos[ch]], color=BLACK, lw=LW, solid_capstyle="butt")
        ax.plot([x1, depth[ch]], [ypos[ch], ypos[ch]], color=BLACK, lw=LW, solid_capstyle="butt")
xmax = max(depth.values())
FS = 4.6
for c in tips:
    x, y = depth[c] + xmax * 0.008, ypos[c]
    n = c.name
    if n.startswith("ASV_"):
        m = re.match(r"ASV_([0-9a-f]+)_reads(\d+)", n)
        ax.text(x, y, f"ASV {m.group(1)} ({int(m.group(2)):,} reads)", color=VERM, fontsize=FS, va="center")
    elif n.startswith("UNNAMED_"):
        ax.text(x, y, "MK019123.1 (unnamed environmental sequence, ‘Cryptococcus sp. OTU796’)",
                color=BLUE, fontsize=FS + 0.6, va="center", fontweight="bold")
    else:
        genus, acc = n.split("_")[1], n.split("_")[-1]
        ax.text(x, y, f"$\\it{{{genus}}}$ {acc}", color=BLACK, fontsize=FS, va="center")
# ultrafast bootstrap values of 70 or more, left of the node
for cl in t.get_nonterminals():
    if cl is t.root or cl.confidence is None: continue
    if float(cl.confidence) >= 70:
        ax.text(depth[cl] - xmax * 0.004, ypos[cl] + 0.35, f"{int(round(float(cl.confidence)))}",
                fontsize=3.8, color=GREY, ha="right", va="bottom")
# brackets: the Tremellomycetes-reference clade and the Curvibasidium outgroup
def bracket(clade_tips, label, xb):
    ys = [ypos[c] for c in tips if c.name in clade_tips]
    ax.plot([xb, xb], [min(ys) - 0.3, max(ys) + 0.3], color=BLACK, lw=0.8)
    ax.text(xb + xmax * 0.012, (min(ys) + max(ys)) / 2, label, rotation=270, fontsize=6, va="center", ha="left")
XB = xmax * 1.62
bracket(inside, f"Clade containing all Tremellomycetes references (ultrafast bootstrap {int(float(trem_clade.confidence))}); "
                f"{n_in} of {len(asvs)} study ASVs", XB)
bracket(set(curv), "Outgroup", XB)
# scale bar
sb = 0.5 if xmax > 2 else 0.1
y0 = len(tips) + 1.5
ax.plot([xmax * 0.3, xmax * 0.3 + sb], [y0, y0], color=BLACK, lw=0.8)
ax.text(xmax * 0.3 + sb + xmax * 0.01, y0, f"{sb} substitutions per site", fontsize=5.5, ha="left", va="center")
# key
from matplotlib.lines import Line2D
key = [Line2D([], [], color=VERM, lw=0, marker="s", ms=4, label=f"Study ASVs of the unresolved lineage (n = {len(asvs)})"),
       Line2D([], [], color=BLACK, lw=0, marker="s", ms=4, label="Named reference sequences (genus, GenBank accession)"),
       Line2D([], [], color=BLUE, lw=0, marker="s", ms=4, label="Unnamed nearest match in NCBI nt")]
ax.legend(handles=key, loc="lower right", bbox_to_anchor=(0.9, 0.08), frameon=False, fontsize=5.5)
ax.set_xlim(-xmax * 0.01, xmax * 1.75)
ax.set_ylim(len(tips) + 3, -4)
ax.axis("off")
fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.005)
fig.savefig(f"{F}/Figure3_tremellomycetes_placement_71.png", dpi=300)
plt.close(fig)
sup = [float(c.confidence) for c in t.get_nonterminals() if c.confidence is not None]
print(f"wrote Figure3_tremellomycetes_placement_71.png | tips {len(tips)} | ASVs inside {n_in}/{len(asvs)} | "
      f"Tremellomycetes clade support {trem_clade.confidence} | max root-to-tip {xmax:.3f} | nodes >=70: "
      f"{sum(s >= 70 for s in sup)} of {len(sup)}")
