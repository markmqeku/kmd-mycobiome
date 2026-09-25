#!/usr/bin/env python
"""64g. Figure 9: abundance distribution of all ASVs with the 1 percent threshold marked, and
assignment rate by abundance class. Same standards as the existing figures: Okabe-Ito palette,
individual points, no boxplots or violins."""
import csv, json, os, random, zipfile
from collections import Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"; F = f"{C}/figures"
os.makedirs(F, exist_ok=True)
DEPTH, THRESH = 57124, 0.01
CB = {"orange":"#E69F00","blue":"#0072B2","green":"#009E73","vermillion":"#D55E00",
      "purple":"#CC79A7","black":"#000000","grey":"#999999"}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.dpi":300,"savefig.dpi":300,"savefig.bbox":"tight"})

meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
lines = [l.rstrip("\n") for l in open(f"{C}/feature-table-clean.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [x for x in hdr[1:] if x.strip()]
per = {s: {} for s in cols}
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]
rng = random.Random(42)
def rarefy(c, d):
    pool = []
    for a, v in c.items(): pool.extend([a]*v)
    return Counter(rng.sample(pool, d))
rar = {s: rarefy(per[s], DEPTH) for s in core}

RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}
def ranks_of(t):
    out = {r: "" for _, r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part = part.strip()
        for pre, rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk == "species" and (b.endswith("_sp") or b.endswith("_sp.")))): v = ""
                out[rk] = v
    return out
final = {}
for r in csv.reader(open(f"{O}/final_taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): final[r[0]] = r[1]

# maximum relative abundance each ASV reaches in any core sample
maxrel = {}
for s in core:
    for a, v in rar[s].items():
        maxrel[a] = max(maxrel.get(a, 0.0), v/DEPTH)
asvs = sorted(maxrel, key=lambda a: -maxrel[a])

fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))

# --- panel a: rank abundance, every ASV as one point ---
ax = axes[0]
xs = list(range(1, len(asvs)+1))
ys = [100*maxrel[a] for a in asvs]
colf = [CB["vermillion"] if maxrel[a] >= THRESH else CB["blue"] for a in asvs]
ax.scatter(xs, ys, s=7, facecolor="none", edgecolor=colf, linewidth=0.7, zorder=3)
ax.axhline(1.0, color=CB["black"], linestyle="--", linewidth=1.1, zorder=4)
ax.set_yscale("log")
n_ab = sum(1 for a in asvs if maxrel[a] >= THRESH)
ax.annotate(f"1 percent threshold\n{n_ab} ASVs above, {len(asvs)-n_ab} below",
            (len(asvs)*0.42, 1.35), fontsize=8)
ax.set_xlabel("ASVs, ranked by maximum relative abundance in any core sample")
ax.set_ylabel("Maximum relative abundance (%, log scale)")
ax.set_title("a. Abundance distribution of all ASVs in the core 16", fontsize=10, loc="left")

# --- panel b: assignment rate by abundance class ---
BINS = [(0.0, 0.01, "below\n0.01%"), (0.01, 0.1, "0.01 to\n0.1%"),
        (0.1, 1.0, "0.1 to\n1%"), (1.0, 10.0, "1 to\n10%"), (10.0, 101.0, "above\n10%")]
ax = axes[1]
for rank, colour, mark in (("family", CB["green"], "o"), ("genus", CB["purple"], "^")):
    xs2, ys2, ns = [], [], []
    for i, (lo, hi, lab) in enumerate(BINS):
        grp = [a for a in asvs if lo <= 100*maxrel[a] < hi]
        if not grp: continue
        rate = 100*sum(1 for a in grp if ranks_of(final.get(a, ""))[rank])/len(grp)
        xs2.append(i); ys2.append(rate); ns.append(len(grp))
    ax.plot(xs2, ys2, color=colour, linewidth=1.2, zorder=2, alpha=.85)
    ax.scatter(xs2, ys2, s=46, facecolor="white", edgecolor=colour, marker=mark,
               linewidth=1.5, zorder=3, label=f"resolved to {rank}")
for i, (lo, hi, lab) in enumerate(BINS):
    grp = [a for a in asvs if lo <= 100*maxrel[a] < hi]
    if grp: ax.annotate(f"n={len(grp)}", (i, 4.5), ha="center", fontsize=7.5, color=CB["grey"])
ax.axvline(2.5, color=CB["black"], linestyle="--", linewidth=1.1)
ax.annotate("1 percent", (2.55, 92), fontsize=8)
ax.set_xticks(range(len(BINS))); ax.set_xticklabels([b[2] for b in BINS], fontsize=8)
ax.set_ylim(0, 100); ax.set_ylabel("ASVs resolved to the rank (%)")
ax.set_xlabel("Maximum relative abundance class")
ax.set_title("b. Taxonomic resolution by abundance class", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="upper left", bbox_to_anchor=(0.0, 0.86))   # KMD PROMPT 4, R6b: clear of the n labels

fig.suptitle("Figure 9. The low-abundance fraction of the fungal community, 16 core samples "
             "rarefied to 57,124 reads.\nEach point in panel a is one ASV; no distribution "
             "summaries are used.", fontsize=10, x=.01, ha="left")
fig.tight_layout(rect=[0, 0, 1, 0.93])
fig.savefig(f"{F}/Figure9_rare_fraction.png"); plt.close(fig)
print(f"wrote {F}/Figure9_rare_fraction.png")
print(f"   ASVs plotted: {len(asvs)}   at or above 1%: {n_ab}   below: {len(asvs)-n_ab}")
for lo, hi, lab in BINS:
    grp = [a for a in asvs if lo <= 100*maxrel[a] < hi]
    if grp:
        g = 100*sum(1 for a in grp if ranks_of(final.get(a,""))["genus"])/len(grp)
        f_ = 100*sum(1 for a in grp if ranks_of(final.get(a,""))["family"])/len(grp)
        print(f"   {lab.replace(chr(10),' '):16} n={len(grp):>4}  family {f_:5.1f}%  genus {g:5.1f}%")
