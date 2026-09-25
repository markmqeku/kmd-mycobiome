#!/usr/bin/env python
"""17b: Figure S1 on the corrected core-16 table at the re-locked depth 57,124."""
import csv, os, re, glob, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
C = f"{O}/clean"; G = f"{C}/gap_closure"; FC = f"{C}/figures"; os.makedirs(FC, exist_ok=True)
DEPTH = 57124
CB = {"blue":"#0072B2","green":"#009E73","black":"#000000"}
SITE_C = {"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M = {"asymptomatic":"o","symptomatic":"^"}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "savefig.dpi":300,"savefig.bbox":"tight"})
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
f = [p for p in glob.glob(f"{G}/rarefaction_exp/*.csv") if "observed" in os.path.basename(p).lower()]
assert f, "no observed_features csv produced"
rows = list(csv.reader(open(f[0]))); hdr = rows[0]
cols = {}
for i, hh in enumerate(hdr):
    m = re.match(r"depth-(\d+)_iter-\d+", hh)
    if m: cols.setdefault(int(m.group(1)), []).append(i)
depths = sorted(cols)
fig, ax = plt.subplots(figsize=(6.4, 4.2))
n_plot = 0
for r in rows[1:]:
    sid = r[0]
    if sid not in meta: continue
    xs, ys = [], []
    for d in depths:
        vals = [float(r[i]) for i in cols[d] if i < len(r) and r[i].strip()]
        if vals: xs.append(d); ys.append(st.mean(vals))
    if not xs: continue
    ax.plot(xs, ys, linewidth=1.1, alpha=.85,
            color=SITE_C.get(meta[sid]["location"], CB["black"]),
            marker=COND_M.get(meta[sid]["condition_std"], "o"), markersize=3, markevery=4)
    n_plot += 1
ax.axvline(DEPTH, color=CB["black"], linestyle="--", linewidth=1)
ax.annotate(f"rarefaction depth\n{DEPTH:,}", (DEPTH, ax.get_ylim()[1]*0.25),
            textcoords="offset points", xytext=(-64, 0), fontsize=7.5)
ax.set_xlabel("Sequencing depth (reads)"); ax.set_ylabel("Observed ASVs")
ax.set_title(f"Figure S1. Rarefaction curves, {n_plot} core samples (corrected fungi-only table)",
             fontsize=10, loc="left")
leg = [Line2D([0],[0],color=SITE_C[s],lw=1.5,label=s) for s in ("Christiana","Pretoria")]
leg += [Line2D([0],[0],color=CB["black"],marker=COND_M[c],lw=0,label=c) for c in COND_M]
ax.legend(handles=leg, frameon=False, fontsize=8)
fig.savefig(f"{FC}/FigureS1_rarefaction_core16.png"); plt.close(fig)
print(f"wrote {FC}/FigureS1_rarefaction_core16.png  ({n_plot} samples, depth marked at {DEPTH:,})")
