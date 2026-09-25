#!/usr/bin/env python
"""17b. Draw Figures 5 and 6 on the corrected table at the re-locked depth 57,124 (DEC-1).
These were blocked while the depth was undecided. Same code paths and same seed as
make_figures.py / make_fig8_fig5.py, so old and new differ only by the correction."""
import csv, math, os, random
from collections import Counter
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"; FC = f"{C}/figures"
os.makedirs(FC, exist_ok=True)
DEPTH = 57124
CB = {"orange":"#E69F00","skyblue":"#56B4E9","green":"#009E73","yellow":"#F0E442",
      "blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7","black":"#000000","grey":"#999999"}
SITE_C = {"Bloemfontein":CB["orange"],"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M = {"reference_site":"s","asymptomatic":"o","symptomatic":"^"}
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
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
assert len(core) == 16, f"expected 16 core samples, found {len(core)}"
assert min(sum(per[s].values()) for s in core) >= DEPTH, "a core sample is below the depth"
gg = [("Christiana","asymptomatic"),("Christiana","symptomatic"),
      ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]

# ---------------- Figure 5 ----------------
rows = list(csv.DictReader(open(f"{C}/hill_numbers_per_sample.tsv"), delimiter="\t"))
fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.7))
titles = {"0":"q = 0  (richness)","1":"q = 1  (exp Shannon)","2":"q = 2  (inverse Simpson)"}
for q, ax in zip(("0","1","2"), axes):
    sub = [r for r in rows if r["hill_q"] == q]
    for i, g in enumerate(gg):
        v = [(r["sample_id"], float(r["qD"])) for r in sub
             if r["site"] == g[0] and r["condition"] == g[1]]
        if not v: continue
        ax.scatter([i+(hash(s) % 100-50)/420 for s, _ in v], [y for _, y in v], s=30,
                   facecolor="white", edgecolor=SITE_C[g[0]], marker=COND_M[g[1]],
                   linewidth=1.2, zorder=3)
        mu = sum(y for _, y in v)/len(v)
        ax.plot([i-.22, i+.22], [mu, mu], color=CB["black"], linewidth=1.6, zorder=4)
    ax.set_xticks(range(len(gg)))
    ax.set_xticklabels([f"{g[0][:4]}\n{g[1][:5]}\nn={sum(1 for r in sub if r['site']==g[0] and r['condition']==g[1])}"
                        for g in gg], fontsize=7.5)
    ax.set_title(titles[q], fontsize=9)
axes[0].set_ylabel("Effective number of ASVs")
fig.suptitle(f"Figure 5. Hill numbers at depth {DEPTH:,}. Points are samples, bars are group means.",
             fontsize=10, x=.02, ha="left")
fig.savefig(f"{FC}/Figure5_hill_numbers.png"); plt.close(fig)
print(f"wrote {FC}/Figure5_hill_numbers.png")

# ---------------- Figure 6 ----------------
rng = random.Random(42)
def rarefy(c, d):
    pool = []
    for a, v in c.items(): pool.extend([a]*v)
    return Counter(pool if len(pool) < d else rng.sample(pool, d))
rar = {s: rarefy(per[s], DEPTH) for s in core}
def bray(a, b):
    ks = set(a) | set(b); num = sum(abs(a.get(k,0)-b.get(k,0)) for k in ks)
    den = sum(a.values())+sum(b.values()); return num/den if den else 0.
def jac(a, b):
    A, Bs = set(a), set(b); return 1-len(A & Bs)/len(A | Bs) if (A | Bs) else 0.
def pcoa(D, n):
    A = [[-0.5*D[i][j]**2 for j in range(n)] for i in range(n)]
    rm = [sum(r)/n for r in A]; cm = [sum(A[i][j] for i in range(n))/n for j in range(n)]
    gm = sum(rm)/n
    G = [[A[i][j]-rm[i]-cm[j]+gm for j in range(n)] for i in range(n)]
    def ax_(G):
        v = [rng.random() for _ in range(n)]
        for _ in range(600):
            w = [sum(G[i][j]*v[j] for j in range(n)) for i in range(n)]
            nn = math.sqrt(sum(x*x for x in w)) or 1; v = [x/nn for x in w]
        lam = sum(v[i]*sum(G[i][j]*v[j] for j in range(n)) for i in range(n)); return lam, v
    l1, v1 = ax_(G); G2 = [[G[i][j]-l1*v1[i]*v1[j] for j in range(n)] for i in range(n)]
    l2, v2 = ax_(G2); tot = sum(G[i][i] for i in range(n))
    return ([v1[i]*math.sqrt(abs(l1)) for i in range(n)],
            [v2[i]*math.sqrt(abs(l2)) for i in range(n)], 100*l1/tot, 100*l2/tot)
fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.2))
for ax, (nm, fnc) in zip(axes, [("Bray-Curtis", bray), ("Jaccard", jac)]):
    n = len(core); D = [[fnc(rar[core[i]], rar[core[j]]) for j in range(n)] for i in range(n)]
    x, y, p1, p2 = pcoa(D, n)
    for i, s in enumerate(core):
        ax.scatter(x[i], y[i], s=52, facecolor="white", edgecolor=SITE_C[meta[s]["location"]],
                   marker=COND_M[meta[s]["condition_std"]], linewidth=1.4, zorder=3)
    ax.axhline(0, color="#EEEEEE", zorder=0); ax.axvline(0, color="#EEEEEE", zorder=0)
    ax.set_xlabel(f"PCo1 ({p1:.1f}%)"); ax.set_ylabel(f"PCo2 ({p2:.1f}%)"); ax.set_title(nm, fontsize=9)
    print(f"   {nm}: PCo1 {p1:.1f}%, PCo2 {p2:.1f}%")
leg = [Line2D([0],[0],marker="o",color="w",markerfacecolor="w",markeredgecolor=SITE_C[s],
              markeredgewidth=1.4,markersize=7,label=s) for s in ("Christiana","Pretoria")]
leg += [Line2D([0],[0],marker=COND_M[c],color="w",markerfacecolor="w",markeredgecolor=CB["black"],
                markeredgewidth=1.2,markersize=7,label=c) for c in ("asymptomatic","symptomatic")]
axes[1].legend(handles=leg, frameon=False, fontsize=8, bbox_to_anchor=(1.02,1), loc="upper left")
fig.suptitle("Figure 6. Ordination of the 16 core samples. Colour = site, symbol = condition.",
             fontsize=10, x=.02, ha="left")
fig.savefig(f"{FC}/Figure6_ordination.png"); plt.close(fig)
print(f"wrote {FC}/Figure6_ordination.png")
