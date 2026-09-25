#!/usr/bin/env python
"""13d. Regenerate ONLY the figures the change audit flagged, from the corrected table.
Styling is kept identical to make_figures.py (Okabe-Ito, 300 dpi, individual points only,
no distribution summaries). Output goes to clean/figures/ so the originals survive for
side-by-side comparison; nothing under figures/ is overwritten."""
import csv, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"; FC = f"{C}/figures"
os.makedirs(FC, exist_ok=True)
CB = {"orange":"#E69F00","skyblue":"#56B4E9","green":"#009E73","yellow":"#F0E442",
      "blue":"#0072B2","vermillion":"#D55E00","purple":"#CC79A7","black":"#000000","grey":"#999999"}
TISSUE_M = {"Leaves":"o","Twigs":"^","Inflorescence":"s","Seeds":"D","Malformation":"P"}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "figure.dpi":300,"savefig.dpi":300,"savefig.bbox":"tight"})
changed = {r["figure"] for r in csv.DictReader(open(f"{C}/figure_change_audit.tsv"), delimiter="\t")
           if r["plotted_values_change"] == "YES"}
print(f"figures flagged as changed: {sorted(changed)}")
drop = {l.strip() for l in open(f"{O}/full_screen/nonfungal_asvs.txt") if l.strip()}

def load_table(p):
    lines = [l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hd = lines[0].lstrip("#").split("\t"); cols = [x for x in hd[1:] if x.strip()]
    per = {s: {} for s in cols}
    for l in lines[1:]:
        q = l.split("\t")
        for s, v in zip(cols, q[1:1+len(cols)]):
            v = int(float(v))
            if v: per[s][q[0]] = v
    return cols, per
cols_all, per_all = load_table(f"{O}/table_exp/feature-table.tsv")
cols_f, per_f = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
per_f = {s: {a: v for a, v in per_f[s].items() if a not in drop} for s in per_f}
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
pct = {s: 100*sum(v for a, v in per_all[s].items() if a in host)/max(1, sum(per_all[s].values()))
       for s in cols_all}

if "Figure 2" in changed:
    bl = [s for s in cols_all if meta.get(s) and meta[s]["location"] == "Bloemfontein"]
    nasv = {s: len(per_f.get(s, {})) for s in cols_all}
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for t, m in TISSUE_M.items():
        ss = [s for s in bl if meta[s]["tissue"] == t]
        if not ss: continue
        ax.scatter([pct[s] for s in ss], [nasv.get(s, 0) for s in ss], marker=m, s=52,
                   facecolor="white", edgecolor=CB["black"], linewidth=1.1, label=t, zorder=3)
    _bl = sorted(bl, key=lambda s: (-pct[s], -nasv.get(s, 0)))
    _off = [(5,4),(5,-9),(-16,5),(-16,-9),(5,12),(-16,12)]
    for i, s in enumerate(_bl):
        ax.annotate(s, (pct[s], nasv.get(s, 0)), textcoords="offset points",
                    xytext=(_off[i % len(_off)] if pct[s] > 90 else (5, 4)),
                    fontsize=6.8, color=CB["grey"])
    ax.set_xlabel("Host DNA in sample (%)"); ax.set_ylabel("Fungal ASVs recovered")
    ax.set_title("Figure 2. Bloemfontein: apparent fungal richness tracks host load, not tissue\n"
                 "(non-fungal ASVs removed)", fontsize=10, loc="left")
    ax.legend(title="Tissue", frameon=False, fontsize=8, title_fontsize=8, loc="upper right")
    fig.savefig(f"{FC}/Figure2_bloem_hostload_vs_ASVs.png"); plt.close(fig)
    print(f"wrote {FC}/Figure2_bloem_hostload_vs_ASVs.png")

import zipfile
from collections import Counter
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUT = {"species":97.,"genus":95.,"family":90.,"order":85.,"class":80.,"phylum":80.,"kingdom":80.}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}
SITE_C = {"Bloemfontein":CB["orange"],"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M = {"reference_site":"s","asymptomatic":"o","symptomatic":"^"}
YL = "Unresolved Tremellomycetes\nlineage"   # KMD PROMPT 4, R3a: growth form not established
tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
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
def lab(a):
    if a in yeast: return YL
    rk = ranks_of(tax.get(a, "")); p = pid.get(a, 0.); keep = {}
    for _, r in RANKS:
        if rk[r] and p >= CUT[r]: keep[r] = rk[r]
        else: break
    last = (None, None)
    for _, r in RANKS:
        if keep.get(r): last = (r, keep[r])
    return f"{last[1]} ({last[0]})" if last[0] else "unassigned"
core = [s for s in cols_f if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
gg = [("Christiana","asymptomatic"),("Christiana","symptomatic"),
      ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]

if "Figure 4" in changed:
    groups = gg   # KMD PROMPT 4, R5: Bloemfontein is presence and absence only; core 16 libraries
    comp, ns = {}, {}
    for g in groups:
        ss = [s for s in cols_f if meta.get(s) and (meta[s]["location"], meta[s]["condition_std"]) == g]
        ns[g] = len(ss); c = Counter(); t = 0
        for s in ss:
            for a, v in per_f[s].items(): c[lab(a)] += v; t += v
        comp[g] = {k: 100*v/t for k, v in c.items()} if t else {}
    top = set()
    for g in groups: top |= {k for k, _ in Counter(comp[g]).most_common(7)}
    top = sorted(top, key=lambda k: -sum(comp[g].get(k, 0) for g in groups))
    palette = [CB["vermillion"],CB["blue"],CB["green"],CB["orange"],CB["purple"],CB["skyblue"],
               CB["yellow"],"#7F7F7F","#4D4D4D","#B2DF8A","#FB9A99","#CAB2D6"]
    colmap = {k: palette[i % len(palette)] for i, k in enumerate(top)}
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    for i, g in enumerate(groups):
        bottom = 0
        for k in top:
            v = comp[g].get(k, 0)
            if v <= 0: continue
            ax.bar(i, v, bottom=bottom, width=.62, color=colmap[k], edgecolor="white", linewidth=.4)
            bottom += v
        if 100-bottom > 0.01:
            ax.bar(i, 100-bottom, bottom=bottom, width=.62, color="#DDDDDD", edgecolor="white", linewidth=.4)
    ax.set_xticks(range(len(groups)))
    ax.set_xticklabels([f"{g[0]}\n{g[1]}\n(n={ns[g]})" for g in groups], fontsize=8)
    ax.set_ylabel("Library-level relative abundance\n(% of fungal reads)"); ax.set_ylim(0, 100)
    ax.set_title("Figure 4. Taxonomic composition of the 16 core libraries by site and condition (non-fungal ASVs removed)",
                 fontsize=10, loc="left")
    h = [plt.Rectangle((0,0),1,1,color=colmap[k]) for k in top]+[plt.Rectangle((0,0),1,1,color="#DDDDDD")]
    ax.legend(h, [k.replace("\n"," ") for k in top]+["other"], frameon=False, fontsize=7,
              bbox_to_anchor=(1.01,1), loc="upper left")
    fig.savefig(f"{FC}/Figure4_composition.png"); plt.close(fig)
    print(f"wrote {FC}/Figure4_composition.png")

if "Figure 7" in changed:
    def genus_of(a):
        t = tax.get(a, "")
        if not t or pid.get(a, 0) < 95: return ""
        for part in t.split(";"):
            part = part.strip()
            if part.startswith("g__"):
                v = part[3:].strip()
                return "" if (v.lower() in PH or v.lower().endswith("_incertae_sedis")) else v
        return ""
    curv = {a for a in tax if genus_of(a) == "Curvibasidium"}
    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    for i, g in enumerate(gg):
        ss = [s for s in core if (meta[s]["location"], meta[s]["condition_std"]) == g]
        vals = []
        for s in ss:
            t = sum(per_f[s].values()); c = sum(v for a, v in per_f[s].items() if a in curv)
            vals.append(100*c/t if t else 0)
        ax.scatter([i+(hash(s) % 100-50)/420 for s in ss], vals, s=40, facecolor="white",
                   edgecolor=SITE_C[g[0]], marker=COND_M[g[1]], linewidth=1.3, zorder=3)
        sv = sorted(vals)
        med = sv[len(sv)//2] if len(sv) % 2 else (sv[len(sv)//2-1]+sv[len(sv)//2])/2
        ax.plot([i-.22, i+.22], [med, med], color=CB["black"], linewidth=1.8, zorder=4)
        ax.text(i, -1.6, f"median {med:.2f}%", ha="center", fontsize=7.2)
    ax.set_xticks(range(len(gg)))
    ax.set_xticklabels([f"{g[0]}\n{g[1]}\n(n={sum(1 for s in core if (meta[s]['location'],meta[s]['condition_std'])==g)})"
                        for g in gg], fontsize=8)
    ax.set_ylabel("Curvibasidium (% of fungal reads)"); ax.set_ylim(-3, 22)
    ax.set_title("Figure 7. Per-sample Curvibasidium relative abundance (non-fungal ASVs removed)",
                 fontsize=10, loc="left")
    fig.savefig(f"{FC}/Figure7_curvibasidium.png"); plt.close(fig)
    print(f"wrote {FC}/Figure7_curvibasidium.png")

# Figures 5, 6 and S1 and Table S8 are NOT regenerated here. All three depend on rarefaction
# at the locked depth 57,141, and the correction drops sample 27-C to 57,131. Choosing a new
# depth changes what the manuscript reports, so it is Mark's decision, not this script's.
if "Fig 5/6/S1" in changed:
    newmin = min((sum(per_f[s].values()), s) for s in core)
    print("\nBLOCKED, not regenerated: Figures 5, 6, S1 and Table S8.")
    print(f"  Reason: they rarefy to the locked depth 57,141. After correction the smallest core")
    print(f"  sample is {newmin[1]} at {newmin[0]:,} reads, so 57,141 no longer retains all 16.")
    print("  Needs a decision on the depth before any of them can be redrawn.")

if "Figure S4" in changed:
    gd = list(csv.DictReader(open(f"{C}/guild_assignments_per_ASV.tsv"), delimiter="\t"))
    yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
             if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
    grand = sum(int(r["total_reads"]) for r in gd)
    ft_a = sum(1 for r in gd if r["FungalTraits_primary_lifestyle"])
    ft_r = sum(int(r["total_reads"]) for r in gd if r["FungalTraits_primary_lifestyle"])
    fg_a = sum(1 for r in gd if r["FUNGuild_guild"])
    fg_r = sum(int(r["total_reads"]) for r in gd if r["FUNGuild_guild"])
    fg_c = sum(int(r["total_reads"]) for r in gd if r["FUNGuild_guild"] and r["ASV_ID_md5"] not in yeast)
    vals = [100*ft_a/len(gd), 100*ft_r/grand, 100*fg_a/len(gd), 100*fg_r/grand, 100*fg_c/grand]
    labels = ["FungalTraits\n% ASVs","FungalTraits\n% reads","FUNGuild\n% ASVs",
              "FUNGuild\n% reads\n(as reported)","FUNGuild\n% reads\n(excl. misassigned\nTremellomycetes lineage)"]
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    ax.bar(range(5), vals, color=[CB["green"],CB["green"],CB["purple"],CB["purple"],CB["blue"]], width=.6)
    for i, v in enumerate(vals): ax.text(i, v+1.5, f"{v:.1f}", ha="center", fontsize=8)
    ax.set_xticks(range(5)); ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("Coverage (%)"); ax.set_ylim(0, 100)
    ax.set_title(f"Figure S4. Guild assignment coverage by tool (n = {len(gd):,} fungal ASVs)",
                 fontsize=10, loc="left")
    fig.savefig(f"{FC}/FigureS4_guild_coverage.png"); plt.close(fig)
    print(f"wrote {FC}/FigureS4_guild_coverage.png")
    print("   values: " + ", ".join(f"{l.splitlines()[0]}..={v:.2f}" for l, v in zip(labels, vals)))
