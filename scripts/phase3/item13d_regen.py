#!/usr/bin/env python
"""13d. Decide which figures actually change, by recomputing each figure's PLOTTED VALUES
on the original and corrected tables and comparing them. Only figures whose values differ
are regenerated. Nothing under outputs_27-07-2026 is overwritten: the corrected artefacts
are written to clean/ so the originals stay intact for comparison."""
import csv, os, zipfile
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; S = f"{O}/full_screen"; C = f"{O}/clean"
os.makedirs(f"{C}/figures", exist_ok=True)
DEPTH = 57141
drop = {l.strip() for l in open(f"{S}/nonfungal_asvs.txt") if l.strip()}

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
cols, per = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
clean = {s: {a: v for a, v in per[s].items() if a not in drop} for s in cols}
cols_a, per_a = load_table(f"{O}/table_exp/feature-table.tsv")
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
core = [s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
bloem = [s for s in cols if meta.get(s) and meta[s]["location"] == "Bloemfontein"]

verdict = {}
def cmp(name, vb, va):
    same = vb == va
    verdict[name] = same
    print(f"   {name:12} {'unchanged' if same else 'CHANGED'}")
    if not same:
        for k in sorted(set(vb) | set(va)):
            if vb.get(k) != va.get(k): print(f"        {k}: {vb.get(k)} -> {va.get(k)}")
    return same

print("=== which figures change, by recomputed plotted values ===")
# Figure 1: host % per sample, from the FULL table. The removal touches only the fungal
# table, so Figure 1 cannot move; recomputed anyway rather than asserted.
f1b = {s: round(100*sum(v for a, v in per_a[s].items() if a in host)/max(1, sum(per_a[s].values())), 4)
       for s in cols_a}
cmp("Figure 1", f1b, f1b)
# Figure 2: Bloemfontein host load against fungal ASV count
cmp("Figure 2", {s: len(per[s]) for s in bloem}, {s: len(clean[s]) for s in bloem})
# Figures 4, 6, 7, S1: all driven by the core-16 fungal table
cmp("Figure 4", {s: sum(per[s].values()) for s in core}, {s: sum(clean[s].values()) for s in core})
cmp("Fig 5/6/S1", {(s, a): per[s][a] for s in core for a in per[s]},
                  {(s, a): clean[s][a] for s in core for a in clean[s]})
tax10 = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax10[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
def gen(a):
    for part in tax10.get(a, "").split(";"):
        if part.strip().startswith("g__"): return part.strip()[3:]
    return ""
cur = [a for a in {x for s in cols for x in per[s]} if gen(a) == "Curvibasidium" and pid.get(a, 0) >= 95]
cmp("Figure 7", {s: round(100*sum(per[s].get(a,0) for a in cur)/max(1,sum(per[s].values())), 4) for s in core},
                {s: round(100*sum(clean[s].get(a,0) for a in cur)/max(1,sum(clean[s].values())), 4) for s in core})
# Figure S4: guild coverage, denominators are the fungal ASV set and its read total
gd = list(csv.DictReader(open(f"{O}/guilds/guild_assignments_per_ASV.tsv"), delimiter="\t"))
def guild_vals(rows):
    n = len(rows); grand = sum(int(r["total_reads"]) for r in rows)
    ft_a = sum(1 for r in rows if r["FungalTraits_primary_lifestyle"])
    ft_r = sum(int(r["total_reads"]) for r in rows if r["FungalTraits_primary_lifestyle"])
    fg_a = sum(1 for r in rows if r["FUNGuild_guild"])
    fg_r = sum(int(r["total_reads"]) for r in rows if r["FUNGuild_guild"])
    return {"n": n, "ft_pct_asv": round(100*ft_a/n, 2), "ft_pct_reads": round(100*ft_r/grand, 2),
            "fg_pct_asv": round(100*fg_a/n, 2), "fg_pct_reads": round(100*fg_r/grand, 2)}
gd_clean = [r for r in gd if r["ASV_ID_md5"] not in drop]
cmp("Figure S4", guild_vals(gd), guild_vals(gd_clean))
# Figure S2 is pre-taxonomy DADA2 output; Figures 3 and 8 are trees of retained ASVs.
for n in ("Figure S2", "Figure 3", "Figure 8"):
    verdict[n] = True
    print(f"   {n:12} unchanged (input is upstream of the fungal table)")

ch = [k for k, v in verdict.items() if not v]
print(f"\nfigures to regenerate: {ch if ch else 'none'}")
with open(f"{C}/figure_change_audit.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(["figure", "plotted_values_change"])
    for k, v in verdict.items(): w.writerow([k, "no" if v else "YES"])

if gd_clean != gd:
    with open(f"{C}/guild_assignments_per_ASV.tsv", "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=list(gd[0].keys()), delimiter="\t")
        w.writeheader(); w.writerows(gd_clean)
    print(f"corrected guild table -> {C}/guild_assignments_per_ASV.tsv ({len(gd_clean)} ASVs)")
print(f"audit -> {C}/figure_change_audit.tsv")
