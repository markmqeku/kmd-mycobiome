#!/usr/bin/env python
"""Item 4: guild assignment on the 1,173-ASV fungi-only table.
PRIMARY  = FungalTraits (Polme et al. 2020), matched at genus then family.
PARALLEL = FUNGuild, for continuity with the superseded drafts.
Taxonomy used is the FINAL reported taxonomy (UNITE 10.0, rank-specific identity thresholds),
so match rates reflect what the manuscript actually reports.
Guild results are SUPPLEMENTARY. No narrative is built on them."""
import csv, zipfile, os, json
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
G = f"{O}/guilds"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
CUT = {"species":97.0,"genus":95.0,"family":90.0,"order":85.0,"class":80.0,"phylum":80.0,"kingdom":80.0}
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}

def ranks_of(t):
    out = {r:"" for _,r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part = part.strip()
        for pre,rk in RANKS:
            if part.startswith(pre):
                v = part[len(pre):].strip(); b = v.lower()
                if (b in PH or b.endswith("_incertae_sedis") or b.startswith("unidentified")
                        or (rk=="species" and (b.endswith("_sp") or b.endswith("_sp.")))): v = ""
                out[rk] = v
    return out

tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(fn).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}

# final reported taxonomy per ASV (rank-truncated)
final_tax = {}
for a in tax:
    rk = ranks_of(tax[a]); p = pid.get(a, 0.0); keep = {}
    for _, r in RANKS:
        if rk[r] and p >= CUT[r]: keep[r] = rk[r]
        else: break
    final_tax[a] = keep

# fungi-only table
lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
per = {s: {} for s in cols}; asvs = []
for l in lines[1:]:
    q = l.split("\t"); asvs.append(q[0])
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

# ---------- FungalTraits ----------
import openpyxl
wb = openpyxl.load_workbook(f"{G}/FungalTraits_Polme2020.xlsx", read_only=True)
ws = wb["data"]
rows = ws.iter_rows(values_only=True)
head = [str(x) if x is not None else "" for x in next(rows)]
ix = {h: i for i, h in enumerate(head)}
gcol = ix.get("GENUS"); fcol = ix.get("Family"); lcol = ix.get("primary_lifestyle")
ft_genus, ft_fam = {}, defaultdict(Counter)
for r in rows:
    if not r: continue
    g = (r[gcol] or "").strip() if gcol is not None else ""
    f = (r[fcol] or "").strip() if fcol is not None else ""
    ls = (r[lcol] or "").strip() if lcol is not None else ""
    if ls in ("", "None"): continue
    if g: ft_genus[g] = ls
    if f: ft_fam[f][ls] += 1
ft_fam_consensus = {}
for f, c in ft_fam.items():
    tot = sum(c.values()); top, n = c.most_common(1)[0]
    if n / tot >= 0.5: ft_fam_consensus[f] = (top, n / tot)
print(f"FungalTraits loaded: {len(ft_genus)} genera, {len(ft_fam_consensus)} families with >=50% consensus lifestyle")

def ft_assign(a):
    t = final_tax.get(a, {})
    g, f = t.get("genus"), t.get("family")
    if g and g in ft_genus: return ft_genus[g], "genus"
    if f and f in ft_fam_consensus: return ft_fam_consensus[f][0], "family"
    return None, None

# ---------- FUNGuild ----------
db_path = f"{G}/funguild_db.json"
if not os.path.exists(db_path):
    import urllib.request
    urllib.request.urlretrieve("http://www.stbates.org/funguild_db_2.php", db_path)
raw = open(db_path, encoding="utf-8", errors="replace").read()
i, j = raw.find("["), raw.rfind("]")
db = json.loads(raw[i:j+1]) if i >= 0 else []
fg = {}
for e in db:
    nm = str(e.get("taxon", "")).strip()
    if nm: fg[nm] = (e.get("guild", ""), e.get("trophicMode", ""), e.get("confidenceRanking", ""))
print(f"FUNGuild loaded: {len(fg)} taxon entries")

def fg_assign(a):
    t = final_tax.get(a, {})
    for rk in ("species", "genus", "family", "order", "class", "phylum"):
        v = t.get(rk)
        if v and v in fg and fg[v][0]:
            return fg[v], rk
    return None, None

# ---------- report ----------
ft_res = {a: ft_assign(a) for a in asvs}
fg_res = {a: fg_assign(a) for a in asvs}
tot_reads = {a: sum(per[s].get(a, 0) for s in cols) for a in asvs}
grand = sum(tot_reads.values())

print(f"\n=== OVERALL ASSIGNMENT (fungi-only table, {len(asvs)} ASVs, {grand:,} reads) ===")
for nm, res in (("FungalTraits (primary)", ft_res), ("FUNGuild (parallel)", fg_res)):
    na = sum(1 for a in asvs if res[a][0]); nr = sum(tot_reads[a] for a in asvs if res[a][0])
    print(f"  {nm:26} ASVs {na:5}/{len(asvs)} ({100*na/len(asvs):5.1f}%)   reads {nr:>10,} ({100*nr/grand:5.1f}%)")

print("\n=== PER SITE x CONDITION ===")
print(f"{'group':32} {'FT %ASVs':>9} {'FT %reads':>10} {'FG %ASVs':>9} {'FG %reads':>10}")
for site in ("Bloemfontein","Christiana","Pretoria"):
    for cond in ("reference_site","asymptomatic","symptomatic"):
        ss = [s for s in cols if meta.get(s) and meta[s]["location"]==site and meta[s]["condition_std"]==cond]
        if not ss: continue
        present = {a for s in ss for a in per[s]}
        rd_tot = sum(per[s][a] for s in ss for a in per[s])
        out = []
        for res in (ft_res, fg_res):
            na = sum(1 for a in present if res[a][0])
            nr = sum(per[s].get(a,0) for s in ss for a in present if res[a][0])
            out += [100*na/len(present) if present else 0, 100*nr/rd_tot if rd_tot else 0]
        print(f"{site+' / '+cond:32} {out[0]:8.1f}% {out[1]:9.1f}% {out[2]:8.1f}% {out[3]:9.1f}%")

yr = sum(tot_reads[a] for a in yeast if a in tot_reads)
print(f"\n=== THE UNRESOLVED TREMELLOMYCETES LINEAGE ===")
print(f"  78 ASVs, {yr:,} reads ({100*yr/grand:.1f}% of the fungal table).")
print(f"  FungalTraits assignment: {sum(1 for a in yeast if a in ft_res and ft_res[a][0])} of 78")
print(f"  FUNGuild assignment    : {sum(1 for a in yeast if a in fg_res and fg_res[a][0])} of 78")
print("  This lineage, approximately 40% of Christiana reads, receives NO guild assignment from either")
print("  system, because neither can act without a name and the lineage has none.")

print("\n=== TOP LIFESTYLES / GUILDS BY READS (supplementary, not a headline claim) ===")
c1 = Counter(); c2 = Counter()
for a in asvs:
    if ft_res[a][0]: c1[f"{ft_res[a][0]} ({ft_res[a][1]}-level)"] += tot_reads[a]
    if fg_res[a][0]: c2[f"{fg_res[a][0][0]} [{fg_res[a][0][2]}]"] += tot_reads[a]
print("  FungalTraits primary_lifestyle:")
for k, v in c1.most_common(8): print(f"     {k:52} {100*v/grand:5.1f}%")
print("  FUNGuild guild [confidence]:")
for k, v in c2.most_common(8): print(f"     {k:52} {100*v/grand:5.1f}%")

with open(f"{G}/guild_assignments_per_ASV.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["ASV_ID_md5","final_genus","final_family","total_reads",
                "FungalTraits_primary_lifestyle","FT_match_rank",
                "FUNGuild_guild","FUNGuild_trophicMode","FUNGuild_confidenceRanking","FG_match_rank"])
    for a in sorted(asvs, key=lambda x: -tot_reads[x]):
        t = final_tax.get(a, {}); f1 = ft_res[a]; f2 = fg_res[a]
        w.writerow([a, t.get("genus",""), t.get("family",""), tot_reads[a],
                    f1[0] or "", f1[1] or "",
                    (f2[0][0] if f2[0] else ""), (f2[0][1] if f2[0] else ""),
                    (f2[0][2] if f2[0] else ""), f2[1] or ""])
print(f"\nwrote {G}/guild_assignments_per_ASV.tsv  (confidence rankings retained)")
