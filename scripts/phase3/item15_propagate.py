#!/usr/bin/env python
"""Item 15. Propagate the 13b correction into every derived artefact, writing to clean/ so
nothing under outputs_27-07-2026 is overwritten. Rebuilds the representative sequences, the
ASV catalogue (Table S1), the Bloemfontein inventory (Table S2) and the threshold table
(Table S3), then writes a manifest and lists the prose changes the correction forces.
Prose is LISTED, never edited: every one of those changes what a manuscript claims."""
import csv, os, subprocess, sys, zipfile
from collections import Counter, defaultdict
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"; C = f"{O}/clean"; P3 = f"{B}/phase3"
os.makedirs(C, exist_ok=True)
drop = {l.strip() for l in open(f"{O}/full_screen/nonfungal_asvs.txt") if l.strip()}
print(f"removing {len(drop)} confirmed non-fungal ASVs\n")

# ---- 1. filtered representative sequences ----
seqs, sid, buf = {}, None, []
for line in open(f"{O}/fungi_only/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line)
if sid: seqs[sid] = "".join(buf)
with open(f"{C}/rep_fungi_clean.fasta", "w") as fh:
    for a, s in seqs.items():
        if a not in drop: fh.write(f">{a}\n{s}\n")
kept = len(seqs) - len(drop & set(seqs))
print(f"1. representative sequences: {len(seqs)} -> {kept}")

# ---- 2. ASV catalogue (Table S1), rebuilt by the same builder, on the corrected inputs ----
cat_out = f"{C}/asv_catalogue"
r = subprocess.run([sys.executable, f"{P3}/build_asv_catalogue.py",
                    f"{O}/final_taxonomy.tsv", f"{C}/rep_fungi_clean.fasta",
                    f"{C}/feature-table-clean.tsv", "NONE", cat_out, "UNITE 10.0 (rank thresholds)"],
                   capture_output=True, text=True)
print(f"2. ASV catalogue -> {cat_out}  ({'ok' if r.returncode == 0 else 'FAILED'})")
if r.returncode != 0: print("   " + (r.stderr or r.stdout)[-500:])

# ---- 3 and 4. Tables S2 and S3 on the corrected table ----
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
def load_tax(p):
    d = {}
    for r_ in csv.reader(open(p), delimiter="\t"):
        if r_ and r_[0] != "Feature ID" and not r_[0].startswith("#"): d[r_[0]] = r_[1]
    return d
final = load_tax(f"{O}/final_taxonomy.tsv")
yeast = {r_[0] for r_ in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r_ and r_[0] != "ASV_ID_md5" and r_[2] == "1"}
def label(a):
    if a in yeast: return "Tremellomycetes yeast lineage (ITS2-unresolved)"
    rk = ranks_of(final.get(a, "")); last = ("none", "")
    for _, r_ in RANKS:
        if rk[r_]: last = (r_, rk[r_])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"
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
cols_c, per_c = load_table(f"{C}/feature-table-clean.tsv")
cols_a, per_a = load_table(f"{O}/table_exp/feature-table.tsv")
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c_: i for i, c_ in enumerate(h)}
    for r_ in rd:
        if r_ and r_[0].strip(): meta[r_[0]] = {c_: r_[di[c_]] for c_ in h}
sel = []
for s in [x for x in cols_a if meta.get(x) and meta[x]["location"] == "Bloemfontein"]:
    tot = sum(per_a[s].values()); hs = sum(v for a, v in per_a[s].items() if a in host)
    hp = 100*hs/tot if tot else 0
    if hp <= 10: sel.append((s, hp, sum(per_c.get(s, {}).values())))
sel.sort(key=lambda x: x[1])
occ, dis = Counter(), defaultdict(set)
for s, _, _ in sel:
    for a in per_c.get(s, {}): occ[label(a)] += 1; dis[label(a)].add(a)
with open(f"{C}/TableS2_bloemfontein_inventory.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["# Bloemfontein presence-absence inventory. Low-host samples only (host <=10%)."])
    w.writerow(["# Confirmed non-fungal ASVs removed. No abundances, no diversity metrics, "
                "no cross-site comparison."])
    w.writerow(["# n_ASVs counts DISTINCT ASVs; n_occurrences counts (sample, ASV) pairs."])
    w.writerow(["sample","tissue","age","host_pct","fungal_reads"])
    for s, hp, fr in sel: w.writerow([s, meta[s]["tissue"], meta[s]["age"], f"{hp:.2f}", fr])
    w.writerow([]); w.writerow(["taxon_present","n_ASVs","n_occurrences"])
    for k, v in sorted(dis.items(), key=lambda kv: -len(kv[1])):
        w.writerow([k, len(v), occ[k]])
print(f"3. Table S2 -> {len(sel)} samples, {len(dis)} taxa, "
      f"{sum(len(v) for v in dis.values())} distinct ASVs")

t82 = load_tax(f"{B}/phase2/taxonomy/exp/taxonomy.tsv")
t10 = load_tax(f"{O}/taxonomy_unite10/exp/taxonomy.tsv")
fungi = {a for s in cols_c for a in per_c[s]}
with open(f"{C}/TableS3_taxonomy_thresholds.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow([f"# ASVs assigned (%) at each rank, corrected fungi-only table (n={len(fungi)})"])
    w.writerow(["rank","UNITE_8.2_permissive","UNITE_10.0_permissive","UNITE_10.0_rank_threshold"])
    for _, rk in RANKS:
        a_ = sum(1 for x in fungi if ranks_of(t82.get(x, ""))[rk])
        b_ = sum(1 for x in fungi if ranks_of(t10.get(x, ""))[rk])
        c_ = sum(1 for x in fungi if ranks_of(final.get(x, ""))[rk])
        w.writerow([rk, f"{100*a_/len(fungi):.1f}", f"{100*b_/len(fungi):.1f}", f"{100*c_/len(fungi):.1f}"])
print(f"4. Table S3 -> denominator now n={len(fungi)}")

# ---- 5. GenBank submission file: one of the 78 is not fungal ----
import re as _re
gb_in = f"{O}/genbank_submission/Tremellomycetes_ASVs.fasta"
kept_gb, dropped_gb, keep = [], [], True
for line in open(gb_in):
    if line.startswith(">"):
        m = _re.search(r"md5=([0-9a-f]+)", line)
        keep = not (m and m.group(1) in drop)
        (kept_gb if keep else dropped_gb).append(line.strip())
    if keep: pass
lines = open(gb_in).read().splitlines()
out_l, keep = [], True
for line in lines:
    if line.startswith(">"):
        m = _re.search(r"md5=([0-9a-f]+)", line)
        keep = not (m and m.group(1) in drop)
    if keep: out_l.append(line)
os.makedirs(f"{C}/genbank_submission", exist_ok=True)
open(f"{C}/genbank_submission/Tremellomycetes_ASVs.fasta", "w").write("\n".join(out_l) + "\n")
print(f"\n5. GenBank submission: {len(kept_gb)+len(dropped_gb)} -> {len(kept_gb)} sequences")
for h in dropped_gb:
    sid = h[1:].split()[0]; md5 = _re.search(r"md5=([0-9a-f]+)", h).group(1)
    print(f"   REMOVED {sid} (md5 {md5[:12]}), labelled 'Tremellomycetes sp.' but non-fungal")
print("   NOTHING WAS EVER SUBMITTED, so no correction to GenBank is required.")

# ---- 6. prose changes forced by the correction: LISTED, NOT APPLIED ----
gd = list(csv.DictReader(open(f"{C}/guild_assignments_per_ASV.tsv"), delimiter="\t"))
grand = sum(int(x["total_reads"]) for x in gd)
ft_a = sum(1 for x in gd if x["FungalTraits_primary_lifestyle"])
ft_r = sum(int(x["total_reads"]) for x in gd if x["FungalTraits_primary_lifestyle"])
fg_a = sum(1 for x in gd if x["FUNGuild_guild"])
_core = [s for s in cols_c if meta.get(s) and meta[s]["location"] in ("Christiana", "Pretoria")]
_mn = min((sum(per_c[s].values()), s) for s in _core)
todo = [
 ("FINAL_NUMBERS 10", "Distinct taxa recorded", "192",
  f"{len(dis)} (final taxonomy basis; see item 14)"),
 ("FINAL_NUMBERS 10", "Most frequent row", "counts labelled 'ASVs' are occurrences",
  "report distinct ASVs and occurrences separately"),
 ("FINAL_NUMBERS 10", "Reason row, P4/P5 ASV counts", "329/326",
  f"{len(per_c.get('P4',{}))}/{len(per_c.get('P5',{}))}"),
 ("DRAFT_RESULTS 3.8", "taxa in the six low-host samples", "192", str(len(dis))),
 ("DRAFT_RESULTS 3.9", "FungalTraits ASVs", "542 of 1,173 (46.2 percent)",
  f"{ft_a} of {len(gd)} ({100*ft_a/len(gd):.1f} percent)"),
 ("DRAFT_RESULTS 3.9", "FungalTraits read coverage", "23.7 percent", f"{100*ft_r/grand:.1f} percent"),
 ("DRAFT_RESULTS 3.9", "FUNGuild ASVs", "692 ASVs (59.0 percent)",
  f"{fg_a} ASVs ({100*fg_a/len(gd):.1f} percent)"),
 ("DRAFT_RESULTS 3.10", "final paragraph", "'remains in the fungal table ... should be removed'",
  "it has been removed, along with the other confirmed non-fungal ASVs"),
 ("DRAFT_METHODS 2.4", "non-target screening", "host screen only",
  "add the whole-table NCBI nt screen, its criterion and its two failure modes"),
 ("DRAFT_RESULTS 3.3 / 3.9", "size of the unresolved yeast lineage", "78 ASVs",
  f"{len(kept_gb)} ASVs; {len(dropped_gb)} were green algae"),
 ("FINAL_NUMBERS 11", "GenBank material", "78 sequences",
  f"{len(kept_gb)} sequences; {len(dropped_gb)} withdrawn "
  f"({', '.join(h[1:].split()[0] for h in dropped_gb)}), never submitted"),
 ("Figure 3", "placement tree of the yeast lineage", "78 tips",
  f"{len(kept_gb)} tips once the {len(dropped_gb)} algal ASVs are dropped"),
 ("SESSION_STATE", "locked rarefaction depth 57,141", "16/16 retained",
  f"smallest core sample {_mn[1]} now holds {_mn[0]:,}; DECISION REQUIRED, nothing changed"),
]
with open(f"{C}/prose_changes_required.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["location","item","current","corrected"]); w.writerows(todo)
print(f"\n5. prose changes required (NOT applied), {len(todo)} items:")
for loc, it, cur, new in todo:
    print(f"   {loc:20} {it[:38]:38} {cur[:30]:30} -> {new[:44]}")
print(f"\n-> {C}/prose_changes_required.tsv")
