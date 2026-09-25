#!/usr/bin/env python
"""Item 14. Reconcile the frozen 192-taxon Bloemfontein inventory (H-5) against the
199-taxon rebuild (Table S2). Recompute BOTH on one run so the sample set is provably
identical and the only remaining variable is the labelling basis."""
import csv, zipfile
from collections import Counter
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
RANKS = [("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),
         ("f__","family"),("g__","genus"),("s__","species")]
PH = {"unidentified","unclassified","incertae_sedis","unknown",""}
YEAST = "Tremellomycetes yeast lineage (ITS2-unresolved)"

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

def load_tax(p):
    d = {}
    for r in csv.reader(open(p), delimiter="\t"):
        if r and r[0] != "Feature ID" and not r[0].startswith("#"): d[r[0]] = r[1]
    return d

tax10 = load_tax(f"{O}/taxonomy_unite10/exp/taxonomy.tsv")
final = load_tax(f"{O}/final_taxonomy.tsv")
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

def load_table(p):
    lines = [l.rstrip("\n") for l in open(p) if l.strip() and not l.startswith("# Constructed")]
    hd = lines[0].lstrip("#").split("\t"); cols = [x for x in hd[1:] if x.strip()]
    per = {s:{} for s in cols}
    for l in lines[1:]:
        q = l.split("\t")
        for s, v in zip(cols, q[1:1+len(cols)]):
            v = int(float(v))
            if v: per[s][q[0]] = v
    return cols, per
cols_f, per_f = load_table(f"{O}/fungi_only/table_exp/feature-table.tsv")
cols_a, per_a = load_table(f"{O}/table_exp/feature-table.tsv")
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

# ---- A. sample set, computed once, used by both labellings ----
bl = [s for s in cols_a if meta.get(s) and meta[s]["location"] == "Bloemfontein"]
sel = []
for s in bl:
    tot = sum(per_a[s].values()); hs = sum(v for a,v in per_a[s].items() if a in host)
    hp = 100*hs/tot if tot else 0
    if hp <= 10.0: sel.append((s, hp))
sel.sort(key=lambda x: x[1])
FROZEN = ["P5","P32","P4","P18","P19","P25"]
got = [s for s,_ in sel]
print("=== A. SAMPLE SET ===")
print(f"Bloemfontein samples in table: {len(bl)}")
print(f"low-host (<=10%) selected:     {len(sel)}  -> {got}")
print(f"FINAL_NUMBERS section 10 lists: {FROZEN}")
same = set(got) == set(FROZEN)
print(f"SAMPLE SETS IDENTICAL: {same}")
for s, hp in sel: print(f"   {s:5} host={hp:5.2f}%  fungal_ASVs={len(per_f.get(s,{})):>4}")
if not same:
    print("\nSAMPLE SETS DIFFER. Not updating FINAL_NUMBERS. Stopping here per item 14.")
    raise SystemExit(1)

# ---- B. the two labellings on that identical sample set ----
def label_frozen(a):                      # H-5 basis: UNITE 10.0 raw + flat 95/85 cuts
    if a in yeast: return YEAST
    rk = ranks_of(tax10.get(a,""))
    if rk["genus"] and pid.get(a,0) >= 95: return rk["genus"]
    for _, r in reversed(RANKS[:5]):
        if rk[r] and pid.get(a,0) >= 85: return f"{rk[r]} ({r})"
    return "unassigned"

def label_final(a):                       # Table S2 basis: final rank-threshold taxonomy
    if a in yeast: return YEAST
    rk = ranks_of(final.get(a,"")); last = ("none","")
    for _, r in RANKS:
        if rk[r]: last = (r, rk[r])
    return f"{last[1]} ({last[0]})" if last[1] else "unassigned"

asvs = sorted({a for s,_ in sel for a in per_f.get(s, {})})
fro, fin = Counter(), Counter()
for a in asvs:
    fro[label_frozen(a)] += 1; fin[label_final(a)] += 1
print(f"\n=== B. LABELLING ===")
print(f"ASVs present across the 6 samples: {len(asvs)}")
print(f"distinct taxa, frozen basis  (UNITE 10.0 raw, flat 95/85): {len(fro)}")
print(f"distinct taxa, final basis   (final rank-threshold taxonomy): {len(fin)}")

only_fro = sorted(set(fro) - set(fin)); only_fin = sorted(set(fin) - set(fro))
print(f"\nlabels only under the frozen basis ({len(only_fro)}):")
for k in only_fro: print(f"   - {k}  ({fro[k]} ASVs)")
print(f"\nlabels only under the final basis ({len(only_fin)}):")
for k in only_fin: print(f"   + {k}  ({fin[k]} ASVs)")

changed = [(a, label_frozen(a), label_final(a)) for a in asvs if label_frozen(a) != label_final(a)]
print(f"\nASVs whose label changes: {len(changed)}")
for a, x, y in changed[:40]: print(f"   {a[:12]}  {x[:44]:44} -> {y}")
if len(changed) > 40: print(f"   ... {len(changed)-40} more")
print(f"\nyeast lineage ASVs in this set: {sum(1 for a in asvs if a in yeast)} "
      f"(identical under both bases, so the yeast relabelling is NOT the cause)")
