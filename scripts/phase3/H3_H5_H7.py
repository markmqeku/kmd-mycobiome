#!/usr/bin/env python
"""H-3 resolution sensitivity (genus >=95/93/90); H-5 Bloemfontein inventory; H-7 GenBank prep."""
import csv, zipfile, os, datetime
from collections import defaultdict, Counter
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

tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
f = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv",".blast6"))][0]
pid = {}
for line in z.read(f).decode("utf-8","replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v

yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}

lines = [l.rstrip("\n") for l in open(f"{O}/fungi_only/table_exp/feature-table.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [h for h in hdr[1:] if h.strip()]
per = {s: {} for s in cols}; totreads = Counter()
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v; totreads[q[0]] += v
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}

def genus_label(a, thr):
    if a in yeast: return YEAST
    rk = ranks_of(tax.get(a,""))
    if rk["genus"] and pid.get(a,0) >= thr: return rk["genus"]
    for _, r in reversed(RANKS[:5]):        # family..kingdom fallback
        if rk[r] and pid.get(a,0) >= 85: return f"{rk[r]} ({r})"
    return "unassigned"

print("=== H-3: RESOLUTION SENSITIVITY - composition at genus thresholds 95 / 93 / 90 ===")
print("Reported side by side. The threshold is NOT chosen to maximise the story.\n")
core_groups = [("Christiana","asymptomatic"),("Christiana","symptomatic"),
               ("Pretoria","asymptomatic"),("Pretoria","symptomatic")]
res = {}
for thr in (95, 93, 90):
    for g in core_groups:
        ss = [s for s in cols if meta.get(s) and (meta[s]["location"],meta[s]["condition_std"])==g]
        c = Counter(); t = 0
        for s in ss:
            for a, v in per[s].items():
                c[genus_label(a, thr)] += v; t += v
        res[(thr,g)] = (c, t)
for g in core_groups:
    print(f"--- {g[0]} / {g[1]} ---")
    names = set()
    for thr in (95,93,90): names |= {k for k,_ in res[(thr,g)][0].most_common(8)}
    print(f"    {'taxon':46} {'>=95%':>8} {'>=93%':>8} {'>=90%':>8}")
    rowsout = []
    for n in names:
        vals = [100*res[(thr,g)][0].get(n,0)/res[(thr,g)][1] for thr in (95,93,90)]
        rowsout.append((max(vals), n, vals))
    for _, n, vals in sorted(rowsout, reverse=True)[:9]:
        flag = "  THRESHOLD-DEPENDENT" if (max(vals)-min(vals)) > 1.0 else "  robust"
        print(f"    {n:46} {vals[0]:7.1f}% {vals[1]:7.1f}% {vals[2]:7.1f}%{flag}")
    print()

print("=== H-3: do the historical patterns survive? ===")
for taxon in ("Didymella","Curvibasidium","Mycosphaerella"):
    print(f"\n  {taxon}:")
    for g in core_groups:
        vals = []
        for thr in (95,93,90):
            c,t = res[(thr,g)]
            v = sum(x for k,x in c.items() if k.split(" (")[0] == taxon)
            vals.append(100*v/t if t else 0)
        print(f"    {g[0]:11} {g[1]:13} >=95: {vals[0]:5.2f}%  >=93: {vals[1]:5.2f}%  >=90: {vals[2]:5.2f}%")

# ---------------- H-5 ----------------
print("\n\n=== H-5: BLOEMFONTEIN INVENTORY (presence/absence, low-host samples only) ===")
host = {l.strip() for l in open(f"{O}/nontarget_screen/host_asv_ids_ALL.txt") if l.strip()}
full = [l.rstrip("\n") for l in open(f"{O}/table_exp/feature-table.tsv")
        if l.strip() and not l.startswith("# Constructed")]
fh2 = full[0].lstrip("#").split("\t"); fcols = [h for h in fh2[1:] if h.strip()]
fper = {s: {} for s in fcols}
for l in full[1:]:
    q = l.split("\t")
    for s, v in zip(fcols, q[1:1+len(fcols)]):
        v = int(float(v))
        if v: fper[s][q[0]] = v
LOWHOST = 10.0
bl = [s for s in cols if meta.get(s) and meta[s]["location"]=="Bloemfontein"]
sel = []
for s in bl:
    tot = sum(fper[s].values()); hst = sum(v for a,v in fper[s].items() if a in host)
    hp = 100*hst/tot if tot else 0
    if hp <= LOWHOST: sel.append((s, hp, sum(per[s].values())))
print(f"low-host Bloemfontein samples (host <= {LOWHOST}%): {len(sel)}")
for s, hp, fr in sorted(sel, key=lambda x: x[1]):
    print(f"   {s:6} {meta[s]['tissue']:14} {meta[s]['age']:6} host={hp:5.2f}%  fungal_reads={fr:,}")
inv = Counter()
for s, hp, fr in sel:
    for a in per[s]:
        inv[genus_label(a, 95)] += 1
print(f"\npresence/absence taxon inventory across those {len(sel)} samples (no abundances, no diversity metrics):")
print(f"   distinct taxa: {len(inv)}")
for k, v in inv.most_common(25):
    print(f"     {k:50} present in {v} ASVs")
print("   NOTE: inventory only. Not comparable to Christiana/Pretoria (different host load, run and chemistry).")

# ---------------- H-7 ----------------
print("\n\n=== H-7: GenBank submission prep ===")
G = f"{O}/genbank_submission"; os.makedirs(G, exist_ok=True)
seqs, sid, buf = {}, None, []
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line = line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid] = "".join(buf)
        sid, buf = line[1:].split()[0], []
    else: buf.append(line.strip())
if sid: seqs[sid] = "".join(buf)
site_reads = {a: Counter() for a in yeast}
nsamp = Counter()
for s in cols:
    loc = meta[s]["location"]
    for a, v in per[s].items():
        if a in yeast: site_reads[a][loc] += v; nsamp[a] += 1
order = sorted(yeast, key=lambda a: -totreads[a])
with open(f"{G}/Tremellomycetes_ASVs.fasta", "w") as fh:
    for i, a in enumerate(order, 1):
        sr = site_reads[a]
        fh.write(f">KMD_TremASV{i:03d} [organism=Tremellomycetes sp.] [host=Searsia lancea] "
                 f"[country=South Africa] [note=ITS2, primers ITS3/ITS4; md5={a}; "
                 f"reads={totreads[a]}; samples={nsamp[a]}]\n{seqs[a]}\n")
with open(f"{G}/Tremellomycetes_ASV_metadata.tsv", "w", newline="") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["submission_id","ASV_ID_md5","length_bp","total_reads","n_samples",
                "reads_Bloemfontein","reads_Christiana","reads_Pretoria",
                "organism","host","marker","primers","country","identification_note"])
    for i, a in enumerate(order, 1):
        sr = site_reads[a]
        w.writerow([f"KMD_TremASV{i:03d}", a, len(seqs[a]), totreads[a], nsamp[a],
                    sr["Bloemfontein"], sr["Christiana"], sr["Pretoria"],
                    "Tremellomycetes sp.", "Searsia lancea", "ITS2", "ITS3/ITS4", "South Africa",
                    "Lineage absent from UNITE 10.0; misassigned by it to Agaricales (Homophron). "
                    "Nearest named Tremellomycetes reference 85.5% over 276 bp; nearest NCBI match "
                    "MK019123.1 Cryptococcus sp. isolate OTU796 (unnamed environmental sequence) "
                    "98-99% over 363 bp. Genus not assigned; LSU D1/D2 required."])
print(f"   wrote {G}/Tremellomycetes_ASVs.fasta  ({len(order)} sequences)")
print(f"   wrote {G}/Tremellomycetes_ASV_metadata.tsv")
print("   organism field = 'Tremellomycetes sp.' - no genus assigned. NOTHING SUBMITTED.")
