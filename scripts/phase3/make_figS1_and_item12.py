#!/usr/bin/env python
"""Figure S1 (rarefaction, core 16, fungi-only table) + item 12 deliverables:
12b figure/table manifest, 12d SRA metadata, 12e INSERT checklist."""
import csv, os, re, glob, statistics as st
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
B = "<KMD_ROOT>/__reanalysis_2026-06"
O = f"{B}/phase3/outputs_27-07-2026"
F = f"{O}/figures"; G = f"{O}/gap_closure"; A = f"{O}/assembly"; os.makedirs(A, exist_ok=True)
CB = {"blue":"#0072B2","green":"#009E73","black":"#000000"}
SITE_C = {"Christiana":CB["blue"],"Pretoria":CB["green"]}
COND_M = {"asymptomatic":"o","symptomatic":"^"}
plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False,
                     "savefig.dpi":300,"savefig.bbox":"tight"})

meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd=csv.reader(fh,delimiter="\t"); h=next(rd); next(rd); di={c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]]={c:r[di[c]] for c in h}

# ---------- Figure S1 ----------
f = [p for p in glob.glob(f"{G}/rarefaction_exp/*.csv") if "observed" in os.path.basename(p).lower()]
if f:
    rows = list(csv.reader(open(f[0])))
    hdr = rows[0]
    cols = {}
    for i,hh in enumerate(hdr):
        m = re.match(r"depth-(\d+)_iter-\d+", hh)
        if m: cols.setdefault(int(m.group(1)), []).append(i)
    depths = sorted(cols)
    fig, ax = plt.subplots(figsize=(6.4,4.2))
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
    ax.axvline(57141, color=CB["black"], linestyle="--", linewidth=1)
    ax.annotate("rarefaction depth\n57,141", (57141, ax.get_ylim()[1]*0.25),
                textcoords="offset points", xytext=(-64,0), fontsize=7.5)
    ax.set_xlabel("Sequencing depth (reads)"); ax.set_ylabel("Observed ASVs")
    ax.set_title(f"Figure S1. Rarefaction curves, {n_plot} core samples (fungi-only table)",
                 fontsize=10, loc="left")
    from matplotlib.lines import Line2D
    leg = [Line2D([0],[0],color=SITE_C[s],lw=1.5,label=s) for s in ("Christiana","Pretoria")]
    leg += [Line2D([0],[0],color=CB["black"],marker=COND_M[c],lw=0,label=c) for c in COND_M]
    ax.legend(handles=leg, frameon=False, fontsize=8)
    fig.savefig(f"{F}/FigureS1_rarefaction_core16.png"); plt.close(fig)
    print(f"wrote FigureS1_rarefaction_core16.png  ({n_plot} samples)")
else:
    print("FigureS1: no observed_features csv found in gap_closure/rarefaction_exp")

# ---------- 12d SRA metadata ----------
RR = "<KMD_ROOT>/__restructured_25-07-2026/raw_reads"
files = {}
for p in glob.glob(f"{RR}/*_R1_001.fastq*"):
    b = os.path.basename(p)
    if b.endswith(".filepart"): continue
    sid = re.sub(r"_S\d+_L001_R1_001\.fastq(\.gz)?$", "", b)
    files[sid] = (b, b.replace("_R1_", "_R2_"))
COORD = {"Bloemfontein":("29.11064 S 26.18103 E","2017"),
         "Christiana":("27.90692 S 25.16358 E","2017"),
         "Pretoria":("25.74725 S 28.25878 E","2017")}
RUN = {"Bloemfontein":("run3","ARBL2","2x250"),
       "Christiana":("run13","BBY9M","2x301"),"Pretoria":("run13","BBY9M","2x301")}
with open(f"{A}/SRA_metadata.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t")
    w.writerow(["sample_name","library_ID","title","library_strategy","library_source",
                "library_selection","library_layout","platform","instrument_model",
                "design_description","filetype","filename","filename2",
                "organism","host","isolation_source","geo_loc_name","lat_lon","collection_date",
                "tissue","condition","sequencing_run","flowcell","read_configuration"])
    n=0
    for sid in sorted(meta):
        if sid not in files: continue
        m=meta[sid]; loc=m["location"]; ll,yr=COORD[loc]; run,fc,rc=RUN[loc]
        w.writerow([sid, f"KMD_{sid}",
            f"ITS2 amplicon metabarcoding of Searsia lancea {m['tissue'].lower()}, {loc}",
            "AMPLICON","GENOMIC","PCR","paired","ILLUMINA","Illumina MiSeq",
            "ITS2 amplified with ITS3/ITS4 (White et al. 1990) carrying Nextera overhangs; "
            "Illumina 16S metagenomic library preparation adapted for fungi",
            "fastq", files[sid][0], files[sid][1],
            "fungal metagenome","Searsia lancea", m["tissue"],
            f"South Africa: {loc}", ll, yr, m["tissue"], m["condition_std"], run, fc, rc])
        n+=1
print(f"wrote {A}/SRA_metadata.tsv  ({n} samples)")
print("  NOTE: lat_lon converted to decimal degrees; collection_date is year only "
      "(month not recorded in any source). NOTHING SUBMITTED.")

# ---------- 12e INSERT checklist ----------
ins=[]
for p in ("<KMD_ROOT>/DRAFT_METHODS.md",
          "<KMD_ROOT>/DRAFT_RESULTS.md"):
    for i,l in enumerate(open(p,encoding="utf-8",errors="replace"),1):
        for m in re.finditer(r"\\?\[INSERT:?([^\]]*)\]", l):
            ins.append((os.path.basename(p), i, m.group(1).strip()[:150]))
with open(f"{A}/INSERT_checklist.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["file","line","outstanding_item"])
    w.writerows(ins)
print(f"\nwrote {A}/INSERT_checklist.tsv  ({len(ins)} outstanding)")
for f_,i,t in ins: print(f"   {f_}:{i}  {t}")
