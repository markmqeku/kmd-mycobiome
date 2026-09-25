#!/usr/bin/env python
"""KMD PROMPT 4, R6d. Supplementary Figure S5: study design and workflow, from sampling at three sites to
community analysis. Every value is taken from ANSWERS.md revision 4 (sampling design), the Methods and
FINAL_NUMBERS; nothing is computed here. Line art, drawn with matplotlib."""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
F = f"{O}/figures"; os.makedirs(F, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "ps.fonttype": 42})
BLUE, GREEN, ORANGE, GREY, BLACK = "#0072B2", "#009E73", "#E69F00", "#F2F2F2", "#000000"

fig, ax = plt.subplots(figsize=(6.85, 7.4))
ax.set_xlim(0, 100); ax.axis("off")


def box(x, y, w, h, title, body, edge=BLACK, face=GREY):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                linewidth=0.8, edgecolor=edge, facecolor=face))
    ax.text(x + w / 2, y + h - 1.2, title, ha="center", va="top", fontsize=6.6, fontweight="bold")
    ax.text(x + w / 2, y + h - 3.9, body, ha="center", va="top", fontsize=5.6, linespacing=1.35)


def arrow(x1, y1, x2, y2):
    ax.annotate("", xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle="-|>", lw=0.8, color=BLACK, shrinkA=0, shrinkB=0))


SITE = [(ORANGE, "#FDF3DC"), (BLUE, "#E1EEF7"), (GREEN, "#DFF3EC")]
PLAIN = (BLACK, GREY)
POOL = (BLACK, "#EEF5F2")
# rows top to bottom: (height, [(x, width, title, body, edge, face), ...])
ROWS = [
    (13, [(1, 31, "Bloemfontein", "campus grounds\nreference site with no\nrecorded KMD occurrence\n26 April 2016\n8 trees", *SITE[0]),
          (34.5, 31, "Christiana", "municipal street trees\nKMD-affected site\n10 May 2017\n24 symptomatic trees\n12 asymptomatic trees", *SITE[1]),
          (68, 31, "Pretoria", "farm\nKMD-affected site\n18 May 2017\n24 symptomatic trees\n12 asymptomatic trees", *SITE[2])]),
    (11, [(1, 31, "One tree per library", "each tree supplied several\nlibraries (tissue types,\ntwo age classes)\n32 libraries", *SITE[0]),
          (34.5, 64.5, "Four trees pooled per library", "tissue from four trees of the same site and condition pooled\n"
           "Christiana: 6 symptomatic and 3 asymptomatic libraries\nPretoria: 6 symptomatic and 3 asymptomatic libraries", *POOL)]),
    (8, [(1, 98, "Surface sterilisation, lyophilisation and DNA extraction",
          "sodium hypochlorite, distilled water, 96 percent ethanol, distilled water; lyophilised; TissueLyser II\n"
          "NucleoSpin Plant II Midi, 40 mg lyophilised tissue; NanoDrop Lite quantification", *PLAIN)]),
    (8, [(1, 98, "ITS2 amplification and library preparation",
          "ITS3 and ITS4 with Nextera overhang adapters, 25 cycles; AMPure XP clean-up;\n"
          "Nextera XT indexing, 8 cycles; libraries normalised to 4 nM", *PLAIN)]),
    (8, [(1, 48, "Library check", "Bloemfontein only: Agilent 2100 Bioanalyzer,\nHigh Sensitivity DNA assay, 15 and 16 August 2016", *SITE[0]),
         (51, 48, "Library check", "Christiana and Pretoria:\nno records exist", *POOL)]),
    (8, [(1, 48, "Illumina MiSeq run 3 (2016)", "2 x 250 bp paired-end reads\nBloemfontein", *SITE[0]),
         (51, 48, "Illumina MiSeq run 13 (2017)", "2 x 301 bp paired-end reads\nChristiana and Pretoria", *POOL)]),
    (8, [(1, 98, "Primer removal and run-specific DADA2 denoising",
          "cutadapt; truncation 228 and 200 (run 3), 280 and 230 (run 13); tables merged\n"
          "45 usable libraries (Bloemfontein 29, Christiana 9, Pretoria 7); 4,010,617 reads; 1,322 ASVs", *PLAIN)]),
    (8, [(1, 98, "Host and non-fungal screening",
          "149 host ASVs removed (1,487,006 reads); 31 non-fungal ASVs removed (UNITE 10.0 alignment, NCBI nt)\n"
          "1,142 fungal ASVs, 2,517,177 reads", *PLAIN)]),
    (8, [(1, 31, "Taxonomy", "VSEARCH consensus, UNITE 10.0,\nrank-specific identity thresholds", *PLAIN),
         (34.5, 31, "Phylogenetic placement", "unresolved lineage, 71 ASVs;\nMAFFT and IQ-TREE", *PLAIN),
         (68, 31, "Community analysis", "16 core libraries, 57,124 reads;\nBloemfontein presence and absence", *PLAIN)]),
]
GAP = 4.0
total = sum(h for h, _ in ROWS) + GAP * (len(ROWS) - 1)
ax.set_ylim(-1, total + 1)
y = total
prev = None
for h, boxes in ROWS:
    y0 = y - h
    for (x, w, ti, bo, ed, fa) in boxes:
        box(x, y0, w, h, ti, bo, edge=ed, face=fa)
    if prev is not None:
        pboxes, py0 = prev
        if len(boxes) >= len(pboxes):      # one to one or fan out: an arrow into each new box
            xs = [x + w / 2 for (x, w, *_r) in boxes]
        else:                              # fan in: an arrow out of each previous box
            xs = [x + w / 2 for (x, w, *_r) in pboxes]
        for c in xs:
            arrow(c, py0 - 0.6, c, y + 0.6)
    prev = (boxes, y0)
    y = y0 - GAP
fig.subplots_adjust(left=0.01, right=0.99, top=0.995, bottom=0.005)
fig.savefig(f"{F}/FigureS5_workflow.png", dpi=300)
plt.close(fig)
print("wrote FigureS5_workflow.png")
