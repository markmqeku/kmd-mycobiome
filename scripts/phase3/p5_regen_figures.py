#!/usr/bin/env python
"""KMD PROMPT 5, D0c (from the Prompt 4 R6 script). Regenerate every main and supplementary figure from its own source script, in the
format Fungal Ecology requires (TIFF or EPS, 300 dpi or more), into SUBMISSION_PACKAGE/figures. Run in WSL with the qiime2-amplicon-2024.10
environment, which holds the matplotlib and Biopython the scripts were written for.

Each source script runs unchanged. Only Figure.savefig is intercepted. Before saving, clean_figure() removes
the figure-level title line ("Figure N. ...", drawn as an axes title or a suptitle; the caption in the manuscript
carries the title) and changes every "sample" or "samples" in the figure's text to "library" or "libraries",
since every plotted unit is a sequencing library. Panel labels (a., b., q = 0, Bray-Curtis, neighbourhood names)
are kept. Then: a figure the script saves under
a target name is written as a 600 dpi TIFF (LZW) and as a vector PDF with TrueType fonts embedded,
under its figure-order name. Anything else the script saves (superseded versions of other figures) is
discarded, and nothing is written back into the analysis outputs. Resolution comes from re-rendering
the vector figure at 600 dpi, never from upscaling a PNG.

SAJB (COWORK_FINDINGS/JOURNAL_SAJB.md 13): charts with lines and text are "combinations bitmapped
line/half-tone", minimum 500 dpi as TIFF or JPEG, or vector EPS or PDF with fonts embedded."""
import os, runpy, sys
import matplotlib
matplotlib.use("Agg")
matplotlib.rcParams["pdf.fonttype"] = 42          # embed TrueType fonts in the PDFs
matplotlib.rcParams["ps.fonttype"] = 42
from matplotlib.figure import Figure

ROOT = "<KMD_ROOT>"
P3 = f"{ROOT}/__reanalysis_2026-06/phase3"
OUT = f"{ROOT}/SUBMISSION_PACKAGE/figures"
os.makedirs(OUT, exist_ok=True)
assert os.path.realpath(OUT).startswith(os.path.realpath(ROOT) + "/")

# the script that drew each embedded figure, and the file name it saves it under
PLAN = [
    ("make_figures.py",      {"Figure1_host_by_site.png": "Figure1", "FigureS2_dada2_retention.png": "FigureS2",
                              "FigureS3_assignment_by_rank.png": "FigureS3"}),
    ("item13d_make_figs.py", {"Figure2_bloem_hostload_vs_ASVs.png": "Figure2", "Figure4_composition.png": "Figure4",
                              "Figure7_curvibasidium.png": "Figure7", "FigureS4_guild_coverage.png": "FigureS4"}),
    ("make_fig3_71.py",      {"Figure3_tremellomycetes_placement_71.png": "Figure3"}),
    ("item17b_figs.py",      {"Figure5_hill_numbers.png": "Figure5", "Figure6_ordination.png": "Figure6"}),
    ("make_fig8_fig5.py",    {"Figure8_placement_by_neighbourhood.png": "Figure8"}),
    ("item64g_figure.py",    {"Figure9_rare_fraction.png": "Figure9"}),
    ("item17b_figS1.py",     {"FigureS1_rarefaction_core16.png": "FigureS1"}),
    ("make_figS5_workflow.py", {"FigureS5_workflow.png": "FigureS5"}),
]

_orig_savefig = Figure.savefig

import re as _re
EMB = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/figures/docx_embed_p5"
os.makedirs(EMB, exist_ok=True)
cleaned = {}


def _lib(t):
    t2 = _re.sub(r"\bsamples\b", "libraries", t)
    t2 = _re.sub(r"\bSamples\b", "Libraries", t2)
    t2 = _re.sub(r"\bsample\b", "library", t2)
    return _re.sub(r"\bSample\b", "Library", t2)


def clean_figure(fig, target):
    log = []
    if fig._suptitle is not None and fig._suptitle.get_text().startswith("Figure"):
        log.append("suptitle removed: " + fig._suptitle.get_text().split("\n")[0][:60]); fig._suptitle.set_text("")
    for ax in fig.axes:
        for loc in ("left", "center", "right"):
            t = ax.get_title(loc=loc)
            if t.startswith("Figure"):
                log.append("title removed: " + t.split("\n")[0][:60]); ax.set_title("", loc=loc)
        texts = [ax.xaxis.label, ax.yaxis.label, ax.title, ax._left_title, ax._right_title] + list(ax.texts)
        leg = ax.get_legend()
        if leg is not None: texts += list(leg.get_texts()) + [leg.get_title()]
        texts += list(ax.get_xticklabels()) + list(ax.get_yticklabels())
        for tx in texts:
            old = tx.get_text()
            if _re.search(r"\b[Ss]amples?\b", old):
                tx.set_text(_lib(old)); log.append(f"text: {old[:50]!r} -> {_lib(old)[:50]!r}")
    for tx in fig.texts:
        old = tx.get_text()
        if _re.search(r"\b[Ss]amples?\b", old): tx.set_text(_lib(old)); log.append(f"fig text: {old[:50]!r}")
    cleaned[target] = log

written, discarded = [], []
current = {}


def savefig(self, fname, *args, **kwargs):
    base = os.path.basename(str(fname))
    target = current.get(base)
    if target is None:
        discarded.append(base)
        return None
    kw = {k: v for k, v in kwargs.items() if k not in ("dpi", "format", "pil_kwargs")}
    kw.setdefault("bbox_inches", "tight")
    tif = f"{OUT}/{target}.tif"
    clean_figure(self, target)
    _orig_savefig(self, f"{EMB}/{target}.png", dpi=300, format="png", **kw)

    _orig_savefig(self, tif, dpi=600, format="tiff", pil_kwargs={"compression": "tiff_lzw"}, **kw)
    if target == "FigureS5":                       # line art: also EPS
        _orig_savefig(self, f"{OUT}/{target}.eps", format="eps", **kw)
    written.append((target, base))
    return None


Figure.savefig = savefig
os.chdir(P3)
sys.path.insert(0, P3)
for script, targets in PLAN:
    current = targets
    print(f"--- running {script}", flush=True)
    runpy.run_path(f"{P3}/{script}", run_name="__main__")
    matplotlib.pyplot.close("all")

print("\ncleaned:")
for k in sorted(cleaned):
    print("  ", k, cleaned[k])
print("\nwritten:")
for t, b in sorted(written):
    print(f"  {t:9} from {b}")
print("discarded (superseded versions drawn by the same scripts, not written):", sorted(set(discarded)))
missing = sorted({t for _, m in PLAN for t in m.values()} - {t for t, _ in written})
print("MISSING:", missing if missing else "none")
