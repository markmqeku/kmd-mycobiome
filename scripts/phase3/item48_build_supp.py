#!/usr/bin/env python
"""Item 48. Supplementary document for review: the four supplementary figures and the eleven
supplementary tables, numbered to match the citations in the main text, each with its caption.

Large tables are truncated in this review copy, with the omission stated explicitly on the
page. The underlying .tsv, .xlsx and .fasta files remain the authoritative versions.
"""
import csv, os, re
from docx import Document
from docx.shared import Pt, Cm, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH, WD_LINE_SPACING
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.oxml.ns import qn

K = r"<KMD_ROOT>"
O = os.path.join(K, "__reanalysis_2026-06", "phase3", "outputs_27-07-2026")
C = os.path.join(O, "clean")
OUT = os.path.join(K, "MANUSCRIPT_KMD_supplementary.docx")
MAXROWS = 25

FIGS = [
 ("Figure S1", os.path.join(C, "figures", "FigureS1_rarefaction_core16.png"),
  "Rarefaction curves for the 16 core samples on the corrected fungi-only table. The dashed "
  "line marks the rarefaction depth of 57,124 reads. Cited in Methods 2.8."),
 ("Figure S2", os.path.join(O, "figures", "FigureS2_dada2_retention.png"),
  "DADA2 read retention by stage across all 45 samples. Cited in Methods 2.3."),
 ("Figure S3", os.path.join(O, "figures", "FigureS3_assignment_by_rank.png"),
  "Taxonomy assignment rate by rank under UNITE 8.2, UNITE 10.0 and the rank-threshold "
  "scheme. Cited in Methods 2.5."),
 ("Figure S4", os.path.join(C, "figures", "FigureS4_guild_coverage.png"),
  "Guild assignment coverage by tool, on 1,142 fungal ASVs. Cited in Results 3.9."),
]
TABS = [
 ("Table S1", os.path.join(C, "asv_catalogue", "ASV_catalogue.tsv"),
  "ASV catalogue: 1,142 fungal ASVs with final taxonomy, per-sample counts and representative "
  "sequences. Cited in Methods 2.9."),
 ("Table S2", os.path.join(C, "TableS2_bloemfontein_inventory.tsv"),
  "Bloemfontein presence and absence inventory, 199 taxa across the six samples with host load "
  "at or below 10 percent. Distinct ASVs and (sample, ASV) occurrences are given separately. "
  "Cited in Results 3.8."),
 ("Table S3", os.path.join(C, "TableS3_taxonomy_thresholds.tsv"),
  "Taxonomy assignment rate at each rank under three schemes, n = 1,142. Cited in Methods 2.5."),
 ("Table S4", os.path.join(O, "assembly", "TableS4_dada2_stats.tsv"),
  "DADA2 statistics for all 45 denoised samples. Cited in Methods 2.3."),
 ("Table S5", os.path.join(C, "guild_assignments_per_ASV.tsv"),
  "Guild assignments per ASV with confidence rankings, 1,142 ASVs. Cited in Results 3.9."),
 ("Table S6", os.path.join(O, "placement", "placement_results.tsv"),
  "Placement results for the 28 abundant unnamed ASVs. Cited in Results 3.10."),
 ("Table S7", os.path.join(O, "placement", "reference_provenance.tsv"),
  "Reference provenance for the placement analysis: accession, name and type status. "
  "Cited in Results 3.10."),
 ("Table S8", os.path.join(C, "TableS8_hill_numbers.tsv"),
  "Hill numbers per sample at depth 57,124. Cited in Results 3.5."),
 ("Table S9", os.path.join(O, "full_screen", "screen_calls.tsv"),
  "Non-fungal screen: call, evidence and reason for each of the 277 ASVs sent to NCBI. "
  "Cited in Methods 2.4 and Results 3.11."),
 ("Table S10", os.path.join(O, "full_screen", "screen_basis.tsv"),
  "Screen basis for all 1,173 ASVs, giving the reason each was or was not sent to NCBI. "
  "Cited in Methods 2.4."),
 ("Table S11", os.path.join(O, "full_screen", "inconclusive_asvs.tsv"),
  "The 21 ASVs retained as inconclusive. Cited in Results 3.11."),
]

def set_font(run, **kw):
    run.font.name = "Times New Roman"; run.font.size = Pt(kw.get("size", 12))
    run._element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
    if "bold" in kw: run.bold = kw["bold"]
    if "italic" in kw: run.italic = kw["italic"]

def add_caption(doc, num, text):
    p = doc.add_paragraph(); p.paragraph_format.space_after = Pt(10)
    set_font(p.add_run(f"{num}. "), bold=True)
    # carry markdown emphasis into real italics
    for part in re.split(r"(\*[^*]+\*)", text):
        if not part: continue
        if part.startswith("*") and part.endswith("*"):
            set_font(p.add_run(part[1:-1]), italic=True)
        else:
            set_font(p.add_run(part))

doc = Document()
s = doc.styles["Normal"]; s.font.name = "Times New Roman"; s.font.size = Pt(12)
s.element.rPr.rFonts.set(qn("w:eastAsia"), "Times New Roman")
sec = doc.sections[0]
sec.top_margin = sec.bottom_margin = sec.left_margin = sec.right_margin = Cm(2.0)

h = doc.add_heading(level=1)
r = h.add_run("Supplementary material"); r.font.name = "Times New Roman"
r.font.size = Pt(16); r.font.bold = True; r.font.color.rgb = RGBColor(0, 0, 0)
p = doc.add_paragraph()
set_font(p.add_run(
  "Companion to the manuscript. All items are built on the corrected 1,142-ASV fungal table. "
  "Tables are reproduced here for review; the .tsv, .xlsx and .fasta files supplied alongside "
  "are the authoritative versions and are complete."))

missing = []
doc.add_page_break()
hh = doc.add_heading(level=2); rr = hh.add_run("Supplementary figures")
rr.font.name = "Times New Roman"; rr.font.size = Pt(14); rr.font.bold = True
rr.font.color.rgb = RGBColor(0, 0, 0)
for num, path, cap in FIGS:
    pp = doc.add_paragraph(); pp.alignment = WD_ALIGN_PARAGRAPH.CENTER
    if os.path.exists(path):
        try: pp.add_run().add_picture(path, width=Cm(16))
        except Exception as e:
            set_font(pp.add_run(f"[NOT EMBEDDED: {num}, {type(e).__name__}]"), bold=True)
            missing.append((num, str(e)))
    else:
        set_font(pp.add_run(f"[FILE NOT FOUND: {num}]"), bold=True); missing.append((num, "missing"))
    add_caption(doc, num, cap)

doc.add_page_break()
hh = doc.add_heading(level=2); rr = hh.add_run("Supplementary tables")
rr.font.name = "Times New Roman"; rr.font.size = Pt(14); rr.font.bold = True
rr.font.color.rgb = RGBColor(0, 0, 0)
counts = []
for num, path, cap in TABS:
    add_caption(doc, num, cap)
    if not os.path.exists(path):
        set_font(doc.add_paragraph().add_run(f"[FILE NOT FOUND: {path}]"), bold=True)
        missing.append((num, "missing")); continue
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = [r for r in csv.reader(fh, delimiter="\t") if r and not r[0].startswith("#")]
    total = len(rows) - 1
    show = rows[:MAXROWS+1]
    ncol = min(max(len(r) for r in show), 8)
    t = doc.add_table(rows=len(show), cols=ncol)
    t.style = "Table Grid"; t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for ri, row in enumerate(show):
        for ci in range(ncol):
            cp = t.cell(ri, ci).paragraphs[0]
            cp.paragraph_format.space_after = Pt(0)
            cp.paragraph_format.line_spacing_rule = WD_LINE_SPACING.SINGLE
            val = row[ci] if ci < len(row) else ""
            set_font(cp.add_run(val[:38]), size=8, bold=(ri == 0))
    note = doc.add_paragraph()
    if total > MAXROWS:
        set_font(note.add_run(
            f"Showing the first {MAXROWS} of {total:,} data rows"
            + (f", and the first {ncol} of {max(len(r) for r in rows)} columns" if max(len(r) for r in rows) > ncol else "")
            + ". The complete table is supplied as a separate file."), size=9, italic=True)
    else:
        set_font(note.add_run(f"Complete: {total:,} data rows."), size=9, italic=True)
    counts.append((num, total))
    doc.add_paragraph()

doc.save(OUT)
print(f"wrote {OUT}")
print(f"   supplementary figures: {len(FIGS)}   supplementary tables: {len(TABS)}")
for n, c in counts: print(f"      {n:11} {c:>8,} data rows")
print("   " + ("all present" if not missing else f"MISSING: {missing}"))
