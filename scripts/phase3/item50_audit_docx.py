#!/usr/bin/env python
"""Item 50. Audit the .docx CONTENT, not the markdown. Conversion can silently drop italics,
lose a figure, or introduce dashes, so every check reads the Word file itself. Reports losses;
fixes nothing."""
import os, re, zipfile
from docx import Document
from docx.oxml.ns import qn

K = r"<KMD_ROOT>"
DOCX = os.path.join(K, "MANUSCRIPT_KMD_mycobiome.docx")
d = Document(DOCX)
fails = []

# gather every run with its formatting, from body paragraphs AND table cells
runs = []
def harvest(paras, where):
    for p in paras:
        for r in p.runs:
            runs.append((r.text, bool(r.italic), bool(r.bold), where))
harvest(d.paragraphs, "body")
for t in d.tables:
    for row in t.rows:
        for cell in row.cells:
            harvest(cell.paragraphs, "table")
text = "".join(r[0] for r in runs)
body_text = "\n".join(p.text for p in d.paragraphs)
all_text = text

print("=== 1. the 13 headline numbers ===")
MUST = {"1,142":"fungal ASVs","2,517,177":"fungal reads","57,124":"rarefaction depth",
        "199":"Bloemfontein taxa","71":"yeast lineage","40 of 45":"lineage occurrence",
        "84.3":"host percent","37.1":"host share","1,322":"total ASVs","149":"host ASVs",
        "4,010,617":"total reads","93.2":"missed Searsia","85.1 to 85.5":"nearest named"}
for v, w in MUST.items():
    ok = v in all_text
    print(f"   {'ok  ' if ok else 'LOST':6} {v:12} {w}")
    if not ok: fails.append(f"number lost in conversion: {v} ({w})")

print("\n=== 2. taxa italicised in the .docx ===")
TAXA = ["Searsia lancea","Curvibasidium","Homophron spadiceum","Amaranthus","Fusarium",
        "Filobasidium","Cryptococcus","Didymella","Mycosphaerella","Filobasidiella",
        "Dysphania melanocarpa","Vigna unguiculata","Botrytis cinerea","Gamsiella stylospora",
        "Epicoccum","Botryosphaeria","Paivomyces mimosae","Endoconidioma"]
for t in TAXA:
    ital = sum(1 for tx, it, bd, wh in runs if it and t in tx)
    plain = sum(1 for tx, it, bd, wh in runs if not it and re.search(r"(?<![A-Za-z])"+re.escape(t)+r"(?![A-Za-z])", tx))
    status = "ok" if ital and not plain else ("CHECK" if ital else "NOT ITALIC")
    print(f"   {status:11} {t:24} italic runs {ital:>3}, plain {plain:>3}")
    if ital == 0 and plain: fails.append(f"{t} not italicised anywhere in the docx")

print("\n=== 3. figures and tables ===")
with zipfile.ZipFile(DOCX) as z:
    imgs = [n for n in z.namelist() if n.startswith("word/media/")]
print(f"   embedded image parts: {len(imgs)}  (expected 9)")
if len(imgs) != 9: fails.append(f"expected 9 embedded figures, found {len(imgs)}")
for n in [f"Figure {i}" for i in range(1, 10)]:
    cited = n in all_text
    print(f"   {'ok  ' if cited else 'LOST':6} {n} cited")
    if not cited: fails.append(f"{n} not cited in the docx")
for n in [f"Figure S{i}" for i in range(1, 5)] + [f"Table S{i}" for i in range(1, 12)]:
    kind, ident = n.split()
    # "Tables S6 and S7" cites both; the identifier can sit on either side of "and"
    cited = (re.search(re.escape(n) + r"(?!\d)", all_text)
             or re.search(rf"{kind}s\s+(?:S\d+\s+and\s+)?{ident}(?!\d)", all_text)
             or re.search(rf"{kind}s\s+{ident}(?!\d)\s+and\s+S\d+", all_text))
    print(f"   {'ok  ' if cited else 'LOST':6} {n} cited")
    if not cited: fails.append(f"{n} not cited in the docx")
print(f"   native Word tables: {len(d.tables)}")
if len(d.tables) == 0: fails.append("no native tables found")

print("\n=== 4. section numbering ===")
heads = [p.text.strip() for p in d.paragraphs if p.style.name.startswith("Heading")]
for prefix, label in (("2.", "Methods"), ("3.", "Results"), ("5.", "Discussion")):
    got = [h for h in heads if re.match(r"^" + re.escape(prefix) + r"\d", h)]
    nums = [int(re.match(r"^\d+\.(\d+)", h).group(1)) for h in got]
    seq = nums == list(range(1, len(nums)+1))
    print(f"   {label}: {len(got)} subsections, sequential {'yes' if seq else 'NO ' + str(nums)}")
    if not seq: fails.append(f"{label} numbering not sequential in the docx")

print("\n=== 5. dashes introduced by conversion ===")
em, en = all_text.count("\u2014"), all_text.count("\u2013")
print(f"   em {em}, en {en}   {'ok' if em == en == 0 else 'FAIL'}")
if em or en: fails.append(f"dashes present in the docx: em {em}, en {en}")

print("\n=== 6. INSERTs visible and highlighted ===")
n_ins = 0
for p in d.paragraphs:
    for r in p.runs:
        if "[INSERT" in r.text:
            hl = r._element.rPr is not None and r._element.rPr.find(qn("w:highlight")) is not None
            n_ins += 1
            print(f"   {'ok  ' if (hl and r.bold) else 'CHECK':6} bold={bool(r.bold)} "
                  f"highlight={'yes' if hl else 'NO'}  {r.text[:82]}")
            if not (hl and r.bold): fails.append("an INSERT is not both bold and highlighted")
print(f"   INSERT markers found: {n_ins}")
if n_ins < 6: fails.append(f"expected at least 6 INSERTs, found {n_ins}")

print("\n=== 6b. title, author line and the C-17 authorship decisions ===")
for what, s_ in (("title", "An unnamed Tremellomycetes lineage dominates the karee mycobiome"),
                 ("author line", "Mark Mqeku"),
                 ("Cason spelling", "Errol Cason"),
                 ("Kinge surname", "Rosemary Tonjock Kinge"),
                 ("corresponding author", "Mark Mqeku, <email>")):
    ok = s_ in all_text
    print(f"   {'ok  ' if ok else 'MISSING':8} {what}: {s_[:62]}")
    if not ok: fails.append(f"{what} missing from the docx")
for bad, why in (("Soumya Ghosh", "C-17 removed this author"),
                 ("Casson", "superseded spelling of Cason"),
                 ("Casen", "superseded spelling of Cason")):
    if bad in all_text:
        print(f"   PRESENT  {bad}: {why}"); fails.append(f"{bad} still present: {why}")
    else:
        print(f"   ok       {bad} correctly absent")

print("\n=== 7. page setup ===")
s = d.sections[0]
sectPr = s._sectPr
ln = sectPr.find(qn("w:lnNumType"))
print(f"   margins: {s.top_margin.cm:.1f} cm all round")
print(f"   line numbering: {'present, continuous' if ln is not None else 'ABSENT'}")
print(f"   footer page field: {'present' if 'PAGE' in s.footer.paragraphs[0]._element.xml else 'ABSENT'}")
if ln is None: fails.append("line numbering absent")
norm = d.styles["Normal"]
print(f"   Normal font: {norm.font.name} {norm.font.size.pt:.0f} pt")

print("\n=== SUMMARY ===")
print("DOCX AUDIT CLEAN" if not fails else f"{len(fails)} issue(s):")
for f in fails: print(f"   - {f}")
