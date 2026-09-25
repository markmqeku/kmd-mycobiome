#!/usr/bin/env python
"""67c. Account for the difference between the markdown word count (10,687 by a naive split)
and Word's count (10,497). Establish whether any body prose was lost in conversion, or whether
the difference is markdown syntax that Word never sees."""
import io, os, re
from docx import Document
K = r"<KMD_ROOT>"
md = io.open(os.path.join(K, "MANUSCRIPT_KMD_mycobiome.md"), encoding="utf-8").read()
doc = Document(os.path.join(K, "MANUSCRIPT_KMD_mycobiome.docx"))

naive = len(md.split())
# strip what Word never receives as words
s = md
s = re.sub(r"^\s*\|[\s\-:|]+\|\s*$", " ", s, flags=re.M)   # table separator rows
s = re.sub(r"^\s*---\s*$", " ", s, flags=re.M)             # horizontal rules
s = s.replace("|", " ")                                     # table cell pipes
s = re.sub(r"[*`]", "", s)                                  # emphasis markers
s = re.sub(r"^#{1,6}\s*", "", s, flags=re.M)                # heading hashes
cleaned = len(s.split())
print(f"markdown, naive split                : {naive:,}")
print(f"markdown, syntax removed             : {cleaned:,}")
print(f"  accounted for by markdown syntax   : {naive - cleaned:,}")

docx_words = []
for p in doc.paragraphs: docx_words += p.text.split()
for t in doc.tables:
    for row in t.rows:
        for c in row.cells:
            for p in c.paragraphs: docx_words += p.text.split()
print(f"docx, counted the same way           : {len(docx_words):,}")
print(f"Word's own ComputeStatistics figure  : 10,497")
print(f"\nresidual, cleaned markdown vs docx   : {cleaned - len(docx_words):+,}")

# does any sentence of body prose exist in the markdown but not the docx?
dtext = " ".join(docx_words)
dnorm = re.sub(r"\s+", " ", dtext)
body_md = md.split("## References")[0]
body_md = re.sub(r"^\s*\|.*$", "", body_md, flags=re.M)     # tables checked separately
sentences = [x.strip() for x in re.split(r"(?<=[.;:])\s+", re.sub(r"[*`#]", "", body_md))
             if len(x.split()) >= 6]
missing = []
for s_ in sentences:
    probe = " ".join(re.sub(r"\s+", " ", s_).split()[:9])
    if probe and probe not in dnorm: missing.append(s_[:120])
print(f"\nbody sentences of 6+ words checked   : {len(sentences)}")
print(f"not found in the docx                : {len(missing)}")
for m in missing[:20]: print(f"   MISSING: {m}...")
print("\nCONCLUSION: " + ("no body prose was lost" if not missing
      else f"{len(missing)} passage(s) need checking"))
