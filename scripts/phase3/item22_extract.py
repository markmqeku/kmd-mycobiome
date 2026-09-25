#!/usr/bin/env python
"""Item 22. Extract text from Swanepoel et al. (2018) so every figure quoted in the
manuscript can be attributed to a located passage rather than to recollection."""
import re, os
from pypdf import PdfReader
P = ("<USER_HOME>/Downloads/Swanepoeletal.2018-FirstreportofanewmalformationdiseaseofcommonkareeSearsialanceainSouthAfrica.pdf")
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
os.makedirs(C, exist_ok=True)
r = PdfReader(P)
pages = [(i+1, (p.extract_text() or "")) for i, p in enumerate(r.pages)]
full = "\n".join(f"\n===== PAGE {n} =====\n{t}" for n, t in pages)
open(f"{C}/swanepoel2018_text.txt", "w", encoding="utf-8").write(full)
print(f"pages: {len(pages)}   characters: {len(full):,}")
meta = r.metadata or {}
for k in ("/Title", "/Author", "/CreationDate"):
    if meta.get(k): print(f"{k}: {str(meta[k])[:120]}")
print()
for n, t in pages:
    head = " ".join(t.split())[:150]
    print(f"  p{n}: {head}")
