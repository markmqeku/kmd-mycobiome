#!/usr/bin/env python
"""69-ADD continued. The four copies differ by checksum. Compare the author block, the
corresponding-author markers and the sampling and extraction sentences across all four, so the
authoritative version can be identified rather than guessed."""
import os
from docx import Document
K = r"<KMD_ROOT>"
NAME = "Article1_Bloemfontein_Baseline_Fungal_Communities.docx"
DIRS = ("files", "files2", "files3", "files4")

def paras(p):
    return [x.text.strip() for x in Document(p).paragraphs if x.text.strip()]

store = {}
for d in DIRS:
    p = os.path.join(K, d, NAME)
    if os.path.exists(p): store[d] = paras(p)

print("=== author line in each copy ===")
for d, ps in store.items():
    print(f"   {d:7} {ps[1][:180]}")
print("\n=== affiliations in each copy ===")
for d, ps in store.items():
    affs = [x for x in ps[2:8] if x[:1].isdigit()]
    print(f"   {d}:")
    for a in affs: print(f"      {a[:150]}")

print("\n=== marker footnotes: what * and other symbols mean ===")
import re
for d, ps in store.items():
    hits = [x for x in ps if re.match(r"^\s*[*¥†‡#]", x) or re.search(r"[Cc]orrespond", x)]
    print(f"   {d}: {hits[:3] if hits else 'no marker footnote found'}")

print("\n=== sampling and extraction sentences ===")
for key, pat in (("trees sampled", r"branches from|different Searsia lancea trees|trees were collected"),
                 ("extraction mass", r"illigram|grams of pulverized|mg of"),
                 ("sample count", r"32 samples|samples from Bloemfontein")):
    print(f"\n   -- {key} --")
    for d, ps in store.items():
        m = [x for x in ps if re.search(pat, x, re.I)]
        print(f"      {d:7} {(m[0][:165] if m else '(not found)')}")

print("\n=== do the four differ in body text at all? ===")
base = store[DIRS[0]]
for d in list(store)[1:]:
    other = store[d]
    diff = [i for i in range(min(len(base), len(other))) if base[i] != other[i]]
    print(f"   {DIRS[0]} vs {d}: {len(base)} vs {len(other)} paragraphs, "
          f"{len(diff)} differing paragraph(s)"
          + (f", first at index {diff[0]}" if diff else ""))
    for i in diff[:3]:
        print(f"      [{i}] A: {base[i][:110]}")
        print(f"      [{i}] B: {other[i][:110]}")
