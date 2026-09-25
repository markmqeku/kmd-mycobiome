#!/usr/bin/env python
"""Item 65a. Measure how the Discussion's words are actually distributed across themes.
Each paragraph is assigned to the theme it is mainly about, judged by which theme's terms
dominate it. Paragraphs are the unit, so no word is counted twice."""
import io, re
K = "<KMD_ROOT>"
t = io.open(f"{K}/DRAFT_DISCUSSION.md", encoding="utf-8").read()

THEMES = {
 "host DNA":            r"host DNA|host sequence|84\.3|host load|host content|co-amplif|plant DNA",
 "yeast lineage":       (r"Tremellomycetes|yeast|Filobasid|Cryptococcus|unnamed|unresolved lineage"
                         r"|LSU|D1/D2|could not be identified|cannot be resolved|named references"
                         r"|environmental sequence|never been named|reference set|reference database"),
 "community findings":  r"Curvibasidium|diversity|ordination|richness|Didymella|Mycosphaerella|shared|composition|site effect|geographic|populations",
 "other failure modes": r"truncation|read length|denois|identity threshold|rank truncation|guild|FUNGuild|FungalTraits|contamination screen|Anacardiaceae|agaric|maximum reported confidence|wrong class|wrong branch",
 "limitations, general": r"sample size|pooled|inferential|replicate|pruning|age was|singleton|coverage-based|constraints on this study",
 # deliberately narrow: this theme must not absorb paragraphs that merely mention a database
 "framing and future":  r"future work|three lines of work|recommend|prepared for deposition|reproducib|central to it|remains unidentified",
}
paras, cur = [], []
sec = "preamble"
for line in t.splitlines():
    s = line.strip()
    if s.startswith("###"):
        if cur: paras.append((sec, " ".join(cur))); cur = []
        sec = s.lstrip("# ").strip(); continue
    if not s:
        if cur: paras.append((sec, " ".join(cur))); cur = []
        continue
    cur.append(s)
if cur: paras.append((sec, " ".join(cur)))

tot = 0
by_theme, by_sec = {}, {}
rows = []
for sec, p in paras:
    if p.startswith("#") or len(p.split()) < 8: continue
    n = len(p.split()); tot += n
    scores = {k: len(re.findall(v, p, re.I)) for k, v in THEMES.items()}
    best = max(scores, key=scores.get)
    if scores[best] == 0: best = "framing and future"
    by_theme[best] = by_theme.get(best, 0) + n
    by_sec.setdefault(sec, {}).setdefault(best, 0)
    by_sec[sec][best] += n
    rows.append((sec, best, n, p[:70]))

print(f"=== 65a. Discussion, {tot:,} words across {len(rows)} paragraphs ===\n")
print(f"{'theme':24} {'words':>7} {'share':>8}")
for k, v in sorted(by_theme.items(), key=lambda x: -x[1]):
    print(f"   {k:22} {v:>7,} {100*v/tot:>7.1f}%")
print("\nby section:")
for sec in by_sec:
    st = sum(by_sec[sec].values())
    print(f"\n   {sec}  ({st:,} words)")
    for k, v in sorted(by_sec[sec].items(), key=lambda x: -x[1]):
        print(f"      {k:22} {v:>6,} {100*v/st:>6.1f}% of the section")
print("\nparagraph assignments:")
for sec, best, n, snip in rows:
    print(f"   {sec[:26]:28} {best:22} {n:>4}w  {snip}...")
