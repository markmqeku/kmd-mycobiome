#!/usr/bin/env python
"""Item 41 formatting, second pass. Section 3.10 lists reference taxa that are italicised
nowhere in the manuscript, so the conservative first pass could not reach them. Each name
below was read out of the placement results and confirmed to be a genus or binomial, not a
family, order or place name. Families and higher ranks (Didymellaceae, Coniochaetales,
Sclerotiniaceae, Saccotheciaceae, Dothideaceae, Mortierellaceae, Microbotryomycetes,
Dothideomycetes) are correctly left in roman."""
import io, re
M = "<KMD_ROOT>/MANUSCRIPT_KMD_mycobiome.md"
d = io.open(M, encoding="utf-8").read()
head, sep, tail = d.partition("## References")
assert sep

# binomials first, then bare genera used with "sp."
NAMES = ["Neodidymelliopsis tinkyukuku", "Nothophoma brennandiae", "Gamsiella stylospora",
         "Coniochaeta luteoviridis", "Endoconidioma populi", "Coniochaeta africana",
         "Paivomyces mimosae", "Botrytis cinerea",
         "Botryosphaeria", "Endoconidioma", "Epicoccum"]
changes = []
for t in NAMES:
    pat = re.compile(r"(?<![*\w])" + re.escape(t) + r"(?![*\w])")
    def repl(m, t=t):
        changes.append((t, " ".join(head[max(0, m.start()-55):m.start()+55].split())))
        return "*" + m.group(0) + "*"
    head = pat.sub(repl, head)
io.open(M, "w", encoding="utf-8").write(head + sep + tail)
print(f"italicised {len(changes)} further taxon names in the body:")
for t, ctx in changes: print(f"   {t:30} ...{ctx}...")
