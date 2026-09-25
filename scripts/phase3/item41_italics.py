#!/usr/bin/env python
"""Item 41 formatting fix: italicise taxon names in the manuscript BODY only.

Scope is deliberately narrow. Only names already italicised elsewhere in the body are
touched, so this cannot invent a taxon or italicise a common noun. The reference list and
the supplementary manifest are left untouched, because those entries are transcribed
bibliographic text and a file inventory, not prose. Every change is printed."""
import io, re
M = "<KMD_ROOT>/MANUSCRIPT_KMD_mycobiome.md"
d = io.open(M, encoding="utf-8").read()
head, sep, tail = d.partition("## References")
assert sep, "could not locate the reference list boundary"

# longest first so binomials are handled before their genus
TAXA = ["Searsia lancea", "Didymella maydis", "Curvibasidium", "Fusarium"]
changes = []
for t in TAXA:
    # skip anything already inside asterisks, and skip markdown headings
    pat = re.compile(r"(?<![*\w])" + re.escape(t) + r"(?![*\w])")
    out, last, buf = [], 0, head
    def repl(m):
        s = m.start()
        line_start = buf.rfind("\n", 0, s) + 1
        if buf[line_start:line_start+2] in ("# ", "##"): return m.group(0)   # leave headings
        changes.append((t, " ".join(buf[max(0, s-60):s+60].split())))
        return "*" + m.group(0) + "*"
    head = pat.sub(repl, buf)

io.open(M, "w", encoding="utf-8").write(head + sep + tail)
print(f"italicised {len(changes)} occurrences in the body:")
for t, ctx in changes: print(f"   {t:18} ...{ctx}...")
print("\nreference list and supplementary manifest left untouched")
