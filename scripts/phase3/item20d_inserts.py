#!/usr/bin/env python
"""20d. Final INSERT checklist for Mark, rebuilt from the updated drafts."""
import re, os, csv, io
K = "<KMD_ROOT>"
A = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
os.makedirs(A, exist_ok=True)
PAT = re.compile(r"\\?\[INSERT:?\s*([^\]]*)\]")
rows = []
for p in ("DRAFT_METHODS.md", "DRAFT_RESULTS.md"):
    for i, line in enumerate(io.open(f"{K}/{p}", encoding="utf-8", errors="replace"), 1):
        for m in PAT.finditer(line):
            item = " ".join(m.group(1).split())
            # a bare [INSERT] carries no text of its own; take the sentence around it so the
            # checklist entry is actionable without opening the file
            ctx = " ".join(line[max(0, m.start()-150):m.end()+60].split())
            rows.append((p, i, item or "(no text) " + ctx, ctx if not item else ""))
with open(f"{A}/INSERT_checklist.tsv", "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["file", "line", "outstanding_item", "context_if_untitled"]); w.writerows(rows)
print(f"{len(rows)} outstanding INSERTs -> {A}/INSERT_checklist.tsv\n")
for f, i, t, c in rows:
    print(f"{f}:{i}")
    print(f"   {t}\n")
