#!/usr/bin/env python
"""Item 37 to 40. Assemble one manuscript from the five drafts. The five drafts and both
superseded manuscripts are NOT touched; this writes a new file only.

Supplied prose is copied verbatim. The only permitted changes are (a) unescaping the
backslash-escaped markdown that two drafts carry from a docx conversion, and (b) inserting
figure and table callouts. Every callout insertion is anchored to an exact string and the
script FAILS LOUDLY if an anchor is not found, so a callout can never land in the wrong place
or be silently dropped."""
import csv, io, os, re, sys
K = "<KMD_ROOT>"
O = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026"; C = f"{O}/clean"
OUT = f"{K}/MANUSCRIPT_KMD_mycobiome.md"

def unescape(t):
    # the docx conversion escaped markdown syntax; strip only those escapes
    for a, b in ((r"\#", "#"), (r"\*", "*"), (r"\[", "["), (r"\]", "]"), (r"\_", "_")):
        t = t.replace(a, b)
    return t

S = {}
for n in ("TITLE_ABSTRACT", "INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION"):
    S[n] = unescape(io.open(f"{K}/DRAFT_{n}.md", encoding="utf-8").read())

# ---- callout insertions: (section, exact anchor, replacement) ----
INS = [
 ("METHODS",
  "Feature tables and representative sequences were merged after denoising.",
  "Per-sample read retention at each denoising stage is given in Supplementary Table S4 and "
  "Supplementary Figure S2. Feature tables and representative sequences were merged after denoising."),
 ("METHODS",
  "Thirty-one ASVs were confirmed non-fungal and removed:",
  "The basis on which each of the 1,173 ASVs was or was not sent to NCBI is given in "
  "Supplementary Table S10, and the per-ASV calls in Supplementary Table S9. "
  "Thirty-one ASVs were confirmed non-fungal and removed:"),
 ("METHODS",
  "Diversity is reported as Hill numbers (Chao",
  "Rarefaction curves for the 16 core samples are shown in Supplementary Figure S1. "
  "Diversity is reported as Hill numbers (Chao"),
 ("METHODS",
  "All processing parameters, software versions, database releases and per-sample read statistics",
  "The complete ASV catalogue is provided as Supplementary Table S1. "
  "All processing parameters, software versions, database releases and per-sample read statistics"),
 ("RESULTS",
  "Diversity is reported as Hill numbers, which express diversity as an effective number of taxa.",
  "Per-sample values are given in Supplementary Table S8 and shown in Figure 5. "
  "Diversity is reported as Hill numbers, which express diversity as an effective number of taxa."),
 ("RESULTS",
  "(Supplementary Table X)",
  "(Supplementary Table S2)"),
 ("RESULTS",
  "Guild assignment was attempted with FungalTraits and FUNGuild.",
  "Guild assignment was attempted with FungalTraits and FUNGuild (Supplementary Table S5, "
  "Supplementary Figure S4)."),
 ("RESULTS",
  "Twenty-four of the 28 were placed within a neighbourhood.",
  "Twenty-four of the 28 were placed within a neighbourhood (Figure 8; Supplementary Tables S6 "
  "and S7)."),
 ("RESULTS",
  "All 1,173 ASVs remaining after host removal were screened.",
  "All 1,173 ASVs remaining after host removal were screened (Supplementary Tables S9 and S10)."),
 ("RESULTS",
  "These were retained and are listed as inconclusive in the supplementary material.",
  "These were retained and are listed as inconclusive in Supplementary Table S11."),
 ("METHODS",
  "Two classification controls were run.",
  "Assignment rates at each rank under both releases and under the rank-threshold scheme are "
  "given in Supplementary Table S3 and Supplementary Figure S3. Two classification controls were run."),
]
for sec, anchor, repl in INS:
    if anchor not in S[sec]:
        sys.exit(f"ABORT: anchor not found in {sec}: {anchor[:80]!r}")
    if S[sec].count(anchor) != 1:
        sys.exit(f"ABORT: anchor is not unique in {sec} ({S[sec].count(anchor)}x): {anchor[:60]!r}")
    S[sec] = S[sec].replace(anchor, repl, 1)
print(f"callouts inserted: {len(INS)}, every anchor matched exactly once")

# ---- strip the per-draft provenance headers, keep the prose ----
def body(t, drop_first_heading=True):
    lines = t.splitlines()
    out, started = [], False
    for l in lines:
        if not started:
            if l.startswith("# ") or l.strip() == "" or l.startswith("# Source"):
                continue
            started = True
        out.append(l)
    return "\n".join(out).strip()

title_block = S["TITLE_ABSTRACT"]
m_title = re.search(r"## Title\s*\n+(.+?)\n\s*\n", title_block, re.S)
m_run = re.search(r"## Running title\s*\n+(.+?)\n\s*\n", title_block, re.S)
m_abs = re.search(r"## Abstract\s*\n+(.+?)\n\s*\n## Keywords", title_block, re.S)
m_kw = re.search(r"## Keywords\s*\n+(.+)$", title_block, re.S)
for nm, mm in (("Title", m_title), ("Running title", m_run), ("Abstract", m_abs), ("Keywords", m_kw)):
    if not mm: sys.exit(f"ABORT: could not parse {nm} from DRAFT_TITLE_ABSTRACT.md")
TITLE = " ".join(m_title.group(1).split())
RUN = " ".join(m_run.group(1).split())
ABS = " ".join(m_abs.group(1).split())
KW = " ".join(m_kw.group(1).split())

FIGS = []
for r in csv.DictReader(open(f"{C}/figure_table_manifest.tsv"), delimiter="\t"):
    FIGS.append(r)

parts = []
parts.append(f"# {TITLE}\n")
# Author block resolved by Mark with Prof. Gryzenhout (C-17): Soumya Ghosh is not an author on
# this manuscript; the surname is Cason, not Casson or Casen; Rosemary is the given name and
# Tonjock Kinge the surname. Baked in so no rebuild can erase it.
AUTHORS = ("Mark Mqeku(1)*, Errol Cason(2), Rosemary Tonjock Kinge(1,3), "
           "Marieka Gryzenhout(1)")
AFFIL = [
 "(1) Department of Genetics, Faculty of Natural and Agricultural Sciences, University of the "
 "Free State, Bloemfontein, South Africa",
 "(2) Department of Agriculture, Faculty of Natural and Agricultural Sciences, University of the "
 "Free State, Bloemfontein, South Africa",
 "(3) Department of Biological Sciences, Faculty of Sciences, University of Bamenda, P.O. Box 39, "
 "Bambili, North West Region, Cameroon",
]
parts.append("**Authors.** " + AUTHORS + "\n")
for _a in AFFIL: parts.append(_a + "\n")
parts.append("\\* Corresponding author. Mark Mqeku, <email>\n")
parts.append("[INSERT: ORCID identifiers, if the journal requires them. None appears in any "
             "source copy.]\n")
parts.append(f"**Running title.** {RUN}\n")
parts.append("---\n")
parts.append("## Abstract\n\n" + ABS + "\n")
parts.append("**Keywords.** " + KW + "\n")
parts.append("---\n")
parts.append(body(S["INTRODUCTION"]))
parts.append("\n---\n")
parts.append(body(S["METHODS"]))
parts.append("\n---\n")
parts.append(body(S["RESULTS"]))
parts.append("\n---\n")
parts.append(body(S["DISCUSSION"]))

doc = "\n\n".join(parts)
doc = re.sub(r"\n{4,}", "\n\n\n", doc)
io.open(OUT, "w", encoding="utf-8").write(doc + "\n")
print(f"\nassembled -> {OUT}  ({len(doc.splitlines())} lines, {len(doc.split()):,} words)")

# ---- verify every figure and table is cited ----
print("\n=== callout audit ===")
missing = []
for r in FIGS:
    num = r["number"]
    if num == "Data": continue
    # allow the plural form: "Supplementary Tables S6 and S7" cites both
    kind, ident = num.split(None, 1)
    pat = rf"(Supplementary\s+)?{kind}s?\s+(S?\d+\s+and\s+)?{re.escape(ident)}\b"
    cited = re.search(pat, doc)
    if not cited: missing.append(num)
    print(f"   {num:11} {'cited' if cited else 'NOT CITED':10} {r['title'][:56]}")
print(f"\nuncited: {missing if missing else 'none'}")
