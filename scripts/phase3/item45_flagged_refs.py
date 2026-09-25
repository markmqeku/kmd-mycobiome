#!/usr/bin/env python
"""Item 45. Build Mark's manual verification checklist for the flagged references, grouped by
what each one actually needs. Nothing is guessed: where Crossref returned nothing, the row
says so."""
import csv, io, re
K = "<KMD_ROOT>"
C = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
OUT = f"{K}/REFERENCES_TO_VERIFY.md"
rows = list(csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t"))
flag = [r for r in rows if r["status"] != "VERIFIED"]
ver = [r for r in rows if r["status"] == "VERIFIED"]

def group(r):
    """Order matters. A 'wrong record entirely' reason often also mentions the year, so the
    title tests must be evaluated before the year test or everything collapses into A."""
    w = (r["reason_if_flagged"] or "").lower()
    if "secondary" in w: return "C"
    if ("title does not" in w or "no reliable match" in w or "no record" in w
            or "neither year nor title" in w): return "B"
    if "review record" in w or "preprint" in w or "no journal" in w: return "A"
    if "year" in w: return "A"
    return "B"
G = {"A": [], "B": [], "C": []}
for r in flag: G[group(r)].append(r)

L = []
A = L.append
A("# References requiring manual verification\n")
A(f"Working checklist for the {len(flag)} flagged entries in `MANUSCRIPT_KMD_mycobiome.md`. "
  f"The other {len(ver)} are verified and need no action.\n")
A("**Why this list exists.** Automated matching against Crossref is unreliable for this "
  "reference set. A returned DOI is not evidence of a correct match: two entries initially "
  "passed a year-plus-title test while Crossref had actually returned a post-publication "
  "review and a peer-review record rather than the papers themselves. Every entry below "
  "failed a check, and **no DOI is carried for any of them**.\n")
A("---\n")

A("## Group A. Right paper, wrong record. Needs only a DOI lookup\n")
A("Crossref found the correct work but returned a preprint, a review record, or the wrong "
  "year. Confirm the published version and record its DOI.\n")
A("| citation | what Crossref returned | what to do |")
A("|---|---|---|")
for r in sorted(G["A"], key=lambda x: x["citation_key"]):
    t = (r["crossref_title"] or "(nothing)")[:70]
    j = r["crossref_journal"] or ""
    A(f"| **{r['citation_key']}** | {t}{(' — ' + j) if j else ''} | {r['reason_if_flagged']} |".replace(" — ", ", "))
A("")

A("## Group B. Wrong record entirely. Needs the original publication checked\n")
A("Crossref returned something unrelated or nothing at all. Confirm authors, year, title, "
  "journal, volume and pages from the paper itself.\n")
A("| citation | what Crossref returned | why it failed |")
A("|---|---|---|")
for r in sorted(G["B"], key=lambda x: x["citation_key"]):
    t = (r["crossref_title"] or "(no record returned)")[:70]
    A(f"| **{r['citation_key']}** | {t} | {r['reason_if_flagged']} |")
A("")

A("## Group C. Secondary transcription. Needs checking against the original\n")
A("These were transcribed from the reference list of Swanepoel *et al.* (2018), which is where "
  "the Introduction draws them from. The wording is theirs, not the publishers'. **Not checked "
  "against the original publications.**\n")
A("| citation | as transcribed |")
A("|---|---|")
for r in sorted(G["C"], key=lambda x: x["citation_key"]):
    A(f"| **{r['citation_key']}** | {r['crossref_title']} |")
A("")

A("---\n")
A("## Verified, no action needed\n")
for r in sorted(ver, key=lambda x: x["citation_key"]):
    A(f"- **{r['citation_key']}.** {r['crossref_title']}"
      + (f" *{r['crossref_journal']}*." if r["crossref_journal"] else "")
      + (f" doi:{r['doi_if_verified']}" if r["doi_if_verified"] else ""))

io.open(OUT, "w", encoding="utf-8").write("\n".join(L).replace("—", ",").replace("–", "-") + "\n")
print(f"{len(flag)} flagged: A={len(G['A'])} DOI lookup, B={len(G['B'])} check original, "
      f"C={len(G['C'])} secondary transcription")
print(f"{len(ver)} verified")
print(f"-> {OUT}")
