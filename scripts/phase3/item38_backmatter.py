#!/usr/bin/env python
"""Items 38 to 40. Append front and back matter, the reference list and the supplementary
manifest to the assembled manuscript. Nothing is invented: anything not on record is written
as an explicit INSERT rather than guessed."""
import csv, io, os
K = "<KMD_ROOT>"
O = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026"; C = f"{O}/clean"
M = f"{K}/MANUSCRIPT_KMD_mycobiome.md"
doc = io.open(M, encoding="utf-8").read()

refs = list(csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t"))
# Crossref-returned strings can carry en dashes; the manuscript requires none. Normalise the
# character only, in transcribed metadata fields. No wording is altered.
for r in refs:
    for k in ("crossref_title", "crossref_journal", "reason_if_flagged"):
        if r.get(k): r[k] = r[k].replace("–", "-").replace("—", "-")
ver = [r for r in refs if r["status"] == "VERIFIED"]
flag = [r for r in refs if r["status"] != "VERIFIED"]
man = list(csv.DictReader(open(f"{C}/figure_table_manifest.tsv"), delimiter="\t"))

L = []
A = L.append
A("\n---\n")
A("## Author contributions\n")
A("Recorded in CRediT taxonomy. [INSERT: confirm the contribution of each named author "
  "before submission. The allocation below is a proposal based on the project record and "
  "has not been confirmed by the authors.]\n")
A("- **M. Mqeku.** Investigation; Formal analysis; Data curation; Visualization; "
  "Writing, original draft.\n"
  "- **M. Gryzenhout.** Conceptualization; Supervision; Resources; Funding acquisition; "
  "Writing, review and editing.\n"
  "- [INSERT: remaining authors and their CRediT roles, including the contributors "
  "acknowledged on the 2017 conference poster if they are to be listed as authors.]\n")
A("## Funding\n")
A("This work was supported by the DST-NRF Centre of Excellence in Tree Health Biotechnology. "
  "[INSERT: confirm whether the 2017 conference poster's wording should also appear. That "
  "poster credited the Centre of Tree Health Biotechnology (CTHB), hosted by the Forestry "
  "and Agricultural Biotechnology Institute (FABI), University of Pretoria. Confirm whether "
  "both forms are required, and supply grant numbers and any further funding bodies "
  "supporting the sequencing or the reanalysis.]\n")
A("## Acknowledgements\n")
A("The authors thank the Next Generation Sequencing facility at the University of the Free "
  "State for generating the sequence data. [INSERT: confirm named individuals to be "
  "acknowledged, and confirm that the municipal authority at Christiana requires no "
  "acknowledgement for site access.]\n")
A("## Conflict of interest\n")
A("The authors declare that they have no known competing financial interests or personal "
  "relationships that could have appeared to influence the work reported in this paper.\n")
A("## Data availability\n")
A("Raw sequence data have been prepared for deposition in the NCBI Sequence Read Archive "
  "[INSERT accession on submission]. Representative sequences of the unresolved "
  "Tremellomycetes lineage have been prepared for deposition [INSERT accessions on "
  "submission]. **No sequence from this study has yet been deposited.** All processing "
  "parameters, software versions, database releases, per-sample read statistics, the "
  "complete ASV catalogue and the full non-fungal screen are provided as supplementary "
  "material. Analysis scripts are archived with the project record.\n")

A("\n---\n")
A("## References\n")
A(f"Compiled from all citations in the assembled text. {len(refs)} entries: "
  f"**{len(ver)} verified against Crossref, {len(flag)} requiring manual confirmation.** "
  "No DOI is carried for any entry that failed verification, and no bibliographic detail "
  "has been guessed.\n")
A("### Verified\n")
for r in sorted(ver, key=lambda x: x["citation_key"]):
    t = r["crossref_title"]; j = r["crossref_journal"]; y = r["crossref_year"]
    A(f"- **{r['citation_key']}.** {t}." + (f" *{j}*." if j else "") + (f" {y}." if y else "")
      + (f" doi:{r['doi_if_verified']}" if r["doi_if_verified"] else ""))
A("\n### Requiring manual confirmation\n")
A("Each entry below failed an automated year-plus-title match. The reason is given so it can "
  "be checked against the original publication. **A returned DOI is not evidence of a correct "
  "match**, and none is carried forward here.\n")
for r in sorted(flag, key=lambda x: x["citation_key"]):
    A(f"- **{r['citation_key']}.** {r['crossref_title'] or '(no record returned)'}"
      + (f" *{r['crossref_journal']}*." if r["crossref_journal"] else "")
      + f"\n  - *Why flagged:* {r['reason_if_flagged']}")

A("\n---\n")
A("## Supplementary material\n")
A(f"All supplementary items are built on the corrected **1,142-ASV** fungal table.\n")
A("| item | description | file | status |")
A("|---|---|---|---|")
for r in man:
    A(f"| {r['number']} | {r['title']} | `{os.path.basename(r['file'])}` | {r['status']} |")

io.open(M, "w", encoding="utf-8").write(doc.rstrip() + "\n" + "\n".join(L) + "\n")
print(f"appended front and back matter, {len(refs)} references, {len(man)} supplementary items")
print(f"   verified {len(ver)}, flagged {len(flag)}")
print(f"-> {M}")
