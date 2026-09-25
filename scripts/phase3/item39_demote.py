#!/usr/bin/env python
"""Item 39 correction. Two entries passed the year-plus-title test but Crossref had returned
a record ABOUT the paper rather than the paper: a Faculty Opinions post-publication review and
a PeerJ peer-review record. Title and year match in both cases, which is exactly why the test
missed them. Demote both, strip their DOIs, and add a record-type guard so this class of error
is caught rather than re-derived by eye."""
import csv, re
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
p = f"{C}/reference_list.tsv"
rows = list(csv.reader(open(p), delimiter="\t"))
hdr, body = rows[0], rows[1:]
i_status, i_j, i_doi, i_why = 1, 4, 6, 7

BAD_JOURNAL = re.compile(r"faculty opinions|post-publication|peer review|preprint|biorxiv", re.I)
BAD_DOI = re.compile(r"10\.3410/|/reviews?/|v\d+\.\d+/reviews", re.I)
demoted = []
for r in body:
    if r[i_status] != "VERIFIED": continue
    j, doi = r[i_j], r[i_doi]
    why = None
    if BAD_JOURNAL.search(j or ""): why = f"Crossref returned a review record, not the article: '{j}'"
    elif BAD_DOI.search(doi or ""): why = f"the DOI resolves to a review or preprint record, not the article: {doi}"
    elif not (j or "").strip(): why = "Crossref returned no journal for this record, so the match cannot be confirmed"
    if why:
        r[i_status] = "FLAGGED"; r[i_doi] = ""; r[i_why] = why
        demoted.append((r[0], why))
with open(p, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)
v = sum(1 for r in body if r[1] == "VERIFIED")
print(f"demoted {len(demoted)} entries that had passed on year and title alone:")
for k, why in demoted: print(f"   {k}\n      {why}")
print(f"\nreference list now: {len(body)} entries, {v} verified, {len(body)-v} flagged.")
print("A returned DOI is still not evidence of a correct match. Record type must be checked too.")
