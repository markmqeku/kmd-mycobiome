#!/usr/bin/env python
"""29c. Add the Introduction's citations that are absent from the consolidated reference list.
Full citations are taken from the reference list of Swanepoel et al. (2018), which is where the
Introduction draws them from. They are recorded as SECONDARY, meaning the wording comes from
that paper's bibliography and has NOT been checked against the original publication. None is
marked verified, and none carries a DOI."""
import csv, os
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
NEW = [
 ("Coates-Palgrave 2002",
  "Coates-Palgrave, M., 2002. Keith Coates-Palgrave Trees of Southern Africa. 4th Impression, "
  "3rd Edition. Random House Struik, Cape Town.", "2002"),
 ("Crous et al. 2000",
  "Crous, P.W., Phillips, A.J., Baxter, A.P., 2000. Phytopathogenic Fungi from South Africa. "
  "University of Stellenbosch Printers, Department of Plant Pathology Press.", "2000"),
 ("Jami et al. 2013",
  "Jami, F., Slippers, B., Wingfield, M.J., Gryzenhout, M., 2013. Botryosphaeriaceae species "
  "overlap on four unrelated, native South African hosts.", "2013"),
 ("Koekemoer et al. 2013",
  "Koekemoer, M., Steyn, H.M., Bester, S.P., 2013. Guide to Plant Families of Southern Africa. "
  "Strelitzia 31. South African National Biodiversity Institute, Pretoria.", "2013"),
 ("Krishnan et al. 2009",
  "Krishnan, A.G., Nailwal, T.K., Shukla, A., Pant, R.C., 2009. Mango (Mangifera indica L.) "
  "malformation, an unsolved mystery. Researcher 1, 20 to 36.", "2009"),
 ("Lange et al. 2012",
  "Lange, C.A., Kotte, K., Smit, M., Van Deventer, P., Van Rensburg, L., 2012. Effects of "
  "different soil ameliorants on karee trees (Searsia lancea).", "2012"),
 ("Moffett 2007",
  "Moffett, R.O., 2007. Name changes in the Old World Rhus and recognition of Searsia "
  "(Anacardiaceae). Bothalia 37, 165 to 175.", "2007"),
 ("Van Wyk & Gericke 2007",
  "Van Wyk, B.E., Gericke, N., 2007. People's Plants: A Guide to Useful Plants of Southern "
  "Africa. Briza Publications, Pretoria.", "2007"),
]
p = f"{C}/reference_list.tsv"
rows = list(csv.reader(open(p), delimiter="\t"))
hdr, body = rows[0], rows[1:]
have = {r[0] for r in body}
added = []
for key, full, yr in NEW:
    if key in have: continue
    # columns: citation_key, status, crossref_title, crossref_year, crossref_journal,
    #          crossref_first_author, doi_if_verified, reason_if_flagged
    body.append([key, "FLAGGED", full, yr, "", "", "",
                 "SECONDARY: transcribed from the reference list of Swanepoel et al. (2018); "
                 "not checked against the original publication"])
    added.append(key)
with open(p, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)
print(f"added {len(added)} references, all flagged for verification:")
for k in added: print(f"   {k}")
tot = len(body); ver = sum(1 for r in body if r[1] == "VERIFIED")
print(f"\nconsolidated list now holds {tot} references: {ver} verified, {tot-ver} flagged.")
print(f"-> {p}")
