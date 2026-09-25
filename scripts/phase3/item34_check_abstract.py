#!/usr/bin/env python
"""Item 34. Verify DRAFT_TITLE_ABSTRACT.md. READ ONLY."""
import csv, io, re
K = "<KMD_ROOT>"
C = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
a = io.open(f"{K}/DRAFT_TITLE_ABSTRACT.md", encoding="utf-8").read()
res = io.open(f"{K}/DRAFT_RESULTS.md", encoding="utf-8").read()
fin = io.open(f"{K}/__reanalysis_2026-06/phase3/FINAL_NUMBERS_27-07-2026.md", encoding="utf-8").read()
bad = []

print("=== 34a. quantities against FINAL_NUMBERS ===")
CLAIMS = [("1,142", "corrected fungal ASVs"), ("2,517,177", "corrected fungal reads"),
          ("45", "samples sequenced"), ("16", "core analytical samples"),
          ("two to six", "per-group n"), ("71", "yeast lineage ASVs"),
          ("40 of 45", "samples holding the lineage"), ("85.1 to 85.5", "nearest named identity"),
          ("five years", "gap to reference release"), ("84.3", "Bloemfontein host percent"),
          ("16 percent", "reads to the agaric"), ("six", "failure modes")]
for v, what in CLAIMS:
    inA = v in a
    inF = v.replace(",", ",") in fin or v in res
    print(f"   {'ok  ' if inA else 'ABSENT':6} {v:14} {what:36} also in Results/FINAL_NUMBERS: {'yes' if inF else 'CHECK'}")
    if inA and not inF: bad.append(f"'{v}' in abstract but not located in Results/FINAL_NUMBERS")

print("\n=== 34b. statistical and causal language ===")
PAT = [(r"\bp\s*[=<>]\s*0?\.\d+", "p-value"), (r"\bsignifican\w*", "significance"),
       (r"\bPERMANOVA\b|\bANOVA\b|\bt-test\b", "test"),
       (r"\bcaus(e|ed|es|ing)\b", "causal verb"), (r"\bprove[sd]?\b|\bdemonstrat\w*", "proof"),
       (r"\bshow(s|ed)? that\b", "assertive claim"), (r"\bdue to\b|\bbecause of\b", "attribution")]
hits = 0
for pat, lab in PAT:
    for m in re.finditer(pat, a, re.I):
        ctx = " ".join(a[max(0, m.start()-100):m.start()+100].split())
        print(f"   [{lab}] ...{ctx}..."); hits += 1
if not hits: print("   none")

print("\n=== 34c. deposition wording ===")
for m in re.finditer(r"[^.]*\bdeposit\w*[^.]*\.", a):
    print(f"   ABSTRACT: {' '.join(m.group(0).split())}")
for nm, txt in (("RESULTS", res),
                ("DISCUSSION", io.open(f"{K}/DRAFT_DISCUSSION.md", encoding="utf-8").read()),
                ("METHODS", io.open(f"{K}/DRAFT_METHODS.md", encoding="utf-8").read())):
    for m in re.finditer(r"[^.]*\b(has been deposited|have been prepared for deposition|are deposited|never submitted)\b[^.]*\.", txt):
        print(f"   {nm}: {' '.join(m.group(0).split())[:190]}")
if re.search(r"\b(are|is|was|were|has been|have been)\s+deposited\b", a):
    bad.append("abstract asserts completed deposition")
    print("   FAIL: abstract asserts completed deposition")
else:
    print("   ok: abstract makes no completed-deposition claim")

print("\n=== 34d. abstract claims traceable to Results ===")
SUPPORT = {
 "1,142 ASVs / 2,517,177 reads": "1,142 ASVs and 2,517,177 reads",
 "16 core samples": "16 samples",
 "71 variants in 40 of 45": "40 of the 45 samples",
 "85.1 to 85.5 percent": "85.1 to 85.5",
 "Curvibasidium direction": "Curvibasidium",
 "diversity not consistently elevated": "richness",
 "ordination by site": "separated samples by site",
 "two genera not recovered": "Mycosphaerella",
 "84.3 percent host": "84.3",
 "16 percent to an agaric": "16 percent",
 "host screen could not detect other families": "unrelated family",
}
for claim, needle in SUPPORT.items():
    ok = needle.lower() in res.lower()
    print(f"   {'ok  ' if ok else 'CHECK':6} {claim:42} <- Results contains '{needle}'")
    if not ok: bad.append(f"abstract claim not located in Results: {claim}")

print("\n=== 34e. dashes ===")
em, en = a.count("—"), a.count("–")
print(f"   em {em}, en {en}  {'ok' if em == en == 0 else 'FAIL'}")
if em or en: bad.append("dashes present")

print("\n=== SUMMARY ===")
print("no automated issues" if not bad else f"{len(bad)} issue(s):")
for b in bad: print(f"   - {b}")
