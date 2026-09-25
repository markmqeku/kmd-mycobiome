#!/usr/bin/env python
"""Item 32. Cross-check the four draft sections for terminology, sample counts and
consistency of claims. READ ONLY."""
import io, re
K = "<KMD_ROOT>"
S = {n: io.open(f"{K}/DRAFT_{n}.md", encoding="utf-8").read()
     for n in ("INTRODUCTION", "METHODS", "RESULTS", "DISCUSSION", "TITLE_ABSTRACT")}

print("=== terminology ===")
TERMS = ["Tshwane", "Pretoria", "disease-free", "reference-site", "reference site",
         "healthy", "asymptomatic", "symptomatic", "malformed", "witches",
         "amplicon sequence variant", "ASV", "morphospecies"]
print(f"{'term':26} " + " ".join(f"{n[:5]:>6}" for n in S))
for t in TERMS:
    print(f"{t:26} " + " ".join(f"{S[n].lower().count(t.lower()):>6}" for n in S))

print("\n=== quantities shared across sections ===")
NUM = ["1,142", "2,517,177", "57,124", "199", "71", "78", "40 of 45", "84.3", "37.1",
       "1,328", "1,201", "127", "5 of 10", "7 of 10", "16", "36", "93.2"]
print(f"{'value':12} " + " ".join(f"{n[:5]:>6}" for n in S))
for v in NUM:
    print(f"{v:12} " + " ".join(f"{S[n].count(v):>6}" for n in S))

print("\n=== 'healthy' usage: Methods forbids it as a label ===")
for n, t in S.items():
    for m in re.finditer(r"\bhealthy\b", t, re.I):
        ctx = " ".join(t[max(0, m.start()-95):m.start()+95].split())
        print(f"   [{n}] ...{ctx}...")

print("\n=== deposition claims (records say nothing has been submitted) ===")
for n, t in S.items():
    for m in re.finditer(r"\b(deposited|deposit|submission|submitted)\b", t, re.I):
        ctx = " ".join(t[max(0, m.start()-105):m.start()+105].split())
        print(f"   [{n}] ...{ctx}...")

print("\n=== elapsed-time claims ===")
for n, t in S.items():
    for m in re.finditer(r"\b(a decade|ten years|\d+ years)\b", t, re.I):
        ctx = " ".join(t[max(0, m.start()-95):m.start()+95].split())
        print(f"   [{n}] ...{ctx}...")
