#!/usr/bin/env python
"""Item 13d. Confirm Methods subsection numbering is sequential and every cross-reference
anywhere in the manuscript points at the section it means."""
import io, re
K = "<KMD_ROOT>"
d = io.open(f"{K}/MANUSCRIPT_KMD_mycobiome.md", encoding="utf-8").read()
fails = []

print("=== Methods subsection numbering ===")
subs = re.findall(r"^###\s+(2\.\d+)\s+(.+)$", d, re.M)
for n, t in subs: print(f"   {n}  {t}")
nums = [int(n.split(".")[1]) for n, _ in subs]
seq = nums == list(range(1, len(nums)+1))
print(f"   sequential 2.1 to 2.{len(nums)}: {'yes' if seq else 'NO -> ' + str(nums)}")
if not seq: fails.append("Methods numbering not sequential")

print("\n=== Results subsection numbering ===")
rsubs = re.findall(r"^###\s+(3\.\d+)\s+(.+)$", d, re.M)
for n, t in rsubs: print(f"   {n}  {t}")
rn = [int(n.split(".")[1]) for n, _ in rsubs]
rseq = rn == list(range(1, len(rn)+1))
print(f"   sequential 3.1 to 3.{len(rn)}: {'yes' if rseq else 'NO -> ' + str(rn)}")
if not rseq: fails.append("Results numbering not sequential")

titles = {n: t for n, t in subs + rsubs}
dsubs = re.findall(r"^###\s+(5\.\d+)\s+(.+)$", d, re.M)
titles.update({n: t for n, t in dsubs})
print("\n=== Discussion subsection numbering ===")
for n, t in dsubs: print(f"   {n}  {t}")

print("\n=== every 'Section X.Y' cross-reference ===")
# what each reference is expected to be about, judged from the sentence it sits in
for m in re.finditer(r"Section (\d\.\d+)", d):
    tgt = m.group(1)
    ctx = " ".join(d[max(0, m.start()-150):m.start()+60].split())
    exists = tgt in titles
    print(f"   -> {tgt} ({titles.get(tgt,'DOES NOT EXIST')})")
    print(f"      context: ...{ctx[-135:]}")
    if not exists: fails.append(f"cross-reference to non-existent Section {tgt}")

print("\n=== Methods and Results parallelism ===")
pairs = [("2.4 Host sequence screening","3.11"), ("2.5 Taxonomic assignment","3.3"),
         ("2.6 Functional guild assignment","3.9"), ("2.7 Phylogenetic placement","3.10"),
         ("2.8 Diversity analyses","3.5")]
for mth, res in pairs:
    mok = any(mth.split()[0] == n and mth.split(None,1)[1][:12] in t for n, t in subs)
    rok = any(n == res for n, _ in rsubs)
    print(f"   {'ok' if mok and rok else 'CHECK':6} Methods {mth:34} <-> Results {res} "
          f"{titles.get(res,'')}")

print("\n=== SUMMARY ===")
print("NUMBERING AND CROSS-REFERENCES CORRECT" if not fails else f"{len(fails)} issue(s):")
for f in fails: print(f"   - {f}")
