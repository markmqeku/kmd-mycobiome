#!/usr/bin/env python
"""Item 43. Confirm every named software tool carries a citation, and that the reference list
and the in-text citations match each other exactly."""
import csv, io, re
K = "<KMD_ROOT>"
C = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
d = io.open(f"{K}/MANUSCRIPT_KMD_mycobiome.md", encoding="utf-8").read()
body = d.split("## References")[0]
fails = []

print("=== 43a. software tools and their citations ===")
TOOLS = {"QIIME2":"Bolyen","cutadapt":"Martin","DADA2":"Callahan","VSEARCH":"Rognes",
         "UNITE":"Nilsson","MAFFT":"Katoh","IQ-TREE":"Minh","ultrafast bootstrap":"Hoang",
         "iNEXT":"Hsieh","hillR":"Li","FungalTraits":"lme","FUNGuild":"Nguyen",
         "ITS3":"White","Hill numbers":"Chao"}
for tool, author in TOOLS.items():
    m = re.search(re.escape(tool), body, re.I)
    if not m:
        print(f"   {'n/a':8} {tool:22} not named in the text"); continue
    # look for the author within 260 characters after the tool is named
    win = body[m.start():m.start()+260]
    ok = re.search(re.escape(author), win)
    print(f"   {'ok' if ok else 'MISSING':8} {tool:22} cited as {author}")
    if not ok: fails.append(f"{tool} named without a citation")

print("\n=== 43b. list against text ===")
keys = [r["citation_key"] for r in csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t")]
uncited = []
for k in keys:
    sur = re.split(r"[ &]", k)[0]
    if not re.search(r"\b" + re.escape(sur), body): uncited.append(k)
print(f"   reference entries: {len(keys)}")
print(f"   entries never cited in the text: {len(uncited)}")
for k in uncited: print(f"      UNCITED {k}")
if uncited: fails.append(f"{len(uncited)} uncited entries")

surnames = {re.split(r"[ &]", k)[0].lower() for k in keys}
intext = set()
for m in re.finditer(r"([A-Z][A-Za-zÀ-ſ'\-]+)(?:\s+and\s+[A-Z][A-Za-z'\-]+)?\s*\*?et\s+al\.\*?,?\s*\(?\d{4}", body):
    intext.add(m.group(1).lower())
for m in re.finditer(r"\(([A-Z][A-Za-zÀ-ſ'\-]+)(?:\s+and\s+([A-Z][A-Za-z'\-]+))?,\s*\d{4}", body):
    intext.add(m.group(1).lower())
missing = sorted(n for n in intext if n not in surnames
                 and n not in {"figure","table","section","supplementary"})
print(f"   in-text first authors not in the list: {len(missing)}")
for n in missing: print(f"      NOT IN LIST {n}")
if missing: fails.append(f"{len(missing)} in-text citations absent from the list")

print("\n=== 43c. counts ===")
rows = list(csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t"))
v = sum(1 for r in rows if r["status"] == "VERIFIED")
print(f"   {len(rows)} references: {v} verified, {len(rows)-v} flagged for manual confirmation")

print("\n=== SUMMARY ===")
print("CITATIONS COMPLETE" if not fails else f"{len(fails)} issue(s):")
for f in fails: print(f"   - {f}")
