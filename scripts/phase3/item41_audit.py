#!/usr/bin/env python
"""Item 41. Final audit of the assembled manuscript. Reports only; changes nothing."""
import csv, io, re
K = "<KMD_ROOT>"
C = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
d = io.open(f"{K}/MANUSCRIPT_KMD_mycobiome.md", encoding="utf-8").read()
fin = io.open(f"{K}/__reanalysis_2026-06/phase3/FINAL_NUMBERS_27-07-2026.md", encoding="utf-8").read()
fails = []

print("=== 1. numbers against FINAL_NUMBERS ===")
MUST = {"1,142":"fungal ASVs","2,517,177":"fungal reads","57,124":"rarefaction depth",
        "199":"Bloemfontein taxa","71":"yeast lineage","40 of 45":"lineage occurrence",
        "84.3":"host percent","37.1":"host share","1,322":"total ASVs","149":"host ASVs",
        "4,010,617":"total reads","93.2":"missed Searsia","85.1 to 85.5":"nearest named"}
for v, w in MUST.items():
    ok = v in d
    print(f"   {'ok  ' if ok else 'ABSENT':7} {v:12} {w}")
STALE = {"57,141":"superseded depth","192 taxa":"superseded inventory",
         "542 of 1,173":"superseded guild denominator","33.1 percent":"unreproducible ordination",
         "143.0 ASVs":"superseded sharing","1,152":"intermediate ASV count",
         # superseded Hill q0 group means, depth 57,141. These sat undetected in the
         # Results table until the cross-reference pass; they belong in this list.
         "| 80.1 |":"superseded Hill q0, Christiana asymptomatic",
         "| 107.6 |":"superseded Hill q0, Christiana symptomatic",
         "| 116.7 |":"superseded Hill q0, Pretoria symptomatic"}
print("   superseded values that must not appear unlabelled:")
for v, w in STALE.items():
    hits = [m.start() for m in re.finditer(re.escape(v), d)]
    bad = [i for i in hits if not re.search(r"supersed|before that correction|previously|then included|first delimited",
                                            d[max(0,i-400):i+400], re.I)]
    print(f"      {'ok  ' if not bad else 'FAIL':6} {v:16} {w}" + (f"  ({len(bad)} unlabelled)" if bad else ""))
    if bad: fails.append(f"unlabelled superseded value {v}")

print("\n=== 2. taxon italicisation ===")
TAXA = ["Searsia lancea","Curvibasidium","Homophron spadiceum","Amaranthus","Fusarium",
        "Alternaria alternata","Filobasidium","Cryptococcus","Didymella","Mycosphaerella",
        "Filobasidiella","Tremellomycetes","Dysphania melanocarpa","Vigna unguiculata",
        "F. temperatum","Mangifera indica","Muribasidiospora indica"]
GENUS_NOT_ITAL = {"Tremellomycetes"}   # class name, correctly not italicised
for t in TAXA:
    plain = len(re.findall(r"(?<![*\w])" + re.escape(t) + r"(?![*\w])", d))
    ital = len(re.findall(r"\*" + re.escape(t) + r"\*", d))
    if t in GENUS_NOT_ITAL:
        print(f"   {'ok':6} {t:24} class name, not italicised ({plain} plain)")
        continue
    status = "ok" if plain == 0 else "CHECK"
    if plain: fails.append(f"{t} appears unitalicised {plain}x")
    print(f"   {status:6} {t:24} italic {ital:>3}, plain {plain:>3}")

print("\n=== 3. significance and causal language ===")
PAT = [(r"\bp\s*[=<>]\s*0?\.\d+","p-value"),(r"\bsignifican\w*","significance"),
       (r"\bPERMANOVA\b|\bANOVA\b|\bt-test\b","test name")]
n = 0
for pat, lab in PAT:
    for m in re.finditer(pat, d, re.I):
        ctx = " ".join(d[max(0,m.start()-110):m.start()+110].split())
        print(f"   [{lab}] ...{ctx}..."); n += 1; fails.append(f"{lab} present")
if not n: print("   none")

print("\n=== 4. condition terminology ===")
for m in re.finditer(r"\bhealthy\b", d, re.I):
    ctx = " ".join(d[max(0,m.start()-120):m.start()+120].split())
    ok = "Swanepoel" in ctx or "isolates" in ctx or "avoided" in ctx or "culture" in ctx
    print(f"   {'ok  ' if ok else 'CHECK':6} ...{ctx[:150]}...")
    if not ok: fails.append("'healthy' used outside the Swanepoel context")
for t in ("symptomatic","asymptomatic","reference-site"):
    print(f"   '{t}' used {len(re.findall(chr(92)+'b'+t+chr(92)+'b', d))} times")

print("\n=== 5. sample sizes beside group values ===")
for m in re.finditer(r"\bn\s*=\s*\d+", d):
    pass
print(f"   explicit 'n =' statements: {len(re.findall(r'n\s*=\s*\d+', d))}")
print(f"   'two to six samples' framing present: {'yes' if 'two to six' in d else 'NO'}")
print(f"   group table with n column present: {'yes' if re.search(r'\|\s*n\s*\|', d) else 'NO'}")

print("\n=== 6. dashes ===")
em, en = d.count("—"), d.count("–")
print(f"   em {em}, en {en}   {'ok' if em == en == 0 else 'FAIL'}")
if em or en: fails.append("dashes present")

print("\n=== 7. citations resolve ===")
ref = {r["citation_key"] for r in csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t")}
def norm(s): return re.sub(r"[^a-z0-9]","",s.lower())
keys = {norm(k): k for k in ref}
body = d.split("## References")[0]
cited, miss = set(), []
NAME = r"[A-Z][A-Za-zÀ-ſ'\-]+"
# Methods and Results write "White *et al.* (1990)"; the asterisks must not break detection
for m in re.finditer(NAME + r"(?:\s+and\s+" + NAME + r")?(?:\s+\*?et\s+al\.\*?)?,?\s*\(?(\d{4})\)?", body):
    nm = re.sub(r"[,(]?\s*\d{4}\)?$","",m.group(0)).strip()
    if not nm or nm.split()[0].lower() in ("the","in","and","section","figure","table","one","two",
        "three","several","given","because","their","its","from","of","at","under","after","with"): continue
    hit = next((keys[k] for k in keys if norm(nm.split()[0]) in k and m.group(1) in k), None)
    if hit: cited.add(hit)
    else: miss.append(f"{nm} {m.group(1)}")
print(f"   in-text citations resolving to the list: {len(cited)}")
print(f"   in-text citations NOT in the list: {len(set(miss))}")
for x in sorted(set(miss)): print(f"      MISSING {x}")
uncited = sorted(ref - cited)
print(f"   list entries never cited in the body: {len(uncited)}")
for x in uncited: print(f"      UNCITED {x}")

print("\n=== SUMMARY ===")
print("AUDIT CLEAN" if not fails else f"{len(fails)} item(s) to review:")
for f in sorted(set(fails)): print(f"   - {f}")
