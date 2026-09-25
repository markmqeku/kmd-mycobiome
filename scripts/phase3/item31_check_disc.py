#!/usr/bin/env python
"""Item 31. Verify DRAFT_DISCUSSION.md. READ ONLY; never writes to that file."""
import csv, io, re, zipfile
from collections import Counter
K = "<KMD_ROOT>"
O = f"{K}/__reanalysis_2026-06/phase3/outputs_27-07-2026"; C = f"{O}/clean"
d = io.open(f"{K}/DRAFT_DISCUSSION.md", encoding="utf-8").read()
issues = []

print("=== 31a. quantities against FINAL_NUMBERS ===")
EXPECT = {"71":"yeast lineage ASVs","40 of 45":"samples holding the lineage",
          "85.1 to 85.5":"nearest named identity","276":"alignment bp","199":"Bloemfontein taxa",
          "84.3":"Bloemfontein host percent","36":"FUNGuild-assigned yeast ASVs",
          "93.2":"missed Searsia identity","1,201":"Swanepoel malformed isolates",
          "127":"Swanepoel healthy isolates","81":"Homophron identity"}
for v, what in EXPECT.items():
    print(f"   {'ok  ' if v in d else 'ABSENT'}  {v:12} {what}")

# Curvibasidium claim, recomputed
meta = {}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd = csv.reader(fh, delimiter="\t"); h = next(rd); next(rd); di = {c: i for i, c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]] = {c: r[di[c]] for c in h}
lines = [l.rstrip("\n") for l in open(f"{C}/feature-table-clean.tsv")
         if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip("#").split("\t"); cols = [x for x in hdr[1:] if x.strip()]
per = {s: {} for s in cols}
for l in lines[1:]:
    q = l.split("\t")
    for s, v in zip(cols, q[1:1+len(cols)]):
        v = int(float(v))
        if v: per[s][q[0]] = v
tax = {}
for r in csv.reader(open(f"{O}/taxonomy_unite10/exp/taxonomy.tsv"), delimiter="\t"):
    if r and r[0] != "Feature ID" and not r[0].startswith("#"): tax[r[0]] = r[1]
pid = {}
z = zipfile.ZipFile(f"{O}/taxonomy_unite10/search.qza")
fn = [n for n in z.namelist() if "/data/" in n and n.endswith((".tsv", ".blast6"))][0]
for line in z.read(fn).decode("utf-8", "replace").splitlines():
    p = line.split("\t")
    if len(p) < 3: continue
    try: v = float(p[2])
    except ValueError: continue
    if p[0] not in pid or v > pid[p[0]]: pid[p[0]] = v
def gen(a):
    for part in tax.get(a, "").split(";"):
        if part.strip().startswith("g__"): return part.strip()[3:]
    return ""
curv = {a for a in tax if gen(a) == "Curvibasidium" and pid.get(a, 0) >= 95}
def pc(s):
    t = sum(per[s].values()); return 100*sum(v for a, v in per[s].items() if a in curv)/t if t else 0
ca = [s for s in cols if meta.get(s) and meta[s]["location"] == "Christiana" and meta[s]["condition_std"] == "asymptomatic"]
cs = [s for s in cols if meta.get(s) and meta[s]["location"] == "Christiana" and meta[s]["condition_std"] == "symptomatic"]
va, vs = sorted(pc(s) for s in ca), sorted(pc(s) for s in cs)
print(f"\n   Curvibasidium at Christiana")
print(f"      asymptomatic: {[f'{x:.3f}' for x in va]}  max {max(va):.3f}")
print(f"      symptomatic : {[f'{x:.3f}' for x in vs]}")
gt1 = sum(1 for x in vs if x > 1.0)
print(f"      claim 'no asymptomatic sample exceeded 0.35 percent': "
      f"{'TRUE' if max(va) <= 0.35 else 'FALSE'} (max {max(va):.3f})")
print(f"      claim 'five of six symptomatic samples exceeded 1 percent': "
      f"{'TRUE' if gt1 == 5 else 'FALSE'} (actually {gt1} of 6 strictly exceed 1.0)")
if gt1 != 5: issues.append(f"5.2 Curvibasidium: '{gt1} of 6' exceed 1 percent, text says five of six")

print("\n=== 31c. statistical and causal language ===")
PAT = [(r"\bp\s*[=<>]\s*0?\.\d+", "p-value"), (r"\bsignifican\w*", "significance"),
       (r"\bPERMANOVA\b|\bANOVA\b|\bt-test\b|\bcorrelat\w*", "test or correlation"),
       (r"\bcaus(e|ed|es|ing|al)\b", "causal"), (r"\bproves?\b|\bdemonstrat(e|ed|es)\b", "proof"),
       (r"\bassociated with\b", "association"), (r"\bgenuine\b|\brepresents a\b", "assertion")]
for pat, lab in PAT:
    for m in re.finditer(pat, d, re.I):
        ctx = " ".join(d[max(0, m.start()-110):m.start()+110].split())
        print(f"   [{lab}] ...{ctx}...")

print("\n=== 31e. dashes ===")
em, en = d.count("—"), d.count("–")
print(f"   em dashes: {em}   en dashes: {en}   {'ok' if em == en == 0 else 'FAIL'}")
if em or en: issues.append("dashes present")

print("\n=== 31d. citations ===")
ref = {r["citation_key"]: r["status"] for r in csv.DictReader(open(f"{C}/reference_list.tsv"), delimiter="\t")}
def norm(s): return re.sub(r"[^a-z0-9]", "", s.lower())
keys = {norm(k): k for k in ref}
NAME = r"[A-Z][A-Za-zÀ-ſ'\-]+"
cit = set()
for m in re.finditer(NAME + r"(?:\s+and\s+" + NAME + r")?(?:\s+et\s+al\.)?,?\s*\(?((?:\d{4}[a-z]?))\)?", d):
    nm = " ".join(m.group(0).split()).rstrip("(),")
    nm = re.sub(r"[,(]?\s*\d{4}[a-z]?\)?$", "", nm).strip()
    if not nm or nm.split()[0].lower() in ("the","in","and","section","figure","table","one","two","three","several","given","because","their","its"): continue
    cit.add((nm, m.group(1)))
miss = []
for nm, yr in sorted(cit):
    sur = nm.split()[0]
    hit = next((keys[k] for k in keys if norm(sur) in k and yr in k), None)
    if hit: print(f"   ok      {nm} {yr} -> {hit} [{ref[hit]}]")
    else: miss.append(f"{nm} {yr}"); print(f"   MISSING {nm} {yr}")
if miss: issues.append(f"citations absent from list: {miss}")

print("\n=== SUMMARY ===")
print("no automated issues" if not issues else f"{len(issues)} issue(s):")
for i in issues: print(f"   - {i}")
