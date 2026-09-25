#!/usr/bin/env python
"""Item 63. Read the exact UNITE release from QIIME2 artifact provenance and the RESCRIPt
fetch log. Provenance is the record of what was actually downloaded; a web page or a search
summary is not. No DOI is supplied unless it is present in the artifact itself."""
import os, re, zipfile
B = "<KMD_ROOT>/__reanalysis_2026-06/phase3"
O = f"{B}/outputs_27-07-2026"

print("=== 63a. RESCRIPt action recorded in the artifact provenance ===")
for q in ("unite_current/unite10_seqs.qza", "unite_current/unite10_tax.qza"):
    p = os.path.join(O, q)
    z = zipfile.ZipFile(p)
    act = [n for n in z.namelist() if n.endswith("action/action.yaml")][0]
    txt = z.read(act).decode("utf-8", "replace")
    print(f"\n--- {q} ---")
    for line in txt.splitlines():
        s = line.strip()
        if not s or s.startswith("#"): continue
        if re.match(r"^(action|plugin|type|version|parameters|-\s*\w+:)", s) or ":" in s:
            print("   " + line.rstrip())

print("\n\n=== DOI search across every provenance file in both artifacts ===")
found = []
for q in ("unite_current/unite10_seqs.qza", "unite_current/unite10_tax.qza",
          "taxonomy_unite10/search.qza"):
    z = zipfile.ZipFile(os.path.join(O, q))
    for n in z.namelist():
        if not (n.endswith(".yaml") or n.endswith(".txt") or n.endswith(".json")): continue
        t = z.read(n).decode("utf-8", "replace")
        for m in re.finditer(r"(10\.\d{4,9}/[^\s\"',;)\]]+|doi[^\n]{0,80})", t, re.I):
            found.append((q, n.split("/")[-1], m.group(0)[:110]))
if found:
    for q, n, d in found[:40]: print(f"   {q:34} {n:22} {d}")
else:
    print("   NO DOI of any kind appears anywhere in the provenance of any artifact.")

print("\n\n=== RESCRIPt fetch logs on disk ===")
logs = []
for root, _, files in os.walk(B):
    for f in files:
        if re.search(r"(unite|rescript|G4|G5)", f, re.I) and f.endswith((".log", ".txt", ".sh")):
            logs.append(os.path.join(root, f))
for l in sorted(set(logs))[:12]:
    print(f"\n--- {os.path.relpath(l, B)} ---")
    txt = open(l, encoding="utf-8", errors="replace").read()
    for line in txt.splitlines():
        if re.search(r"unite|version|doi|10\.15156|singleton|dynamic|eukaryote|fungi|--p-", line, re.I):
            print("   " + line.strip()[:150])
