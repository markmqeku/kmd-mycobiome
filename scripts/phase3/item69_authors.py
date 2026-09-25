#!/usr/bin/env python
"""69-ADD. Extract the author list, affiliations, corresponding author and any ORCID from
Article1_Bloemfontein_Baseline_Fungal_Communities.docx, exactly as recorded. Four copies exist
in files/ to files4/; all are compared so a divergence between them cannot go unnoticed."""
import hashlib, os, re
from docx import Document
K = r"<KMD_ROOT>"
NAME = "Article1_Bloemfontein_Baseline_Fungal_Communities.docx"
paths = [os.path.join(K, d, NAME) for d in ("files", "files2", "files3", "files4")]
paths = [p for p in paths if os.path.exists(p)]

print("=== copies found ===")
digests = {}
for p in paths:
    h = hashlib.sha256(open(p, "rb").read()).hexdigest()[:16]
    digests.setdefault(h, []).append(os.path.relpath(p, K))
    print(f"   {os.path.relpath(p, K):74} sha256:{h}")
print(f"   distinct versions: {len(digests)}"
      + ("  (all identical)" if len(digests) == 1 else "  <-- THEY DIFFER"))

d = Document(paths[0])
paras = [p.text.strip() for p in d.paragraphs]
print("\n=== first 30 non-empty paragraphs, verbatim ===")
n = 0
for t in paras:
    if not t: continue
    print(f"   [{n:>2}] {t[:200]}")
    n += 1
    if n >= 30: break

print("\n=== targeted extraction ===")
joined = "\n".join(paras)
for label, pat in (
    ("ORCID", r"ORCID[^\n]{0,120}|\b\d{4}-\d{4}-\d{4}-\d{3}[\dX]\b"),
    ("email / corresponding", r"[\w.\-]+@[\w.\-]+|[Cc]orrespond\w*[^\n]{0,140}"),
    ("affiliation markers", r"^\s*\d\s*[A-Z][^\n]{10,140}$"),
    ("University / Institute / ARC lines", r"^[^\n]*(University|Institute|Department|Faculty|ARC|Agricultural Research)[^\n]*$"),
):
    hits = re.findall(pat, joined, re.M)
    seen, out = set(), []
    for h in hits:
        s = " ".join(h.split())
        if s and s not in seen: seen.add(s); out.append(s)
    print(f"\n   {label}:")
    if out:
        for s in out[:12]: print(f"      {s[:170]}")
    else:
        print("      none found")
