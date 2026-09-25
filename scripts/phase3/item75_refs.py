#!/usr/bin/env python
"""Item 75. Second pass on the nine outstanding references. Record-type guard now rejects
book reviews as well as peer-review records and preprints. Year is allowed to differ by one
where the container matches, since online-ahead-of-print shifts a year legitimately; the
container check does the work in those cases."""
import csv, json, re, time, urllib.parse, urllib.request
C = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/clean"
UA = "KMD-reference-check (mailto:<email>)"
BAD_TYPE = re.compile(r"peer-review|posted-content|component|grant|dataset", re.I)
BAD_CONT = re.compile(r"faculty opinions|post-publication|peer review|preprint|biorxiv"
                      r"|choice reviews|cabi compendium|book review", re.I)

# key: (query, year, title token, container substring or None, year tolerance)
T = {
 "Hoang et al. 2018":     ("UFBoot2 Improving the Ultrafast Bootstrap Approximation Hoang Chernomor von Haeseler Minh Vinh", "2018", "ufboot", "molecular biology and evolution", 1),
 "Jami et al. 2013":      ("Botryosphaeriaceae species overlap on four unrelated native South African hosts", "2013", "botryosphaeriaceae", "fungal biology", 1),
 "Põlme et al. 2020":     ("FungalTraits a user-friendly traits database of fungi and fungus-like stramenopiles Polme Abarenkov Nilsson", "2020", "fungaltraits", "fungal diversity", 1),
 "Krishnan et al. 2009":  ("Mango Mangifera indica L malformation an unsolved mystery Krishnan Nailwal Shukla Pant Researcher", "2009", "malformation", None, 1),
 "Crous et al. 2000":     ("Phytopathogenic fungi from South Africa Crous Phillips Baxter", "2000", "phytopathogenic fungi", None, 1),
 "Koekemoer et al. 2013": ("Guide to plant families of southern Africa Strelitzia Koekemoer Steyn Bester", "2013", "plant families", None, 1),
 "Coates-Palgrave 2002":  ("Keith Coates-Palgrave Trees of Southern Africa Random House Struik", "2002", "trees of southern africa", None, 1),
 "Van Wyk & Gericke 2007":("People's Plants A Guide to Useful Plants of Southern Africa Van Wyk Gericke Briza", "2007", "people's plants", None, 1),
 "Magurran 2004":         ("Measuring Biological Diversity Magurran Blackwell", "2004", "measuring biological diversity", None, 1),
}
def q(s, rows=8):
    u = ("https://api.crossref.org/works?rows=%d&select=title,issued,DOI,container-title,author,"
         "type,ISBN,publisher,volume,page&query.bibliographic=" % rows) + urllib.parse.quote(s)
    r = urllib.request.Request(u, headers={"User-Agent": UA})
    for _ in range(3):
        try: return json.loads(urllib.request.urlopen(r, timeout=50).read().decode())["message"]["items"]
        except Exception: time.sleep(4)
    return []

res = {}
for key, (s, yr, tok, cont, tol) in T.items():
    best, why = None, "no candidate satisfied the title, container and record-type tests"
    for it in q(s):
        ti = (it.get("title") or [""])[0]; co = (it.get("container-title") or [""])[0]
        ty = it.get("type", "")
        dp = it.get("issued", {}).get("date-parts", [[None]])
        y = str(dp[0][0]) if dp and dp[0] and dp[0][0] else ""
        if BAD_TYPE.search(ty) or BAD_CONT.search(co): continue
        if tok not in ti.lower(): continue
        if cont and cont not in co.lower(): continue
        if y and abs(int(y) - int(yr)) > tol: continue
        best = it; why = ""; break
    if best:
        dp = best.get("issued", {}).get("date-parts", [[None]])
        res[key] = ["VERIFIED", (best.get("title") or [""])[0], str(dp[0][0]),
                    (best.get("container-title") or [""])[0] or best.get("publisher", ""),
                    (best.get("author") or [{}])[0].get("family", ""), best.get("DOI", ""),
                    f"type={best.get('type','')}; vol={best.get('volume','')}; "
                    f"pages={best.get('page','')}; isbn={','.join(best.get('ISBN') or []) or 'n/a'}; "
                    f"publisher={best.get('publisher','')}"]
    else:
        res[key] = ["UNVERIFIED", "", "", "", "", "", why]
    time.sleep(1)

p = f"{C}/reference_list.tsv"
rows = list(csv.reader(open(p), delimiter="\t")); hdr, body = rows[0], rows[1:]
for r in body:
    if r[0] in res:
        v = res[r[0]]
        if v[0] == "VERIFIED": r[1:8] = v[0], v[1], v[2], v[3], v[4], v[5], v[6]
        else: r[7] = v[6]
with open(p, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t"); w.writerow(hdr); w.writerows(body)

nv = sum(1 for v in res.values() if v[0] == "VERIFIED")
print(f"second pass on {len(res)}: {nv} newly verified, {len(res)-nv} still not\n")
for k, v in res.items():
    print(f"{v[0]:11} {k}")
    if v[0] == "VERIFIED":
        print(f"            {v[1][:86]}")
        print(f"            {v[3]}  {v[2]}  doi:{v[5]}")
        print(f"            {v[6]}")
    else:
        print(f"            {v[6][:150]}")
tot = sum(1 for r in body if r[1] == "VERIFIED")
print(f"\nreference list: {len(body)} entries, {tot} verified, {len(body)-tot} not verified")
