#!/usr/bin/env python
"""KMD PROMPT 5, D1 and D2: verify every ontology term in the EBI Ontology Lookup Service (OLS4) and every taxid in
NCBI Taxonomy (E-utilities) before it is written to a submission file. Writes p5_verified_terms.json."""
import json, time, urllib.parse, urllib.request
from pathlib import Path
try:
    import truststore; truststore.inject_into_ssl()
except ImportError:
    pass

OUT = Path(__file__).resolve().parent / "p5_verified_terms.json"
TERMS = {"ENVO:01000177": "grassland biome", "ENVO:01000178": "savanna biome", "ENVO:00000467": "university campus",
         "ENVO:01000447": "roadside", "ENVO:00000078": "farm", "ENVO:01001121": "plant matter",
         "PO:0025034": "leaf", "PO:0025073": "branch", "PO:0009049": "inflorescence", "PO:0009010": "seed",
         "PO:0009006": "shoot system"}


def get(u):
    req = urllib.request.Request(u, headers={"User-Agent": "KMD-verify/1.0", "Accept": "application/json"})
    return json.loads(urllib.request.urlopen(req, timeout=60).read())


res = {"ols": {}, "taxonomy": {}}
for cid, label in TERMS.items():
    onto = cid.split(":")[0].lower()
    iri = f"http://purl.obolibrary.org/obo/{cid.replace(':', '_')}"
    u = f"https://www.ebi.ac.uk/ols4/api/ontologies/{onto}/terms/" + urllib.parse.quote(urllib.parse.quote(iri, safe=""), safe="")
    d = get(u)
    ok = d.get("label") == label and not d.get("is_obsolete")
    res["ols"][cid] = {"expected": label, "label": d.get("label"), "obsolete": d.get("is_obsolete"),
                       "definition": (d.get("description") or [""])[0], "url": u, "ok": ok}
    print(f"{'OK ' if ok else 'BAD'} {cid} {d.get('label')!r} obsolete={d.get('is_obsolete')}")
    time.sleep(0.5)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
for name in ("Searsia lancea", "uncultured Tremellomycetes", "uncultured fungus", "plant metagenome"):
    ids = get(E + "esearch.fcgi?db=taxonomy&retmode=json&term=" + urllib.parse.quote(f'"{name}"[Scientific Name]'))["esearchresult"]["idlist"]
    s = get(E + "esummary.fcgi?db=taxonomy&retmode=json&id=" + ",".join(ids))["result"] if ids else {}
    rec = [{"taxid": i, "name": s[i].get("scientificname"), "rank": s[i].get("rank"), "division": s[i].get("division")} for i in ids]
    res["taxonomy"][name] = rec
    print(name, rec)
    time.sleep(0.5)
json.dump(res, open(OUT, "w", encoding="utf-8"), indent=1, ensure_ascii=False)
