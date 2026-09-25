#!/usr/bin/env python
"""G-3e: check accepted status of candidate genera against Index Fungorum / MycoBank web services.
No name is asserted from NCBI's legacy label. Reports failure honestly if the services are unreachable."""
import urllib.request, urllib.parse, json

GENERA = ["Filobasidium", "Cryptococcus", "Naganishia", "Papiliotrema", "Vishniacozyma", "Curvibasidium"]

def get(url, timeout=40):
    req = urllib.request.Request(url, headers={"User-Agent": "KMD-reanalysis/1.0"})
    with urllib.request.urlopen(req, timeout=timeout) as r:
        return r.read().decode("utf-8", "replace")

print("=== Index Fungorum / MycoBank lookup ===")
for g in GENERA:
    done = False
    # Index Fungorum REST
    try:
        u = "http://www.indexfungorum.org/ixfwebservice/fungus.asmx/NameSearch?SearchText=" \
            + urllib.parse.quote(g) + "&AnywhereInText=false&MaxNumber=3"
        x = get(u)
        if "<" in x and g.lower() in x.lower():
            import re
            names = re.findall(r"<NAME_x0020_OF_x0020_FUNGUS>(.*?)</", x)[:3]
            auth  = re.findall(r"<AUTHORS>(.*?)</", x)[:3]
            cur   = re.findall(r"<CURRENT_x0020_NAME>(.*?)</", x)[:3]
            print(f"  {g:16} IF names={names} authors={auth} current={cur}")
            done = True
    except Exception as e:
        print(f"  {g:16} IndexFungorum unreachable: {type(e).__name__}")
    if not done:
        try:
            u = "https://api.gbif.org/v1/species/match?name=" + urllib.parse.quote(g)
            d = json.loads(get(u))
            print(f"  {g:16} GBIF backbone: status={d.get('status')} rank={d.get('rank')} "
                  f"accepted={d.get('canonicalName')} class={d.get('class')} order={d.get('order')} "
                  f"family={d.get('family')} conf={d.get('confidence')}")
        except Exception as e:
            print(f"  {g:16} GBIF unreachable too: {type(e).__name__}")
