#!/usr/bin/env python
"""G-3a: fetch reference ITS sequences spanning Tremellomycetes segregate genera.
No name is assigned from these; they are only used to SHOW placement on a tree."""
import urllib.parse, urllib.request, time, os
D = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/yeast_placement"
os.makedirs(D, exist_ok=True)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

def get(url):
    req = urllib.request.Request(url, headers={"User-Agent": "KMD-reanalysis/1.0"})
    with urllib.request.urlopen(req, timeout=90) as r:
        return r.read().decode("utf-8", "replace")

# Tremellomycetes genera, including the segregates of Cryptococcus sensu lato
GENERA = ["Cryptococcus", "Naganishia", "Papiliotrema", "Vishniacozyma", "Filobasidium",
          "Kwoniella", "Hannaella", "Saitozyma", "Solicoccozyma", "Bullera",
          "Dioszegia", "Tremella", "Curvibasidium"]
out = os.path.join(D, "tremello_ref.fasta")
open(out, "w").close()
tally = {}
for g in GENERA:
    q = f'{g}[Organism] AND (internal transcribed spacer[All Fields] OR ITS1[All Fields] OR ITS2[All Fields])'
    try:
        xml = get(f"{E}/esearch.fcgi?db=nuccore&retmax=25&term=" + urllib.parse.quote(q))
    except Exception as e:
        print(f"  {g}: esearch failed {e}"); continue
    ids = [x.split("<")[0] for x in xml.split("<Id>")[1:]]
    if not ids:
        print(f"  {g}: 0"); tally[g] = 0; continue
    try:
        fa = get(f"{E}/efetch.fcgi?db=nuccore&id={','.join(ids)}&rettype=fasta&retmode=text")
    except Exception as e:
        print(f"  {g}: efetch failed {e}"); continue
    # tag headers with the queried genus so tree tips are readable
    lines = []
    for ln in fa.splitlines():
        lines.append(f">REF_{g}_{ln[1:40].split()[0]}" if ln.startswith(">") else ln)
    with open(out, "a") as fh: fh.write("\n".join(lines) + "\n")
    tally[g] = fa.count(">")
    print(f"  {g}: {tally[g]}")
    time.sleep(0.4)
print(f"\ntotal reference sequences: {sum(tally.values())} -> {out}")
