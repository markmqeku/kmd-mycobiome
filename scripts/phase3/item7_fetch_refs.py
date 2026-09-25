#!/usr/bin/env python
"""Item 7a: retrieve curated references per neighbourhood, preferring ex-type and
voucher-backed accessions. Records accession, name and type status for every reference kept."""
import urllib.request, urllib.parse, time, os, re, json

P = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/placement"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NEIGH = ["Didymellaceae", "Microbotryomycetes", "Mortierellaceae", "Dothideomycetes",
         "Coniochaetales", "Dothideaceae", "Saccotheciaceae", "Sclerotiniaceae"]

def get(u, t=90):
    return urllib.request.urlopen(
        urllib.request.Request(u, headers={"User-Agent": "KMD-reanalysis/1.0"}), timeout=t
    ).read().decode("utf-8", "replace")

def esearch(term, retmax):
    x = get(f"{E}/esearch.fcgi?db=nuccore&retmax={retmax}&term=" + urllib.parse.quote(term))
    return [i.split("<")[0] for i in x.split("<Id>")[1:]]

rows = []
for nb in NEIGH:
    # ex-type / voucher-backed first, then general ITS, deduplicated
    q_type = (f'{nb}[Organism] AND (internal transcribed spacer[All Fields]) '
              f'AND (type[All Fields] OR "ex-type"[All Fields] OR voucher[All Fields])')
    q_any = f'{nb}[Organism] AND internal transcribed spacer[All Fields]'
    ids = []
    for q, n in ((q_type, 30), (q_any, 40)):
        try:
            for i in esearch(q, n):
                if i not in ids: ids.append(i)
        except Exception as e:
            print(f"  {nb}: esearch failed {type(e).__name__}")
        time.sleep(0.4)
    ids = ids[:55]
    if not ids:
        print(f"  {nb}: NO REFERENCES FOUND"); continue
    try:
        fa = get(f"{E}/efetch.fcgi?db=nuccore&id={','.join(ids)}&rettype=fasta&retmode=text")
        sm = json.loads(get(f"{E}/esummary.fcgi?db=nuccore&id={','.join(ids)}&retmode=json"))
    except Exception as e:
        print(f"  {nb}: efetch failed {type(e).__name__}"); continue
    titles = {v.get("accessionversion", k): v.get("title", "")
              for k, v in sm.get("result", {}).items() if k != "uids"}
    kept, cur, buf = {}, None, []
    for line in fa.splitlines():
        if line.startswith(">"):
            if cur: kept[cur] = "".join(buf)
            cur, buf = line[1:], []
        else: buf.append(line.strip())
    if cur: kept[cur] = "".join(buf)
    out = f"{P}/{nb}_refs.fasta"
    n_written = 0
    with open(out, "w") as fh:
        for hdr, seq in kept.items():
            if len(seq) < 200 or len(seq) > 2000: continue
            acc = hdr.split()[0]
            title = titles.get(acc, hdr)
            is_type = bool(re.search(r"\btype\b|ex-type|isotype|holotype|neotype", title, re.I))
            is_vouch = bool(re.search(r"voucher|strain|isolate|culture", title, re.I))
            name = " ".join(title.split()[1:3]) if len(title.split()) > 2 else title[:40]
            safe = re.sub(r"[^A-Za-z0-9._]", "", f"{name}_{acc}")[:48]
            fh.write(f">REF_{safe}\n{seq}\n"); n_written += 1
            rows.append([nb, acc, name, "type" if is_type else ("voucher/strain" if is_vouch else "none"), len(seq)])
    print(f"  {nb}: {n_written} references written")
    time.sleep(0.5)

with open(f"{P}/reference_provenance.tsv", "w", newline="", encoding="utf-8") as fh:
    import csv
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["neighbourhood", "accession", "name_from_title", "type_status", "length_bp"])
    w.writerows(rows)
print(f"\nreference provenance recorded for {len(rows)} accessions -> {P}/reference_provenance.tsv")
n_type = sum(1 for r in rows if r[3] == "type")
print(f"  of these, {n_type} carry type wording and {sum(1 for r in rows if r[3]=='voucher/strain')} are voucher or strain backed")
