import csv, json, urllib.request, urllib.parse
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
accs = set()
for r in csv.reader(open(f"{O}/nontarget_screen/blast_priority.tsv"), delimiter="\t"):
    if len(r) >= 8 and "Cryptococcus" in r[7]:
        accs.add(r[1])
accs = sorted(accs)[:5]
accs.append("PP844358.1")   # nearest curated Tremellomycetes reference (by bitscore)
print("accessions to date-check:", accs)
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
try:
    x = urllib.request.urlopen(
        f"{E}/esummary.fcgi?db=nuccore&id={urllib.parse.quote(','.join(accs))}&retmode=json", timeout=60
    ).read().decode()
    d = json.loads(x)
    for k, v in d.get("result", {}).items():
        if k == "uids": continue
        print(f"  {v.get('accessionversion','?'):14} created={v.get('createdate','?')} "
              f"{v.get('title','')[:60]}")
except Exception as e:
    print("esummary failed:", type(e).__name__, e)
print("\nUNITE 10.0 released 2024-04-04.")
