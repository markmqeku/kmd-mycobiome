import csv
B = "<KMD_ROOT>/__reanalysis_2026-06"
cat = f"{B}/phase3/outputs_27-07-2026/asv_catalogue/ASV_catalogue.tsv"
tab = f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv"
meta = f"{B}/phase3/outputs_27-07-2026/metadata_v2_27-07-2026.tsv"

rows = list(csv.DictReader(open(cat), delimiter="\t"))
tot_reads = sum(int(r["total_reads"]) for r in rows)

print("=== TOP 15 ASVs BY TOTAL READS ===")
for r in sorted(rows, key=lambda r: -int(r["total_reads"]))[:15]:
    print(f"  {r['display_label']:52} reads={int(r['total_reads']):>9,} ({100*int(r['total_reads'])/tot_reads:5.2f}%)  prev={r['prevalence_n_samples']}/45  conf={r['classifier_confidence'][:6]}")

hom = [r for r in rows if "Homophron" in r["display_label"]]
hr = sum(int(r["total_reads"]) for r in hom)
print(f"\n=== Homophron: {len(hom)} ASVs | {hr:,} reads = {100*hr/tot_reads:.2f}% of dataset ===")
print(f"  max prevalence of any Homophron ASV: {max(int(r['prevalence_n_samples']) for r in hom)}/45")
print(f"  mean sequence length: {sum(int(r['length_bp']) for r in hom)/len(hom):.0f} bp")
print(f"  confidence values (first 5): {[r['classifier_confidence'][:5] for r in hom[:5]]}")

# which samples carry Homophron reads
lines = [l.rstrip('\n') for l in open(tab) if l.strip() and not l.startswith("# Constructed")]
hdr = lines[0].lstrip('#').split('\t'); samples = [h for h in hdr[1:] if h.strip()]
idx = {r["ASV_ID_md5"] for r in hom}
persample = {s: 0 for s in samples}
for l in lines[1:]:
    p = l.split('\t')
    if p[0] in idx:
        for s, v in zip(samples, p[1:1+len(samples)]):
            persample[s] += int(float(v))
md = {r[0]: r for r in csv.reader(open(meta), delimiter="\t")}
top = sorted(persample.items(), key=lambda x: -x[1])[:10]
print("\n  top samples by Homophron reads:")
for s, v in top:
    loc = md[s][1] if s in md else "?"
    print(f"    {s:6} {loc:13} {v:>9,}")
nz = sum(1 for v in persample.values() if v > 0)
print(f"  samples containing any Homophron: {nz}/45")
