#!/usr/bin/env bash
# C-4 taxonomy: classify-consensus-vsearch (RAM-safe; no sklearn training) vs UNITE 8.2 (nearest-available v8; v8.0 not fetchable).
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2
mkdir -p "$WD/taxonomy"
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$WD/merged2/rep.qza" \
  --i-reference-reads "$WD/unite/unite82_seqs.qza" \
  --i-reference-taxonomy "$WD/unite/unite82_tax.qza" \
  --p-threads 4 --p-maxaccepts 10 --p-perc-identity 0.8 \
  --o-classification "$WD/taxonomy/taxonomy.qza" \
  --o-search-results "$WD/taxonomy/search.qza" --verbose 2>"$WD/taxonomy/classify.log"
qiime tools export --input-path "$WD/taxonomy/taxonomy.qza" --output-path "$WD/taxonomy/exp"
echo "TAXONOMY DONE"
# summary: how many ASVs classified to each rank
python - <<PY
import csv
rows=[r for r in csv.reader(open("$WD/taxonomy/exp/taxonomy.tsv"),delimiter="\t")][2:]
tot=len(rows)
def depth(tax):
    parts=[p for p in tax.split(";") if p.strip() and not p.strip().endswith("__")]
    return len(parts)
levels=["kingdom","phylum","class","order","family","genus","species"]
from collections import Counter
c=Counter(min(depth(r[1]),7) for r in rows)
print(f"total ASVs classified: {tot}")
for i,lv in enumerate(levels,1):
    n=sum(v for k,v in c.items() if k>=i); print(f"  to {lv:8}: {n} ({100*n/tot:.0f}%)")
unassigned=sum(1 for r in rows if r[1].strip() in ("Unassigned","") )
print(f"  Unassigned: {unassigned}")
PY
