#!/usr/bin/env bash
# T-0b: re-classify vs UNITE 10.0 at a defensible identity threshold (0.97) for ITS species/genus assignment.
# The 0.80 run is RETAINED as a documented sensitivity comparison. Nothing overwritten.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
R=<KMD_ROOT>/__reanalysis_2026-06
OUT=$R/phase3/outputs_27-07-2026
mkdir -p "$OUT/taxonomy_unite10_id97"

qiime feature-classifier classify-consensus-vsearch \
  --i-query "$R/phase2/merged2/rep.qza" \
  --i-reference-reads "$OUT/unite_current/unite10_seqs.qza" \
  --i-reference-taxonomy "$OUT/unite_current/unite10_tax.qza" \
  --p-threads 4 --p-maxaccepts 10 --p-perc-identity 0.97 --p-min-consensus 0.51 \
  --o-classification "$OUT/taxonomy_unite10_id97/taxonomy.qza" \
  --o-search-results "$OUT/taxonomy_unite10_id97/search.qza" 2>"$OUT/taxonomy_unite10_id97/classify.log"
qiime tools export --input-path "$OUT/taxonomy_unite10_id97/taxonomy.qza" --output-path "$OUT/taxonomy_unite10_id97/exp"
echo "CLASSIFY id97 DONE"

python - <<'PY'
import csv
B="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
RANKS=[("k__","kingdom"),("p__","phylum"),("c__","class"),("o__","order"),("f__","family"),("g__","genus"),("s__","species")]
PH={"unidentified","unclassified","incertae_sedis","unknown",""}
def load(p):
    d={}
    for r in csv.reader(open(p),delimiter="\t"):
        if r and not r[0].startswith("#") and r[0]!="Feature ID": d[r[0]]=r[1]
    return d
def ranks(t):
    out={r:"" for _,r in RANKS}
    if not t or t.lower().startswith("unassigned"): return out
    for part in t.split(";"):
        part=part.strip()
        for pre,rk in RANKS:
            if part.startswith(pre):
                v=part[len(pre):].strip()
                if v.lower() in PH or v.lower().startswith("unidentified"): v=""
                out[rk]=v
    return out
a=load(f"{B}/taxonomy_unite10/exp/taxonomy.tsv")
b=load(f"{B}/taxonomy_unite10_id97/exp/taxonomy.tsv")
common=sorted(set(a)&set(b))
print(f"\n=== RANK RESOLUTION: identity 0.80 (permissive) vs 0.97 (defensible), UNITE 10.0, n={len(common)} ===")
print(f"{'rank':9} {'id=0.80':>12} {'id=0.97':>12}")
for _,rk in RANKS:
    na=sum(1 for x in common if ranks(a[x])[rk]); nb=sum(1 for x in common if ranks(b[x])[rk])
    print(f"{rk:9} {na:>6} ({100*na/len(common):4.1f}%) {nb:>6} ({100*nb/len(common):4.1f}%)")
ha=sum(1 for x in common if "Homophron" in a[x]); hb=sum(1 for x in common if "Homophron" in b[x])
print(f"\nHomophron-labelled ASVs: id0.80={ha}   id0.97={hb}")
PY
echo "T0b COMPLETE"
