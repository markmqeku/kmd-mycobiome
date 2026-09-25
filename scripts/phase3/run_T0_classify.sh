#!/usr/bin/env bash
# T-0 (part 2): classify merged rep-seqs against UNITE 10.0; retain 8.2 as documented sensitivity comparison.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
R=<KMD_ROOT>/__reanalysis_2026-06
OUT=$R/phase3/outputs_27-07-2026
mkdir -p "$OUT/taxonomy_unite10"

echo "### provenance of the fetched UNITE 10.0 artifact (release / version / citation)"
qiime tools peek "$OUT/unite_current/unite10_seqs.qza"
python - <<'PY'
import zipfile,glob,re
p=glob.glob("<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/unite_current/unite10_seqs.qza")[0]
z=zipfile.ZipFile(p)
acts=[n for n in z.namelist() if n.endswith("action.yaml")]
txt="\n".join(z.read(a).decode("utf-8","replace") for a in acts)
for kw in ["version","doi","citation","url","cluster_id","taxon_group","singletons"]:
    for line in txt.splitlines():
        if kw in line.lower() and len(line.strip())<300:
            print("   ",line.strip())
PY

echo; echo "### classify vs UNITE 10.0"
qiime feature-classifier classify-consensus-vsearch \
  --i-query "$R/phase2/merged2/rep.qza" \
  --i-reference-reads "$OUT/unite_current/unite10_seqs.qza" \
  --i-reference-taxonomy "$OUT/unite_current/unite10_tax.qza" \
  --p-threads 4 --p-maxaccepts 10 --p-perc-identity 0.8 \
  --o-classification "$OUT/taxonomy_unite10/taxonomy.qza" \
  --o-search-results "$OUT/taxonomy_unite10/search.qza" 2>"$OUT/taxonomy_unite10/classify.log"
qiime tools export --input-path "$OUT/taxonomy_unite10/taxonomy.qza" --output-path "$OUT/taxonomy_unite10/exp"
echo "CLASSIFY 10.0 DONE"

echo; echo "### sensitivity comparison: UNITE 8.2 vs 10.0 genus-level assignment"
python - <<'PY'
import csv
def load(p):
    d={}
    with open(p) as f:
        rd=csv.reader(f,delimiter="\t"); next(rd)
        for r in rd:
            if r and not r[0].startswith("#"): d[r[0]]=r[1]
    return d
def genus(t):
    if not t or t.lower().startswith("unassigned"): return ""
    for part in t.split(";"):
        part=part.strip()
        if part.startswith("g__"):
            g=part[3:].strip()
            if g.lower() in ("","unidentified","unclassified"): return ""
            return g
    return ""
old=load("<KMD_ROOT>/__reanalysis_2026-06/phase2/taxonomy/exp/taxonomy.tsv")
new=load("<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/taxonomy_unite10/exp/taxonomy.tsv")
common=set(old)&set(new)
same=gain=lost=changed=0; examples=[]
for a in common:
    go,gn=genus(old[a]),genus(new[a])
    if go==gn: same+=1
    elif go and not gn: lost+=1
    elif gn and not go: gain+=1
    else:
        changed+=1
        if len(examples)<15: examples.append((a[:8],go,gn))
print(f"ASVs compared: {len(common)}")
print(f"  genus identical            : {same} ({100*same/len(common):.1f}%)")
print(f"  genus gained (8.2 none->10): {gain}")
print(f"  genus lost   (8.2->none)   : {lost}")
print(f"  genus CHANGED (different)  : {changed}")
print(f"  TOTAL genus-level changes  : {gain+lost+changed} ({100*(gain+lost+changed)/len(common):.1f}%)")
if examples:
    print("  examples (asv, 8.2, 10.0):")
    for e in examples: print("   ",e)
og=sum(1 for a in common if genus(old[a])); ng=sum(1 for a in common if genus(new[a]))
print(f"  ASVs with a genus: 8.2={og} ({100*og/len(common):.1f}%)  10.0={ng} ({100*ng/len(common):.1f}%)")
PY
echo "T0 COMPLETE"
