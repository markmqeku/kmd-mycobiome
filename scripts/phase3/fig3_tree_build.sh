#!/usr/bin/env bash
# Figure 3 tree: 78 Tremellomycetes ASVs + curated named references + the unnamed nearest match
# MK019123.1, aligned with MAFFT, IQ-TREE with 1000 ultrafast bootstrap replicates.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
O=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026
D=$O/yeast_placement; F=$O/figures; mkdir -p "$F"

# fetch the unnamed nearest match
python - <<'PY'
import urllib.request
u="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=MK019123.1&rettype=fasta&retmode=text"
fa=urllib.request.urlopen(u,timeout=60).read().decode()
lines=fa.splitlines()
seq="".join(l for l in lines[1:] if l.strip())
p="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/yeast_placement/MK019123.fasta"
open(p,"w").write(">UNNAMED_Cryptococcus_sp_OTU796_MK019123.1\n"+seq+"\n")
print("MK019123.1 length:",len(seq))
PY

# thin the reference set to <=4 per genus for a legible figure, keep all 78 ASVs
python - <<'PY'
import collections
D="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/yeast_placement"
def rd(p):
    s={};k=None;b=[]
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"):
            if k: s[k]="".join(b)
            k=l[1:];b=[]
        else: b.append(l)
    if k: s[k]="".join(b)
    return s
refs=rd(f"{D}/tremello_ref.fasta"); asvs=rd(f"{D}/yeast_asvs.fasta"); unk=rd(f"{D}/MK019123.fasta")
cnt=collections.Counter(); keep={}
for h,s in refs.items():
    g=h.split("_")[1] if h.startswith("REF_") else "?"
    if cnt[g]<4 and len(s)>200:
        keep[h]=s; cnt[g]+=1
with open(f"{D}/fig3_input.fasta","w") as fh:
    for h,s in {**asvs,**keep,**unk}.items(): fh.write(f">{h}\n{s}\n")
print(f"fig3 input: {len(asvs)} ASVs + {len(keep)} refs + {len(unk)} unnamed = {len(asvs)+len(keep)+len(unk)} tips")
PY

mafft --auto --quiet --thread 4 "$D/fig3_input.fasta" > "$D/fig3_aln.fasta" 2>"$D/fig3_mafft.log"
echo "aligned: $(grep -c '^>' "$D/fig3_aln.fasta")"
cd "$D"
iqtree2 -s fig3_aln.fasta -m GTR+G -B 1000 -T 4 --prefix fig3_tree -quiet -redo 2>fig3_iqtree.log
echo "IQ-TREE done"; ls -la fig3_tree.contree fig3_tree.treefile 2>/dev/null
