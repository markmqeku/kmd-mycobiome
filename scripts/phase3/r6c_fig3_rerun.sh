#!/usr/bin/env bash
# KMD PROMPT 4, R6c. Figure 3 placement rerun with the seven green-algal ASVs excluded: the 71 retained
# lineage ASVs plus exactly the same 52 curated references and the unnamed match MK019123.1 used for the
# original tree (fig3_input.fasta), same MAFFT (--auto) and IQ-TREE settings (GTR+G, 1000 ultrafast
# bootstrap replicates, 4 threads). The original run's seed (276583, from fig3_tree.log) is fixed here so
# the rerun is reproducible. Outputs go to a new folder; the original run is left untouched.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
O=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026
S=$O/yeast_placement; D=$O/yeast_placement_71; mkdir -p "$D"
python - <<'PY'
import re
O="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
def rd(p):
    s={};k=None
    for l in open(p):
        l=l.rstrip("\n")
        if l.startswith(">"): k=l[1:]; s[k]=""
        elif k: s[k]+=l.strip()
    return s
inp=rd(f"{O}/yeast_placement/fig3_input.fasta")
kept={v.upper() for v in rd(f"{O}/clean/genbank_submission/Tremellomycetes_ASVs.fasta").values()}
out={h:s for h,s in inp.items() if not h.startswith("ASV_") or s.upper() in kept}
n_asv=sum(1 for h in out if h.startswith("ASV_")); n_ref=sum(1 for h in out if h.startswith("REF_")); n_un=sum(1 for h in out if h.startswith("UNNAMED_"))
assert n_asv==71 and n_ref==52 and n_un==1, (n_asv,n_ref,n_un)
with open(f"{O}/yeast_placement_71/fig3_71_input.fasta","w") as fh:
    for h,s in out.items(): fh.write(f">{h}\n{s}\n")
print(f"input: {n_asv} ASVs + {n_ref} references + {n_un} unnamed = {len(out)} tips; dropped {len(inp)-len(out)} green-algal ASVs")
PY
mafft --version 2>&1 | head -1
mafft --auto --quiet --thread 4 "$D/fig3_71_input.fasta" > "$D/fig3_71_aln.fasta" 2>"$D/fig3_71_mafft.log"
echo "aligned: $(grep -c '^>' "$D/fig3_71_aln.fasta") sequences"
cd "$D"
iqtree2 -s fig3_71_aln.fasta -m GTR+G -B 1000 -T 4 -seed 276583 --prefix fig3_71_tree -quiet -redo 2>fig3_71_iqtree.log
grep -E "Input data|Model of substitution|Log-likelihood of the tree" fig3_71_tree.iqtree
ls -la fig3_71_tree.contree fig3_71_tree.treefile
