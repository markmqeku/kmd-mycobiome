#!/usr/bin/env bash
# G-3b/c: align the 78 yeast-lineage ASVs with Tremellomycetes references, build a tree, report placement.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
B=<KMD_ROOT>/__reanalysis_2026-06
O=$B/phase3/outputs_27-07-2026
D=$O/yeast_placement

# 78 ASV representatives
python - <<'PY'
import csv
B="<KMD_ROOT>/__reanalysis_2026-06"
D=f"{B}/phase3/outputs_27-07-2026/yeast_placement"
rows=[r for r in csv.reader(open(f"{B}/phase3/outputs_27-07-2026/nontarget_screen/query_abundance.tsv"),delimiter="\t") if r and r[0]!="ASV_ID_md5"]
unk={r[0]:int(r[1]) for r in rows if r[2]=="1"}
seqs,sid,buf={},None,[]
for line in open(f"{B}/phase2/merged2/rep_exp/dna-sequences.fasta"):
    line=line.rstrip("\n")
    if line.startswith(">"):
        if sid: seqs[sid]="".join(buf)
        sid,buf=line[1:].split()[0],[]
    else: buf.append(line.strip())
if sid: seqs[sid]="".join(buf)
with open(f"{D}/yeast_asvs.fasta","w") as fh:
    for a in sorted(unk,key=lambda x:-unk[x]):
        fh.write(f">ASV_{a[:8]}_reads{unk[a]}\n{seqs[a]}\n")
print(f"yeast ASVs written: {len(unk)}")
PY

# nearest references by identity
makeblastdb -in "$D/tremello_ref.fasta" -dbtype nucl -out "$D/tremellodb" >/dev/null
blastn -db "$D/tremellodb" -query "$D/yeast_asvs.fasta" \
  -outfmt "6 qseqid sseqid pident length qcovs evalue" -max_target_seqs 5 -evalue 1e-5 -num_threads 4 \
  -out "$D/yeast_vs_tremello.tsv"
echo "blast hit lines: $(wc -l < "$D/yeast_vs_tremello.tsv")"

# align + tree
cat "$D/yeast_asvs.fasta" "$D/tremello_ref.fasta" > "$D/combined.fasta"
mafft --auto --quiet --thread 4 "$D/combined.fasta" > "$D/combined_aln.fasta" 2>"$D/mafft.log"
echo "aligned: $(grep -c '^>' "$D/combined_aln.fasta") sequences"
fasttree -nt -gtr -quiet "$D/combined_aln.fasta" > "$D/yeast_placement.tree" 2>"$D/fasttree.log"
echo "tree written: $D/yeast_placement.tree"
python <KMD_ROOT>/__reanalysis_2026-06/phase3/G3_analyse_tree.py
