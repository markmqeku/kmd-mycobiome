#!/usr/bin/env bash
# Item 7 b,c: per-neighbourhood MAFFT alignment + IQ-TREE with 1000 ultrafast bootstrap,
# plus nearest-reference identity by BLAST (bitscore-ranked). Never one combined tree.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
P=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/placement
cd "$P"
for nb in Didymellaceae Microbotryomycetes Mortierellaceae Dothideomycetes Coniochaetales Dothideaceae Saccotheciaceae Sclerotiniaceae; do
  [ -f "${nb}_asvs.fasta" ] || { echo "SKIP $nb (no ASV file)"; continue; }
  [ -f "${nb}_refs.fasta" ] || { echo "SKIP $nb (no refs)"; continue; }
  echo "=== $nb ==="
  # nearest reference by bitscore
  makeblastdb -in "${nb}_refs.fasta" -dbtype nucl -out "${nb}_db" >/dev/null 2>&1
  blastn -db "${nb}_db" -query "${nb}_asvs.fasta" \
    -outfmt "6 qseqid sseqid pident length qcovhsp evalue bitscore" \
    -max_target_seqs 20 -evalue 1e-5 -num_threads 4 -out "${nb}_blast.tsv" 2>/dev/null
  echo "  blast lines: $(wc -l < "${nb}_blast.tsv")"
  # align + tree
  cat "${nb}_asvs.fasta" "${nb}_refs.fasta" > "${nb}_combined.fasta"
  mafft --auto --quiet --thread 4 "${nb}_combined.fasta" > "${nb}_aln.fasta" 2>/dev/null
  ntip=$(grep -c '^>' "${nb}_aln.fasta")
  iqtree2 -s "${nb}_aln.fasta" -m GTR+G -B 1000 -T 4 --prefix "${nb}_tree" -quiet -redo >/dev/null 2>&1
  if [ -f "${nb}_tree.contree" ]; then echo "  tree OK ($ntip tips)"; else echo "  TREE FAILED"; fi
done
echo "ITEM7 TREES DONE"
