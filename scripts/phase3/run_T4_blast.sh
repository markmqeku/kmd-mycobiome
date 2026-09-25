#!/usr/bin/env bash
# T-4a: remote blastn of the non-target candidate ASVs against NCBI nt.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
D=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen

echo "### remote blastn: priority set (covers ~59% of dataset reads)"
blastn -remote -db nt \
  -query "$D/nontarget_query_priority.fasta" \
  -out "$D/blast_priority.tsv" \
  -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle" \
  -max_target_seqs 5 -evalue 1e-5 2> "$D/blast_priority.err"
echo "exit=$?  hits: $(wc -l < "$D/blast_priority.tsv" 2>/dev/null || echo 0)"
head -3 "$D/blast_priority.err" 2>/dev/null
echo "T4_BLAST DONE"
