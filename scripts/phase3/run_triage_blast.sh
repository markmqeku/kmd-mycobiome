#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
P=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/placement
blastn -remote -db nt -query "$P/triage_query.fasta" -out "$P/triage_blast.tsv" \
  -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle" \
  -max_target_seqs 10 -evalue 1e-5 2> "$P/triage_blast.err"
echo "triage hits: $(wc -l < "$P/triage_blast.tsv" 2>/dev/null || echo 0)"
