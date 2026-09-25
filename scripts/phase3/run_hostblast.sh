#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
D=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen
makeblastdb -in "$D/hostref/host_ref.fasta" -dbtype nucl -out "$D/hostref/hostdb" >/dev/null
blastn -db "$D/hostref/hostdb" -query "$D/nontarget_query_all.fasta" \
  -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle" \
  -max_target_seqs 3 -evalue 1e-5 -num_threads 4 -out "$D/blast_host.tsv"
echo "host-db hit lines: $(wc -l < "$D/blast_host.tsv")"
python <KMD_ROOT>/__reanalysis_2026-06/phase3/summarise_host.py
