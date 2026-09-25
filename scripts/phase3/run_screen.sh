#!/usr/bin/env bash
# 13a NCBI nt screen of the 277 ASVs not securely fungal by UNITE alignment.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
S=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/full_screen
FMT="6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames sskingdoms stitle"
rm -f "$S/screen_hits.tsv" "$S/screen.err"
awk -v S="$S" 'BEGIN{c=0;f=0}/^>/{if(c%50==0){f++};c++}{print > (S"/_chunk"f".fa")}' "$S/screen_set.fasta"
for CH in "$S"/_chunk*.fa; do
  n=$(grep -c '^>' "$CH")
  for attempt in 1 2 3; do
    blastn -remote -db nt -query "$CH" -outfmt "$FMT" -max_target_seqs 5 -evalue 1e-5 \
      >> "$S/screen_hits.tsv" 2>> "$S/screen.err" && break
    echo "   retry $attempt on $(basename "$CH")"; sleep 30
  done
  echo "$(basename "$CH") ($n seqs) done, cumulative hit lines $(wc -l < "$S/screen_hits.tsv" 2>/dev/null || echo 0), $(date -u +%H:%M:%S)"
  rm -f "$CH"
done
echo "SCREEN COMPLETE $(date -u +%H:%M:%S)"
