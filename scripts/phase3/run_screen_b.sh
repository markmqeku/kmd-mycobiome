#!/usr/bin/env bash
# Second stream for the 13a screen: chunks 4-6, written to its own hit file so the two
# streams cannot interleave writes. Two concurrent NCBI jobs only, which stays within
# acceptable remote BLAST use.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
S=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/full_screen
FMT="6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames sskingdoms stitle"
rm -f "$S/screen_hits_b.tsv"
for i in 4 5 6; do
  CH="$S/_b${i}.fa"
  [ -f "$CH" ] || { echo "missing $CH, skipping"; continue; }
  for attempt in 1 2 3; do
    blastn -remote -db nt -query "$CH" -outfmt "$FMT" -max_target_seqs 5 -evalue 1e-5 \
      >> "$S/screen_hits_b.tsv" 2>> "$S/screen_b.err" && break
    echo "   retry $attempt on $(basename "$CH")"; sleep 30
  done
  echo "$(basename "$CH") done, cumulative $(wc -l < "$S/screen_hits_b.tsv" 2>/dev/null || echo 0) lines, $(date -u +%H:%M:%S)"
  rm -f "$CH"
done
echo "STREAM B COMPLETE $(date -u +%H:%M:%S)"
