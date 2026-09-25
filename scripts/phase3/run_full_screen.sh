#!/usr/bin/env bash
# 13a whole-table non-fungal screen against NCBI nt. Tiers run in risk order and each
# writes its own result file, so partial completion is still usable and auditable.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
S=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/full_screen
FMT="6 qseqid sacc pident length qcovs evalue bitscore staxids sscinames sskingdoms stitle"
for T in 1 2 3 4; do
  n=$(grep -c '^>' "$S/tier${T}.fasta")
  echo "=== tier $T ($n sequences) starting $(date -u +%H:%M:%S) ==="
  # split into chunks of 50 so a single failed submission cannot lose a whole tier
  rm -f "$S/tier${T}_hits.tsv"
  awk -v S="$S" -v T="$T" 'BEGIN{c=0;f=0}/^>/{if(c%50==0){f++};c++}{print > (S"/_t"T"_chunk"f".fa")}' "$S/tier${T}.fasta"
  for CH in "$S"/_t${T}_chunk*.fa; do
    for attempt in 1 2 3; do
      blastn -remote -db nt -query "$CH" -outfmt "$FMT" \
        -max_target_seqs 5 -evalue 1e-5 >> "$S/tier${T}_hits.tsv" 2>> "$S/tier${T}.err" && break
      echo "   retry $attempt on $(basename "$CH")"; sleep 20
    done
    echo "   $(basename "$CH") done, cumulative hits $(wc -l < "$S/tier${T}_hits.tsv" 2>/dev/null || echo 0), $(date -u +%H:%M:%S)"
    rm -f "$CH"
  done
  echo "=== tier $T complete: $(wc -l < "$S/tier${T}_hits.tsv" 2>/dev/null || echo 0) hit lines ==="
done
echo "FULL SCREEN COMPLETE $(date -u +%H:%M:%S)"
