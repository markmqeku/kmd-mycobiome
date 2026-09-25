#!/usr/bin/env bash
# T-1 + T-2: display labels + ASV catalogue workbook, built on the UNITE 10.0 taxonomy (the release actually used).
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
R=<KMD_ROOT>/__reanalysis_2026-06
OUT=$R/phase3/outputs_27-07-2026

# feature table -> tsv
mkdir -p "$OUT/table_exp"
qiime tools export --input-path "$R/phase2/merged2/table.qza" --output-path "$OUT/table_exp" >/dev/null 2>&1
biom convert -i "$OUT/table_exp/feature-table.biom" -o "$OUT/table_exp/feature-table.tsv" --to-tsv
# masked alignment used for the tree -> fasta
mkdir -p "$OUT/aln_exp"
qiime tools export --input-path "$R/phase2/merged2/masked.qza" --output-path "$OUT/aln_exp" >/dev/null 2>&1
ALN="$OUT/aln_exp/aligned-dna-sequences.fasta"; [ -f "$ALN" ] || ALN=NONE

python "$R/phase3/build_asv_catalogue.py" \
  "$OUT/taxonomy_unite10/exp/taxonomy.tsv" \
  "$R/phase2/merged2/rep_exp/dna-sequences.fasta" \
  "$OUT/table_exp/feature-table.tsv" \
  "$ALN" \
  "$OUT/asv_catalogue" \
  "UNITE 10.0 (fungi, dynamic clustering), classify-consensus-vsearch"
echo "T1_T2 COMPLETE"
