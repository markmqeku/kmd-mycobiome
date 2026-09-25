#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
B=<KMD_ROOT>/__reanalysis_2026-06
O=$B/phase3/outputs_27-07-2026
# rebuild ASV catalogue on the FUNGI-ONLY table with the final rank-threshold taxonomy
mkdir -p "$O/asv_catalogue_final"
python "$B/phase3/build_asv_catalogue.py" \
  "$O/taxonomy_unite10_id97/exp/taxonomy.tsv" \
  "$O/fungi_only/rep_exp/dna-sequences.fasta" \
  "$O/fungi_only/table_exp/feature-table.tsv" \
  NONE \
  "$O/asv_catalogue_final" \
  "UNITE 10.0 (fungi, dynamic); classify-consensus-vsearch; rank-specific identity thresholds (species>=97, genus>=95, family>=90, order>=85)"
