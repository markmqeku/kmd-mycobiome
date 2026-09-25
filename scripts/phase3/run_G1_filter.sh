#!/usr/bin/env bash
# G-1: build the fungi-only table by removing ONLY the 149 Searsia lancea host ASVs.
# The 78 Cryptococcus s.l. ASVs are fungi and are RETAINED. Original table untouched.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
B=<KMD_ROOT>/__reanalysis_2026-06
O=$B/phase3/outputs_27-07-2026
D=$O/nontarget_screen
F=$O/fungi_only
mkdir -p "$F"

# metadata file of host IDs to exclude
{ echo -e "Feature ID\thost"; while read -r id; do [ -n "$id" ] && echo -e "$id\thost"; done < "$D/host_asv_ids_ALL.txt"; } > "$F/host_ids.tsv"
echo "host ASVs to exclude: $(($(wc -l < "$F/host_ids.tsv")-1))"

qiime feature-table filter-features --i-table "$B/phase2/merged2/table.qza" \
  --m-metadata-file "$F/host_ids.tsv" --p-exclude-ids \
  --o-filtered-table "$F/table_fungi.qza"
qiime feature-table filter-seqs --i-data "$B/phase2/merged2/rep.qza" \
  --m-metadata-file "$F/host_ids.tsv" --p-exclude-ids \
  --o-filtered-data "$F/rep_fungi.qza"
qiime tools export --input-path "$F/table_fungi.qza" --output-path "$F/table_exp" >/dev/null
biom convert -i "$F/table_exp/feature-table.biom" -o "$F/table_exp/feature-table.tsv" --to-tsv
qiime tools export --input-path "$F/rep_fungi.qza" --output-path "$F/rep_exp" >/dev/null
echo "fungi-only ASVs: $(grep -c '^>' "$F/rep_exp/dna-sequences.fasta")"
echo "G1 FILTER DONE"
