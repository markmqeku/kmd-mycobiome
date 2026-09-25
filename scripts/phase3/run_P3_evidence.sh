#!/usr/bin/env bash
# PHASE 3 EVIDENCE ONLY: alpha-rarefaction curves + Good's coverage + exported table for retention analysis.
# Does NOT run core-metrics. Does NOT choose a depth. Fork reserved for Mark + supervisor.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
WD=<KMD_ROOT>/__reanalysis_2026-06
P2=$WD/phase2; P3=$WD/phase3
mkdir -p "$P3"

# 1) export merged table to TSV (per-sample frequencies + ASV presence for retention/coverage)
qiime tools export --input-path "$P2/merged2/table.qza" --output-path "$P3/table_exp"
biom convert -i "$P3/table_exp/feature-table.biom" -o "$P3/table_exp/feature-table.tsv" --to-tsv
echo "TABLE EXPORTED"

# 2) alpha-rarefaction curves (x-axis range only; NOT a depth choice)
qiime diversity alpha-rarefaction \
  --i-table "$P2/merged2/table.qza" --i-phylogeny "$P2/merged2/rooted-tree.qza" \
  --p-max-depth 50000 --p-steps 20 --p-iterations 10 \
  --m-metadata-file "$WD/phase0/metadata_consolidated.tsv" \
  --o-visualization "$P3/alpha_rarefaction.qzv" 2>"$P3/rarefaction.log"
qiime tools export --input-path "$P3/alpha_rarefaction.qzv" --output-path "$P3/rarefaction_exp"
echo "RAREFACTION CURVES DONE"

# 3) Good's coverage per sample (completeness evidence for coverage-based options)
qiime diversity alpha --i-table "$P2/merged2/table.qza" --p-metric goods_coverage --o-alpha-diversity "$P3/goods_coverage.qza" 2>"$P3/goods.log"
qiime tools export --input-path "$P3/goods_coverage.qza" --output-path "$P3/goods_exp"
# observed features (unrarefied) for context
qiime diversity alpha --i-table "$P2/merged2/table.qza" --p-metric observed_features --o-alpha-diversity "$P3/observed.qza" 2>>"$P3/goods.log"
qiime tools export --input-path "$P3/observed.qza" --output-path "$P3/observed_exp"
echo "COVERAGE DONE"
echo "P3 EVIDENCE COMPLETE"
