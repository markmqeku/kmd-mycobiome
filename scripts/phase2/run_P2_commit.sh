#!/usr/bin/env bash
# PHASE 2 committing denoise — per-run DADA2, merge tables, phylogeny. Parameters below are the ARCHIVED Methods values.
# GROUP A (Bloemfontein, 2x250, n=29): cutadapt ITS3/ITS4 -> DADA2 trunc-len-f 245 trunc-len-r 200
# GROUP B (Christiana+Pretoria, 2x301, n=15): cutadapt ITS3/ITS4 -> DADA2 trunc-len-f 280 trunc-len-r 230
# (OT-P excluded: corrupt FASTQ record; YT-P/P6/P15/P24 written off.)  NO ITSxpress. No rarefaction/core-metrics.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2
ITS3=GCATCGATGAAGAACGCAGC; ITS4=TCCTCCGCTTATTGATATGC
mkdir -p "$WD/groupA" "$WD/groupB" "$WD/merged"

echo "===== GROUP A (Bloem) cutadapt + DADA2 245/200 ====="
qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/P2_0_bloem/demux.qza" \
  --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores 6 \
  --o-trimmed-sequences "$WD/groupA/trimmed.qza" 2>"$WD/groupA/cutadapt.log"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/groupA/trimmed.qza" \
  --p-trunc-len-f 245 --p-trunc-len-r 200 --p-n-threads 6 \
  --o-table "$WD/groupA/table.qza" --o-representative-sequences "$WD/groupA/rep.qza" \
  --o-denoising-stats "$WD/groupA/stats.qza" 2>"$WD/groupA/dada2.log"
qiime tools export --input-path "$WD/groupA/stats.qza" --output-path "$WD/groupA/stats_exp"
echo "GROUP A DONE"

echo "===== GROUP B (Christiana+Pretoria) build manifest + import ====="
MANB="$WD/groupB/manifest.tsv"; printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANB"
for sid in 25-C 26-C 27-C OL-C OM-C OT-C YL-C YM-C YT-C 28-C 29-C OL-P OM-P YL-P YM-P; do
  r1=$(ls "$RR"/${sid}_S*_L001_R1_001.fastq 2>/dev/null | head -1)
  r2=$(ls "$RR"/${sid}_S*_L001_R2_001.fastq 2>/dev/null | head -1)
  [ -n "$r1" ] && printf "%s\t%s\t%s\n" "$sid" "$r1" "$r2" >> "$MANB" || echo "MISSING $sid" >&2
done
echo "Group B samples: $(($(wc -l < "$MANB")-1))"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$MANB" --output-path "$WD/groupB/demux.qza"
qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/groupB/demux.qza" \
  --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores 6 \
  --o-trimmed-sequences "$WD/groupB/trimmed.qza" 2>"$WD/groupB/cutadapt.log"
echo "===== GROUP B DADA2 280/230 ====="
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/groupB/trimmed.qza" \
  --p-trunc-len-f 280 --p-trunc-len-r 230 --p-n-threads 3 \
  --o-table "$WD/groupB/table.qza" --o-representative-sequences "$WD/groupB/rep.qza" \
  --o-denoising-stats "$WD/groupB/stats.qza" 2>"$WD/groupB/dada2.log"
qiime tools export --input-path "$WD/groupB/stats.qza" --output-path "$WD/groupB/stats_exp"
echo "GROUP B DONE"

echo "===== MERGE tables + seqs ====="
qiime feature-table merge --i-tables "$WD/groupA/table.qza" "$WD/groupB/table.qza" --o-merged-table "$WD/merged/table.qza"
qiime feature-table merge-seqs --i-data "$WD/groupA/rep.qza" "$WD/groupB/rep.qza" --o-merged-data "$WD/merged/rep.qza"
qiime feature-table summarize --i-table "$WD/merged/table.qza" --o-visualization "$WD/merged/table_summary.qzv"
qiime tools export --input-path "$WD/merged/table_summary.qzv" --output-path "$WD/merged/table_summary_exp"
qiime tools export --input-path "$WD/merged/rep.qza" --output-path "$WD/merged/rep_exp"
echo -n "MERGED total ASVs: "; grep -c '^>' "$WD/merged/rep_exp/dna-sequences.fasta"

echo "===== PHYLOGENY (MAFFT -> mask -> FastTree -> midpoint root) ====="
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$WD/merged/rep.qza" \
  --p-n-threads 4 \
  --o-alignment "$WD/merged/aligned.qza" --o-masked-alignment "$WD/merged/masked.qza" \
  --o-tree "$WD/merged/unrooted-tree.qza" --o-rooted-tree "$WD/merged/rooted-tree.qza" 2>"$WD/merged/phylogeny.log"
echo "PHYLOGENY DONE"
echo "P2 COMMIT COMPLETE"
