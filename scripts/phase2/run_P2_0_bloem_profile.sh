#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2
mkdir -p "$WD/P2_0_bloem"
MAN="$WD/P2_0_bloem/manifest_bloem29.tsv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MAN"
for n in 1 2 3 4 5 7 8 9 10 11 12 13 14 16 17 18 19 20 21 22 23 25 26 27 28 29 30 31 32; do
  r1=$(ls "$RR"/P${n}_S*_L001_R1_001.fastq.gz 2>/dev/null | head -1)
  r2=$(ls "$RR"/P${n}_S*_L001_R2_001.fastq.gz 2>/dev/null | head -1)
  [ -n "$r1" ] && printf "P%s\t%s\t%s\n" "$n" "$r1" "$r2" >> "$MAN"
done
echo "Bloem samples: $(($(wc -l < "$MAN")-1))"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$MAN" --output-path "$WD/P2_0_bloem/demux.qza"
qiime demux summarize --i-data "$WD/P2_0_bloem/demux.qza" --o-visualization "$WD/P2_0_bloem/demux_sum.qzv"
qiime tools export --input-path "$WD/P2_0_bloem/demux_sum.qzv" --output-path "$WD/P2_0_bloem/exp"
echo "P2-0 DONE"
