#!/usr/bin/env bash
# D-1: import intact provisional set (44 samples) + demux summarize + export quality plots. STOP after.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10

RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase0
mkdir -p "$WD/D1_demux"
MAN="$WD/D1_demux/manifest_intact44.tsv"

# excluded: truncated/empty (P6,P15,P24), broken pair (YT-P), R1!=R2 (OT-P); extras CR,PR,S1-S4 excluded by not being P#/-C/-P design ids
declare -A EXCL=( [P6]=1 [P15]=1 [P24]=1 [YT-P]=1 [OT-P]=1 [CR]=1 [PR]=1 [S1]=1 [S2]=1 [S3]=1 [S4]=1 )

printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MAN"
for r1 in "$RR"/*_R1_001.fastq "$RR"/*_R1_001.fastq.gz; do
  [ -e "$r1" ] || continue
  b=$(basename "$r1")
  case "$b" in *.filepart) continue;; esac
  sid=$(echo "$b" | sed -E 's/_S[0-9]+_L001_R1_001\.fastq(\.gz)?$//')
  [ "${EXCL[$sid]:-0}" = "1" ] && continue
  r2="${r1/_R1_/_R2_}"
  if [ ! -e "$r2" ]; then echo "WARN no R2 for $sid" >&2; continue; fi
  printf "%s\t%s\t%s\n" "$sid" "$r1" "$r2" >> "$MAN"
done
echo "Manifest samples: $(($(wc -l < "$MAN")-1))"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-format PairedEndFastqManifestPhred33V2 \
  --input-path "$MAN" \
  --output-path "$WD/D1_demux/demux.qza"
echo "IMPORT OK"

qiime demux summarize \
  --i-data "$WD/D1_demux/demux.qza" \
  --o-visualization "$WD/D1_demux/demux_summary.qzv"
echo "SUMMARIZE OK"

# export the qzv contents (contains forward/reverse quality plots + per-sample counts)
qiime tools export --input-path "$WD/D1_demux/demux_summary.qzv" --output-path "$WD/D1_demux/demux_summary_export"
echo "EXPORT OK -> $WD/D1_demux/demux_summary_export"
ls "$WD/D1_demux/demux_summary_export" | head -40
echo "D-1 COMPLETE"
