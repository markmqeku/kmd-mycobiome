#!/usr/bin/env bash
# Phase 1 (corrected): per-read-length truncation for Arm B + working ITSxpress Arm A. Reuses sub/ subsamples. STOP after (report only).
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
WD=<KMD_ROOT>/__reanalysis_2026-06/phase1
ITS3=GCATCGATGAAGAACGCAGC; ITS4=TCCTCCGCTTATTGATATGC; TH=4
mkdir -p "$WD/B_bloem" "$WD/B_chrpta" "$WD/A2_itsxpress"

manifest(){ local out=$1; shift; printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$out"; for sid in "$@"; do printf "%s\t%s\t%s\n" "$sid" "$WD/sub/${sid}_R1.fastq.gz" "$WD/sub/${sid}_R2.fastq.gz" >> "$out"; done; }

armB(){ # $1=dir $2=truncF $3=truncR ; rest=samples
  local dir=$1 tf=$2 tr=$3; shift 3
  manifest "$WD/$dir/man.tsv" "$@"
  qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$WD/$dir/man.tsv" --output-path "$WD/$dir/demux.qza"
  qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/$dir/demux.qza" --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores $TH --o-trimmed-sequences "$WD/$dir/trimmed.qza" 2>"$WD/$dir/cutadapt.log"
  qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/$dir/trimmed.qza" --p-trunc-len-f $tf --p-trunc-len-r $tr --p-n-threads $TH --o-table "$WD/$dir/table.qza" --o-representative-sequences "$WD/$dir/rep.qza" --o-denoising-stats "$WD/$dir/stats.qza" 2>"$WD/$dir/dada2.log"
  qiime tools export --input-path "$WD/$dir/stats.qza" --output-path "$WD/$dir/stats_exp"
  echo "ARM B [$dir] trunc $tf/$tr DONE"
}

echo "### ARM B — Bloemfontein group (2x250) trunc 230/200 ###"
armB B_bloem 230 200 P16 P20 P28
echo "### ARM B — Christiana/Pretoria group (2x301) trunc 280/230 ###"
armB B_chrpta 280 230 25-C 28-C OL-C OL-P

echo "### ARM A — ITSxpress (all 7) then DADA2 paired trunc 0 ###"
MANA="$WD/A2_itsxpress/man.tsv"; printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANA"
for sid in P16 P20 P28 25-C 28-C OL-C OL-P; do
  o1="$WD/A2_itsxpress/${sid}_R1.fastq.gz"; o2="$WD/A2_itsxpress/${sid}_R2.fastq.gz"
  itsxpress --fastq "$WD/sub/${sid}_R1.fastq.gz" --fastq2 "$WD/sub/${sid}_R2.fastq.gz" --region ITS2 --taxa Fungi --cluster_id 1.0 --threads $TH --outfile "$o1" --outfile2 "$o2" --log "$WD/A2_itsxpress/${sid}.log" 2>>"$WD/A2_itsxpress/err.log" \
    && printf "%s\t%s\t%s\n" "$sid" "$o1" "$o2" >> "$MANA" || echo "ITSXPRESS FAIL $sid"
  # record how many pairs survived itsxpress
  n=$(zcat "$o1" 2>/dev/null | awk 'END{print NR/4}'); echo "  itsxpress $sid -> $n pairs"
done
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$MANA" --output-path "$WD/A2_itsxpress/demux.qza"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/A2_itsxpress/demux.qza" --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads $TH --o-table "$WD/A2_itsxpress/table.qza" --o-representative-sequences "$WD/A2_itsxpress/rep.qza" --o-denoising-stats "$WD/A2_itsxpress/stats.qza" 2>"$WD/A2_itsxpress/dada2.log"
qiime tools export --input-path "$WD/A2_itsxpress/stats.qza" --output-path "$WD/A2_itsxpress/stats_exp"
echo "ARM A DONE"
echo "P1b COMPLETE"
