#!/usr/bin/env bash
# Phase 1 P1-1/P1-2: primer removal % + ITSxpress-vs-fixed-truncation DADA2 comparison on a representative subset.
# STOP after: report primer %, both DADA2 stats tables (merge retention). No full-45 denoise. Params provisional (for comparison), Mark-approved before committing.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg

RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase1
mkdir -p "$WD/sub" "$WD/armA_itsxpress" "$WD/armB_fixed"
ITS3=GCATCGATGAAGAACGCAGC
ITS4=TCCTCCGCTTATTGATATGC
NSUB=25000   # read pairs per sample (bounds runtime; % retention is scale-invariant)
THREADS=4

# representative subset: sample-id : source-basename (compression auto)
declare -A SRC=(
 [P16]=P16_S74 [P20]=P20_S78 [P28]=P28_S86
 [25-C]=25-C_S53 [OL-C]=OL-C_S41
 [28-C]=28-C_S56 [OL-P]=OL-P_S47 )

echo "### subsampling $NSUB pairs/sample ###"
for sid in "${!SRC[@]}"; do
  base=${SRC[$sid]}
  for RN in R1 R2; do
    src_gz="$RR/${base}_L001_${RN}_001.fastq.gz"; src_fq="$RR/${base}_L001_${RN}_001.fastq"
    out="$WD/sub/${sid}_${RN}.fastq.gz"
    if [ -e "$src_gz" ]; then zcat "$src_gz" | head -n $((NSUB*4)) | gzip > "$out"
    else head -n $((NSUB*4)) "$src_fq" | gzip > "$out"; fi
  done
done
echo "subsampled: $(ls "$WD/sub" | wc -l) files"

# ---------- ARM B: cutadapt primer removal (P1-1) + fixed-truncation DADA2 (paired) ----------
echo "### ARM B: manifest + import ###"
MANB="$WD/armB_fixed/manifest.tsv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANB"
for sid in "${!SRC[@]}"; do printf "%s\t%s\t%s\n" "$sid" "$WD/sub/${sid}_R1.fastq.gz" "$WD/sub/${sid}_R2.fastq.gz" >> "$MANB"; done
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$MANB" --output-path "$WD/armB_fixed/demux.qza"

echo "### P1-1: cutadapt remove ITS3/ITS4 (discard-untrimmed) ###"
qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/armB_fixed/demux.qza" \
  --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores $THREADS \
  --o-trimmed-sequences "$WD/armB_fixed/trimmed.qza" --verbose 2>"$WD/armB_fixed/cutadapt.log"
qiime demux summarize --i-data "$WD/armB_fixed/demux.qza"   --o-visualization "$WD/armB_fixed/demux_sum.qzv"
qiime demux summarize --i-data "$WD/armB_fixed/trimmed.qza" --o-visualization "$WD/armB_fixed/trimmed_sum.qzv"
qiime tools export --input-path "$WD/armB_fixed/demux_sum.qzv"   --output-path "$WD/armB_fixed/demux_exp"
qiime tools export --input-path "$WD/armB_fixed/trimmed_sum.qzv" --output-path "$WD/armB_fixed/trimmed_exp"

echo "### ARM B: DADA2 denoise-paired (provisional trunc F=280 R=230) ###"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/armB_fixed/trimmed.qza" \
  --p-trunc-len-f 280 --p-trunc-len-r 230 --p-n-threads $THREADS \
  --o-table "$WD/armB_fixed/table.qza" --o-representative-sequences "$WD/armB_fixed/repseq.qza" \
  --o-denoising-stats "$WD/armB_fixed/stats.qza" 2>"$WD/armB_fixed/dada2.log"
qiime tools export --input-path "$WD/armB_fixed/stats.qza" --output-path "$WD/armB_fixed/stats_exp"
echo "ARM B DONE"

# ---------- ARM A: ITSxpress (paired, trim 5.8S/LSU flanks) + DADA2 (paired, trunc 0) ----------
echo "### ARM A: itsxpress per sample ###"
MANA="$WD/armA_itsxpress/manifest.tsv"
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$MANA"
for sid in "${!SRC[@]}"; do
  o1="$WD/armA_itsxpress/${sid}_R1.fastq.gz"; o2="$WD/armA_itsxpress/${sid}_R2.fastq.gz"
  itsxpress --fastq "$WD/sub/${sid}_R1.fastq.gz" --fastq2 "$WD/sub/${sid}_R2.fastq.gz" \
    --region ITS2 --taxa Fungi --cluster_id 1.0 --threads $THREADS \
    --outfile "$o1" --outfile2 "$o2" --log "$WD/armA_itsxpress/${sid}.log" 2>>"$WD/armA_itsxpress/itsxpress_err.log" \
    && printf "%s\t%s\t%s\n" "$sid" "$o1" "$o2" >> "$MANA" || echo "ITSXPRESS FAILED for $sid" >&2
done
echo "itsxpress produced: $(($(wc -l < "$MANA")-1)) samples"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$MANA" --output-path "$WD/armA_itsxpress/demux.qza"
echo "### ARM A: DADA2 denoise-paired (trunc 0/0; flanks already trimmed) ###"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/armA_itsxpress/demux.qza" \
  --p-trunc-len-f 0 --p-trunc-len-r 0 --p-n-threads $THREADS \
  --o-table "$WD/armA_itsxpress/table.qza" --o-representative-sequences "$WD/armA_itsxpress/repseq.qza" \
  --o-denoising-stats "$WD/armA_itsxpress/stats.qza" 2>"$WD/armA_itsxpress/dada2.log"
qiime tools export --input-path "$WD/armA_itsxpress/stats.qza" --output-path "$WD/armA_itsxpress/stats_exp"
echo "ARM A DONE"
echo "P1 COMPLETE"
