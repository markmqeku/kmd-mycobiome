#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2/otp_fixed
mkdir -p "$WD"
# 1) drop malformed records (base/qual length mismatch) from each read file
reformat.sh -Xmx3g tossbrokenreads qin=33 in="$RR/OT-P_S48_L001_R1_001.fastq" out="$WD/R1_clean.fastq.gz" overwrite=t 2>"$WD/reformat_r1.log"
reformat.sh -Xmx3g tossbrokenreads qin=33 in="$RR/OT-P_S48_L001_R2_001.fastq" out="$WD/R2_clean.fastq.gz" overwrite=t 2>"$WD/reformat_r2.log"
# 2) re-pair the cleaned files into matched pairs
repair.sh -Xmx3g in="$WD/R1_clean.fastq.gz" in2="$WD/R2_clean.fastq.gz" out="$WD/OT-P_R1.fastq.gz" out2="$WD/OT-P_R2.fastq.gz" outs="$WD/OT-P_singletons.fastq.gz" overwrite=t 2>"$WD/repair2.log"
echo -n "cleaned+matched OT-P pairs (R1): "; zcat "$WD/OT-P_R1.fastq.gz" 2>/dev/null | awk 'END{print NR/4}'
echo -n "                             (R2): "; zcat "$WD/OT-P_R2.fastq.gz" 2>/dev/null | awk 'END{print NR/4}'
