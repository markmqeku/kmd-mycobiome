#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
RR=<KMD_ROOT>/__restructured_16-05-2026/raw_reads
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2/otp_fixed
mkdir -p "$WD"
repair.sh -Xmx3g ignorebadquality qin=33 in="$RR/OT-P_S48_L001_R1_001.fastq" in2="$RR/OT-P_S48_L001_R2_001.fastq" \
  out="$WD/OT-P_R1.fastq.gz" out2="$WD/OT-P_R2.fastq.gz" outs="$WD/OT-P_singletons.fastq.gz" overwrite=t 2>"$WD/repair.log"
echo -n "matched R1 pairs: "; zcat "$WD/OT-P_R1.fastq.gz" 2>/dev/null | awk 'END{print NR/4}'
echo -n "matched R2 pairs: "; zcat "$WD/OT-P_R2.fastq.gz" 2>/dev/null | awk 'END{print NR/4}'
grep -iE "Pairs:|Singletons:|Reads:|Result" "$WD/repair.log" | head
