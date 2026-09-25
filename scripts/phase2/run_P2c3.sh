#!/usr/bin/env bash
# PHASE 2 CORRECTED v3 — single-thread DADA2 (RAM-safe), verbose. Group A 228/200 + OT-P 280/230 + merge + phylogeny.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2
ITS3=GCATCGATGAAGAACGCAGC; ITS4=TCCTCCGCTTATTGATATGC
rm -rf "$WD/groupA2" "$WD/otp" "$WD/merged2"; mkdir -p "$WD/groupA2" "$WD/otp" "$WD/merged2"

echo "===== C-1 Group A DADA2 228/200 (1 thread, verbose) ====="
if qiime dada2 denoise-paired --verbose --i-demultiplexed-seqs "$WD/groupA/trimmed.qza" \
  --p-trunc-len-f 228 --p-trunc-len-r 200 --p-n-threads 3 \
  --o-table "$WD/groupA2/table.qza" --o-representative-sequences "$WD/groupA2/rep.qza" \
  --o-denoising-stats "$WD/groupA2/stats.qza" >"$WD/groupA2/dada2.log" 2>&1; then
  qiime tools export --input-path "$WD/groupA2/stats.qza" --output-path "$WD/groupA2/stats_exp"; echo "GROUP A2 DONE"
else echo "GROUP A2 FAILED — see groupA2/dada2.log"; tail -20 "$WD/groupA2/dada2.log"; exit 1; fi

echo "===== OT-P import + cutadapt + DADA2 280/230 (1 thread) ====="
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$WD/otp/man.tsv"
printf "OT-P\t%s\t%s\n" "$WD/otp_fixed/OT-P_R1.fastq.gz" "$WD/otp_fixed/OT-P_R2.fastq.gz" >> "$WD/otp/man.tsv"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$WD/otp/man.tsv" --output-path "$WD/otp/demux.qza"
qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/otp/demux.qza" --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores 2 --o-trimmed-sequences "$WD/otp/trimmed.qza" 2>"$WD/otp/cutadapt.log"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/otp/trimmed.qza" --p-trunc-len-f 280 --p-trunc-len-r 230 --p-n-threads 1 --o-table "$WD/otp/table.qza" --o-representative-sequences "$WD/otp/rep.qza" --o-denoising-stats "$WD/otp/stats.qza" >"$WD/otp/dada2.log" 2>&1
qiime tools export --input-path "$WD/otp/stats.qza" --output-path "$WD/otp/stats_exp"; echo "OT-P DONE"

echo "===== C-2 merge A2 + B + OT-P ====="
qiime feature-table merge --i-tables "$WD/groupA2/table.qza" "$WD/groupB/table.qza" "$WD/otp/table.qza" --o-merged-table "$WD/merged2/table.qza"
qiime feature-table merge-seqs --i-data "$WD/groupA2/rep.qza" "$WD/groupB/rep.qza" "$WD/otp/rep.qza" --o-merged-data "$WD/merged2/rep.qza"
qiime feature-table summarize --i-table "$WD/merged2/table.qza" --o-visualization "$WD/merged2/table_summary.qzv"
qiime tools export --input-path "$WD/merged2/table_summary.qzv" --output-path "$WD/merged2/table_summary_exp"
qiime tools export --input-path "$WD/merged2/rep.qza" --output-path "$WD/merged2/rep_exp"
echo -n "MERGED total ASVs: "; grep -c '^>' "$WD/merged2/rep_exp/dna-sequences.fasta"

echo "===== C-5 phylogeny on MERGED rep-seqs (2 threads) ====="
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$WD/merged2/rep.qza" --p-n-threads 2 \
  --o-alignment "$WD/merged2/aligned.qza" --o-masked-alignment "$WD/merged2/masked.qza" \
  --o-tree "$WD/merged2/unrooted-tree.qza" --o-rooted-tree "$WD/merged2/rooted-tree.qza" 2>"$WD/merged2/phylogeny.log"
echo "P2c3 COMPLETE"
