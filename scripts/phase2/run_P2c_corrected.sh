#!/usr/bin/env bash
# PHASE 2 CORRECTED. C-0 rule: trunc-len set against POST-PRIMER length, asserted below.
# Group A (Bloem, ~231bp post-primer): DADA2 trunc-f 228 / trunc-r 200.  Group B: UNCHANGED (frozen).
# OT-P (salvaged 110,302 pairs, 2x301): cutadapt -> DADA2 280/230 as its own run, then merged.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2
ITS3=GCATCGATGAAGAACGCAGC; ITS4=TCCTCCGCTTATTGATATGC
mkdir -p "$WD/groupA2" "$WD/otp" "$WD/merged2"

echo "===== C-0 assertion: post-primer read lengths vs trunc ====="
python - <<PY
import zipfile,gzip
def minlen(qza,tag):
    z=zipfile.ZipFile(qza); n=[x for x in z.namelist() if tag in x and x.endswith(".fastq.gz")][0]
    d=gzip.decompress(z.read(n)).decode("latin1").split("\n")
    ls=[len(d[i]) for i in range(1,4*4000,4) if i<len(d) and d[i]]
    return min(ls), sorted(set(ls))[:3], max(ls)
for grp,qza in [("GroupA_Bloem","$WD/groupA/trimmed.qza")]:
    for tag,trunc in [("_R1_",228),("_R2_",200)]:
        mn,lo,mx=minlen(qza,tag); ok = trunc<=mn
        print(f"  {grp} {tag} post-primer len min={mn} max={mx} ; trunc={trunc} -> {'OK' if ok else 'FAIL(trunc>min!)'}")
PY

echo "===== C-1 Group A re-run DADA2 228/200 ====="
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/groupA/trimmed.qza" \
  --p-trunc-len-f 228 --p-trunc-len-r 200 --p-n-threads 6 \
  --o-table "$WD/groupA2/table.qza" --o-representative-sequences "$WD/groupA2/rep.qza" \
  --o-denoising-stats "$WD/groupA2/stats.qza" 2>"$WD/groupA2/dada2.log"
qiime tools export --input-path "$WD/groupA2/stats.qza" --output-path "$WD/groupA2/stats_exp"
echo "GROUP A2 DONE"

echo "===== OT-P import + cutadapt + DADA2 280/230 ====="
printf "sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n" > "$WD/otp/man.tsv"
printf "OT-P\t%s\t%s\n" "$WD/otp_fixed/OT-P_R1.fastq.gz" "$WD/otp_fixed/OT-P_R2.fastq.gz" >> "$WD/otp/man.tsv"
qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-format PairedEndFastqManifestPhred33V2 --input-path "$WD/otp/man.tsv" --output-path "$WD/otp/demux.qza"
qiime cutadapt trim-paired --i-demultiplexed-sequences "$WD/otp/demux.qza" --p-front-f $ITS3 --p-front-r $ITS4 --p-discard-untrimmed --p-cores 4 --o-trimmed-sequences "$WD/otp/trimmed.qza" 2>"$WD/otp/cutadapt.log"
qiime dada2 denoise-paired --i-demultiplexed-seqs "$WD/otp/trimmed.qza" --p-trunc-len-f 280 --p-trunc-len-r 230 --p-n-threads 4 --o-table "$WD/otp/table.qza" --o-representative-sequences "$WD/otp/rep.qza" --o-denoising-stats "$WD/otp/stats.qza" 2>"$WD/otp/dada2.log"
qiime tools export --input-path "$WD/otp/stats.qza" --output-path "$WD/otp/stats_exp"
echo "OT-P DONE"

echo "===== C-2 merge A2 + B + OT-P ====="
qiime feature-table merge --i-tables "$WD/groupA2/table.qza" "$WD/groupB/table.qza" "$WD/otp/table.qza" --o-merged-table "$WD/merged2/table.qza"
qiime feature-table merge-seqs --i-data "$WD/groupA2/rep.qza" "$WD/groupB/rep.qza" "$WD/otp/rep.qza" --o-merged-data "$WD/merged2/rep.qza"
qiime feature-table summarize --i-table "$WD/merged2/table.qza" --o-visualization "$WD/merged2/table_summary.qzv"
qiime tools export --input-path "$WD/merged2/table_summary.qzv" --output-path "$WD/merged2/table_summary_exp"
qiime tools export --input-path "$WD/merged2/rep.qza" --output-path "$WD/merged2/rep_exp"
echo -n "MERGED total ASVs: "; grep -c '^>' "$WD/merged2/rep_exp/dna-sequences.fasta"

echo "===== C-5 phylogeny on MERGED rep-seqs ====="
qiime phylogeny align-to-tree-mafft-fasttree --i-sequences "$WD/merged2/rep.qza" --p-n-threads 4 \
  --o-alignment "$WD/merged2/aligned.qza" --o-masked-alignment "$WD/merged2/masked.qza" \
  --o-tree "$WD/merged2/unrooted-tree.qza" --o-rooted-tree "$WD/merged2/rooted-tree.qza" 2>"$WD/merged2/phylogeny.log"
echo "P2c COMPLETE"
