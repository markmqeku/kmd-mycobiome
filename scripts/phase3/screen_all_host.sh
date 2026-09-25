#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
B=<KMD_ROOT>/__reanalysis_2026-06
D=$B/phase3/outputs_27-07-2026/nontarget_screen
blastn -db "$D/hostref/hostdb" -query "$B/phase2/merged2/rep_exp/dna-sequences.fasta" \
  -outfmt "6 qseqid sacc pident qcovs" -max_target_seqs 3 -evalue 1e-5 -num_threads 4 \
  -out "$D/blast_host_ALL.tsv" 2>/dev/null
echo -n "ASVs (of 1322) with any host hit: "; cut -f1 "$D/blast_host_ALL.tsv" | sort -u | wc -l
awk -F'\t' '$3>=95 && $4>=80 {print $1}' "$D/blast_host_ALL.tsv" | sort -u > "$D/host_asv_ids_ALL.txt"
echo -n "STRONG host ASVs (>=95% id, >=80% qcov): "; wc -l < "$D/host_asv_ids_ALL.txt"
python - <<'PY'
B="<KMD_ROOT>/__reanalysis_2026-06"
D=f"{B}/phase3/outputs_27-07-2026/nontarget_screen"
strong={l.strip() for l in open(f"{D}/host_asv_ids_ALL.txt") if l.strip()}
prev={l.strip() for l in open(f"{D}/host_asv_ids.txt") if l.strip()}
lines=[l.rstrip("\n") for l in open(f"{B}/phase3/outputs_27-07-2026/table_exp/feature-table.tsv") if l.strip() and not l.startswith("# Constructed")]
hdr=lines[0].lstrip("#").split("\t"); ns=len([h for h in hdr[1:] if h.strip()])
tot=0; hst=0
for l in lines[1:]:
    p=l.split("\t"); v=sum(int(float(x)) for x in p[1:1+ns])
    tot+=v
    if p[0] in strong: hst+=v
print(f"host reads (full screen): {hst:,} / {tot:,} = {100*hst/tot:.1f}% of dataset")
print(f"newly found beyond the unassigned-set screen: {len(strong-prev)} ASVs")
PY
