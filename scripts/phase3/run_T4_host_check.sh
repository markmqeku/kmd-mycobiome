#!/usr/bin/env bash
# T-4 complementary test: build a small local reference of host-family (Anacardiaceae / Searsia) ITS
# plus a plant outgroup, and blast the non-target candidates against it. Independent of NCBI BLAST queueing.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
D=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen
mkdir -p "$D/hostref"
E="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

fetch () {  # $1=query label  $2=entrez query
  echo "  fetching: $1"
  ids=$(curl -s "$E/esearch.fcgi?db=nuccore&retmax=120&term=$(printf '%s' "$2" | sed 's/ /+/g')" | grep -oP '(?<=<Id>)[0-9]+' | tr '\n' ',' | sed 's/,$//')
  n=$(printf '%s' "$ids" | tr ',' '\n' | grep -c . || echo 0)
  echo "    ids: $n"
  [ "$n" -gt 0 ] && curl -s "$E/efetch.fcgi?db=nuccore&id=${ids}&rettype=fasta&retmode=text" >> "$D/hostref/host_ref.fasta"
  sleep 1
}
: > "$D/hostref/host_ref.fasta"
fetch "Searsia ITS"        'Searsia[Organism] AND (internal transcribed spacer[Title] OR ITS[Title])'
fetch "Rhus ITS"           'Rhus[Organism] AND (internal transcribed spacer[Title] OR ITS[Title])'
fetch "Anacardiaceae ITS"  'Anacardiaceae[Organism] AND internal transcribed spacer[Title]'
echo "  host reference seqs: $(grep -c '^>' "$D/hostref/host_ref.fasta" 2>/dev/null || echo 0)"

if [ "$(grep -c '^>' "$D/hostref/host_ref.fasta" 2>/dev/null || echo 0)" -gt 0 ]; then
  makeblastdb -in "$D/hostref/host_ref.fasta" -dbtype nucl -out "$D/hostref/hostdb" >/dev/null
  blastn -db "$D/hostref/hostdb" -query "$D/nontarget_query_all.fasta" \
    -outfmt "6 qseqid sacc pident length qcovs evalue bitscore stitle" \
    -max_target_seqs 3 -evalue 1e-5 -num_threads 4 -out "$D/blast_host.tsv"
  echo "  host-db hits: $(wc -l < "$D/blast_host.tsv")"
  python - <<'PY'
import csv
D="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/nontarget_screen"
best={}
for r in csv.reader(open(f"{D}/blast_host.tsv"),delimiter="\t"):
    if len(r)<8: continue
    q,pid,qcov,title=r[0],float(r[2]),float(r[4]),r[7]
    if q not in best or pid>best[q][0]: best[q]=(pid,qcov,title)
ab={r[0]:int(r[1]) for r in csv.reader(open(f"{D}/query_abundance.tsv"),delimiter="\t") if r and r[0]!="ASV_ID_md5"}
hom={r[0] for r in csv.reader(open(f"{D}/query_abundance.tsv"),delimiter="\t") if r and r[0]!="ASV_ID_md5" and r[2]=="1"}
tot=sum(ab.values())
hit=sum(ab.get(q,0) for q in best)
print(f"\n  ASVs with a host-family hit: {len(best)} / {len(ab)}")
print(f"  reads represented          : {hit:,} / {tot:,} of query set ({100*hit/tot:.1f}%)")
h=[q for q in best if q in hom]
print(f"  of the 78 Homophron ASVs   : {len(h)} hit the host reference")
strong=[(q,v) for q,v in best.items() if v[0]>=95 and v[1]>=80]
print(f"  strong hits (>=95% id, >=80% qcov): {len(strong)}")
for q,v in sorted(strong,key=lambda x:-ab.get(x[0],0))[:8]:
    print(f"    {q[:10]} reads={ab.get(q,0):>8,} pid={v[0]:.1f} qcov={v[1]:.0f}  {v[2][:70]}")
PY
fi
echo "T4_HOST DONE"
