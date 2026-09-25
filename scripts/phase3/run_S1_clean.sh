#!/usr/bin/env bash
# 17b: rebuild the core-16 rarefaction on the CORRECTED table, so Figure S1 matches the
# corrected data and the re-locked depth. Same QIIME2 route as the original run.
set -euo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
O=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026
C=$O/clean; G=$C/gap_closure; mkdir -p "$G"
cd "$G"
python - <<'PY'
import csv
O="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
C=f"{O}/clean"
meta={}
with open(f"{O}/metadata_v2_27-07-2026.tsv") as fh:
    rd=csv.reader(fh,delimiter="\t"); h=next(rd); next(rd); di={c:i for i,c in enumerate(h)}
    for r in rd:
        if r and r[0].strip(): meta[r[0]]={c:r[di[c]] for c in h}
lines=[l.rstrip("\n") for l in open(f"{C}/feature-table-clean.tsv") if l.strip()]
hdr=lines[0].lstrip("#").split("\t"); cols=[x for x in hdr[1:] if x.strip()]
core=[s for s in cols if meta.get(s) and meta[s]["location"] in ("Christiana","Pretoria")]
assert len(core)==16, core
idx=[cols.index(s) for s in core]
with open(f"{C}/gap_closure/table_core16_clean.tsv","w",newline="") as fh:
    fh.write("# Constructed from biom file\n")
    fh.write("#OTU ID\t"+"\t".join(core)+"\n")
    for l in lines[1:]:
        q=l.split("\t"); vals=[q[1+i] for i in idx]
        if any(float(v)>0 for v in vals): fh.write(q[0]+"\t"+"\t".join(vals)+"\n")
# metadata restricted to the core 16
with open(f"{C}/gap_closure/metadata_core16.tsv","w",newline="") as fh:
    fh.write("sample-id\tlocation\tcondition_std\ttissue\n")
    fh.write("#q2:types\tcategorical\tcategorical\tcategorical\n")
    for s in core:
        m=meta[s]; fh.write(f"{s}\t{m['location']}\t{m['condition_std']}\t{m['tissue']}\n")
print("wrote core-16 clean table and metadata")
PY
biom convert -i table_core16_clean.tsv -o table_core16_clean.biom --table-type="OTU table" --to-hdf5
qiime tools import --input-path table_core16_clean.biom --type 'FeatureTable[Frequency]' \
  --input-format BIOMV210Format --output-path table_core16_clean.qza
qiime diversity alpha-rarefaction --i-table table_core16_clean.qza \
  --p-max-depth 57124 --p-min-depth 1000 --p-steps 20 --p-iterations 10 \
  --p-metrics observed_features --m-metadata-file metadata_core16.tsv \
  --o-visualization alpha_rarefaction_core16_clean.qzv
qiime tools export --input-path alpha_rarefaction_core16_clean.qzv --output-path rarefaction_exp
ls rarefaction_exp/*.csv
echo "S1 INPUTS READY"
