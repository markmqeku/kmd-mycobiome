#!/usr/bin/env bash
# Item 10a: alpha-rarefaction on the fungi-only table, 16 CORE samples only.
# Item 10b: clean machine-readable iNEXT per-sample Hill export.
# Narrowly scoped recomputation. Nothing in FINAL_NUMBERS changes.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export MPLBACKEND=Agg
B=<KMD_ROOT>/__reanalysis_2026-06
O=$B/phase3/outputs_27-07-2026
G=$O/gap_closure; mkdir -p "$G"

# core-16 metadata subset
python - <<'PY'
import csv
O="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
rows=list(csv.reader(open(f"{O}/metadata_v2_27-07-2026.tsv"),delimiter="\t"))
hdr,types=rows[0],rows[1]; di={c:i for i,c in enumerate(hdr)}
keep=[r for r in rows[2:] if r and r[0].strip() and r[di["location"]] in ("Christiana","Pretoria")]
with open(f"{O}/gap_closure/metadata_core16.tsv","w",newline="") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(hdr); w.writerow(types); w.writerows(keep)
print("core-16 metadata rows:",len(keep))
PY

qiime feature-table filter-samples --i-table "$O/fungi_only/table_fungi.qza" \
  --m-metadata-file "$G/metadata_core16.tsv" --o-filtered-table "$G/table_core16.qza"
qiime feature-table summarize --i-table "$G/table_core16.qza" --o-visualization "$G/table_core16.qzv"
echo "10a: core-16 table built"

qiime diversity alpha-rarefaction --i-table "$G/table_core16.qza" \
  --p-max-depth 57141 --p-steps 20 --p-iterations 10 \
  --m-metadata-file "$G/metadata_core16.tsv" \
  --o-visualization "$G/alpha_rarefaction_core16.qzv" 2>"$G/rarefaction.log"
qiime tools export --input-path "$G/alpha_rarefaction_core16.qzv" --output-path "$G/rarefaction_exp"
echo "10a: rarefaction done"; ls "$G/rarefaction_exp" | head -6
