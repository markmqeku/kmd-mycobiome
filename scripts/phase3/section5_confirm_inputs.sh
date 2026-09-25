#!/usr/bin/env bash
# Section 5 — confirm inputs only. No computation on the data beyond summarising. No edits.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
R=<KMD_ROOT>/__reanalysis_2026-06
OUT=$R/phase3/section5
mkdir -p "$OUT"

echo "### QIIME2 version"; qiime --version | head -2

echo; echo "### FEATURE TABLE summary"
qiime feature-table summarize --i-table "$R/phase2/merged2/table.qza" --o-visualization "$OUT/table_sum.qzv" >/dev/null 2>&1
qiime tools export --input-path "$OUT/table_sum.qzv" --output-path "$OUT/table_sum_exp" >/dev/null 2>&1
python - <<'PY'
import re,os,json
p="<KMD_ROOT>/__reanalysis_2026-06/phase3/section5/table_sum_exp"
# per-sample depths live in sample-frequency-detail.html; parse the embedded values
import glob
html=open(os.path.join(p,"sample-frequency-detail.html"),encoding="utf-8",errors="replace").read()
nums=re.findall(r'"([\w\-]+)"\s*:\s*([0-9]+\.?[0-9]*)',html)
vals={k:float(v) for k,v in nums if not k.isdigit()}
if vals:
    import statistics as st
    d=sorted(vals.values())
    print(f"samples in table: {len(d)}")
    print(f"min depth  : {int(d[0])}")
    print(f"median depth: {int(st.median(d))}")
    print(f"max depth  : {int(d[-1])}")
    print(f"total reads: {int(sum(d))}")
    print("samples below 5000:", sorted([k for k,v in vals.items() if v<5000]))
PY

echo; echo "### REP-SEQS feature count (authoritative ASV count)"
qiime tools export --input-path "$R/phase2/merged2/rep.qza" --output-path "$OUT/rep_exp" >/dev/null 2>&1
echo -n "ASVs in rep-seqs: "; grep -c '^>' "$OUT/rep_exp/dna-sequences.fasta"

echo; echo "### TAXONOMY"
echo -n "taxonomy rows (minus header): "; tail -n +2 "$R/phase2/taxonomy/exp/taxonomy.tsv" | grep -vc '^#' || true

echo; echo "### R / iNEXT / hillR availability"
if command -v R >/dev/null 2>&1; then
  R --version | head -1
  Rscript -e 'for (p in c("iNEXT","hillR","vegan")) cat(p, ifelse(requireNamespace(p,quietly=TRUE),"INSTALLED","MISSING"), "\n")' 2>/dev/null || echo "Rscript present but package check failed"
else
  echo "R: NOT INSTALLED in this env"
fi
echo; echo "### python fallback libs"
python -c "import skbio,pandas,scipy; print('scikit-bio',skbio.__version__,'| pandas',pandas.__version__,'| scipy ok')" 2>/dev/null || echo "skbio check failed"
