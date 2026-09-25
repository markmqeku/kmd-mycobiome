#!/usr/bin/env bash
# G-4a: train naive-Bayes classifier on UNITE 10.0 and classify the fungi-only rep-seqs at confidence 0.7.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
export JOBLIB_TEMP_FOLDER=/tmp
B=<KMD_ROOT>/__reanalysis_2026-06
O=$B/phase3/outputs_27-07-2026
mkdir -p "$O/taxonomy_sklearn"

echo "### training naive-Bayes on UNITE 10.0 (this is the memory-heavy step)"
qiime feature-classifier fit-classifier-naive-bayes \
  --i-reference-reads "$O/unite_current/unite10_seqs.qza" \
  --i-reference-taxonomy "$O/unite_current/unite10_tax.qza" \
  --o-classifier "$O/taxonomy_sklearn/unite10_nb.qza" 2>"$O/taxonomy_sklearn/train.log"
rc=$?
if [ $rc -ne 0 ] || [ ! -f "$O/taxonomy_sklearn/unite10_nb.qza" ]; then
  echo "TRAIN FAILED (rc=$rc) - see train.log"; tail -15 "$O/taxonomy_sklearn/train.log"; exit 1
fi
echo "TRAIN OK"

echo "### classify-sklearn, confidence 0.7"
qiime feature-classifier classify-sklearn \
  --i-classifier "$O/taxonomy_sklearn/unite10_nb.qza" \
  --i-reads "$O/fungi_only/rep_fungi.qza" \
  --p-confidence 0.7 --p-reads-per-batch 2000 --p-n-jobs 1 \
  --o-classification "$O/taxonomy_sklearn/taxonomy.qza" 2>"$O/taxonomy_sklearn/classify.log"
rc=$?
if [ $rc -ne 0 ]; then echo "CLASSIFY FAILED (rc=$rc)"; tail -15 "$O/taxonomy_sklearn/classify.log"; exit 1; fi
qiime tools export --input-path "$O/taxonomy_sklearn/taxonomy.qza" --output-path "$O/taxonomy_sklearn/exp"
echo "G4A DONE"
