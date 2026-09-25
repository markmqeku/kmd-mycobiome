#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
WD=<KMD_ROOT>/__reanalysis_2026-06/phase2/unite
mkdir -p "$WD"
# NOTE: UNITE v8.0 (manuscript version) is NOT fetchable via RESCRIPt in q2-2024.10 (offers 8.2/8.3/9.0/10.0)
# and the PlutoF DOI page is a JS SPA (no static link). Using 8.2 = nearest available v8 release (Nilsson 2019). FLAGGED.
qiime rescript get-unite-data --p-version 8.2 --p-taxon-group fungi --p-cluster-id dynamic \
  --o-taxonomy "$WD/unite82_tax.qza" --o-sequences "$WD/unite82_seqs.qza" --verbose 2>"$WD/fetch.log"
echo "UNITE 8.2 FETCH DONE"; ls -la "$WD"/*.qza
