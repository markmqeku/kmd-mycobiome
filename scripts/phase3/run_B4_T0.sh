#!/usr/bin/env bash
# B4: install iNEXT + hillR from CRAN (STOP on failure, no silent Python fallback).
# T-0: fetch most current UNITE release via RESCRIPt, report exact version, re-classify.
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
R=<KMD_ROOT>/__reanalysis_2026-06
OUT=$R/phase3/outputs_27-07-2026
mkdir -p "$OUT" "$OUT/unite_current"

echo "########## B4: CRAN install ##########"
Rscript -e 'options(repos=c(CRAN="https://cloud.r-project.org"));
 need <- c("iNEXT","hillR");
 for (p in need) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, quiet=TRUE);
 for (p in need) cat(p, ifelse(requireNamespace(p, quietly=TRUE), "INSTALLED", "FAILED"), "\n")' 2>&1 | tail -20

echo; echo "########## T-0: RESCRIPt UNITE versions available ##########"
qiime rescript get-unite-data --help 2>&1 | grep -A6 -- "--p-version" | head -12

echo; echo "########## T-0: fetch UNITE 10.0 (most current offered) ##########"
qiime rescript get-unite-data --p-version 10.0 --p-taxon-group fungi --p-cluster-id dynamic \
  --o-taxonomy "$OUT/unite_current/unite10_tax.qza" --o-sequences "$OUT/unite_current/unite10_seqs.qza" \
  --verbose 2>&1 | tail -15
ls -la "$OUT/unite_current"/*.qza 2>/dev/null
echo "B4_T0 STEP DONE"
