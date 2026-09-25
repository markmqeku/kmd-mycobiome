#!/usr/bin/env bash
# KMD PROMPT 4, R6a. Reconcile the Jaccard axis-1 value (Figure 6 shows 34.1 percent, FINAL_NUMBERS and the
# text gave 34.4 percent). FINAL_NUMBERS cites item17b_ord_clean.py, which still carries DEPTH = 57141, the
# depth superseded by DEC-1. This runs that script unchanged at its recorded depth and, as a copy with only
# the depth line changed, at the locked depth 57,124 (same seed), then runs item23_position.py (already at
# 57,124) for the Christiana distance comparison. Outputs go to outputs_27-07-2026/clean/r6a/.
set -uo pipefail
cd <KMD_ROOT>/__reanalysis_2026-06/phase3
PY=<HOME>/miniconda3/envs/qiime2-amplicon-2024.10/bin/python
D=outputs_27-07-2026/clean/r6a; mkdir -p "$D"
$PY item17b_ord_clean.py > "$D/ord_depth57141_as_recorded.txt" 2>&1
sed 's/^DEPTH = 57141/DEPTH = 57124/' item17b_ord_clean.py > "$D/item17b_ord_clean_depth57124.py"
grep -n "^DEPTH" "$D/item17b_ord_clean_depth57124.py"
$PY "$D/item17b_ord_clean_depth57124.py" > "$D/ord_depth57124.txt" 2>&1
$PY item23_position.py > "$D/item23_position_57124.txt" 2>&1
grep -E "PCoA \(axis1|mean richness|difference" "$D/ord_depth57141_as_recorded.txt" "$D/ord_depth57124.txt"
