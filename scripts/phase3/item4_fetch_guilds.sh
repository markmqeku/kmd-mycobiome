#!/usr/bin/env bash
set -uo pipefail
source ~/miniconda3/etc/profile.d/conda.sh; conda activate qiime2-amplicon-2024.10
G=<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/guilds
mkdir -p "$G"
echo "### FungalTraits (Polme et al. 2020, Fungal Diversity, MOESM4)"
curl -sL -o "$G/FungalTraits_Polme2020.xlsx" \
 "https://static-content.springer.com/esm/art%3A10.1007%2Fs13225-020-00466-2/MediaObjects/13225_2020_466_MOESM4_ESM.xlsx"
ls -la "$G/FungalTraits_Polme2020.xlsx"
echo "### FUNGuild script"
curl -sL -o "$G/Guilds_v1.1.py" "https://raw.githubusercontent.com/UMNFuN/FUNGuild/master/Guilds_v1.1.py"
ls -la "$G/Guilds_v1.1.py"
echo "### FUNGuild database endpoint reachable?"
for u in "http://www.stbates.org/funguild_db_2.php" "http://www.stbates.org/funguild_db.php"; do
  echo "  $(timeout 30 curl -sL -o /dev/null -w '%{http_code} %{size_download}' "$u" 2>/dev/null)  $u"
done
echo "### inspect FungalTraits sheets"
python - <<'PY'
import openpyxl
p="<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/guilds/FungalTraits_Polme2020.xlsx"
wb=openpyxl.load_workbook(p, read_only=True)
print("  sheets:", wb.sheetnames)
ws=wb[wb.sheetnames[0]]
rows=ws.iter_rows(min_row=1,max_row=3,values_only=True)
for i,r in enumerate(rows):
    print(f"  row{i}:", [str(x)[:26] for x in r[:12]])
PY
