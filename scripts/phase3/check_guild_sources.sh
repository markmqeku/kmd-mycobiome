#!/usr/bin/env bash
URLS=(
"https://raw.githubusercontent.com/UMNFuN/FUNGuild/master/Guilds_v1.1.py"
"https://raw.githubusercontent.com/UMNFuN/FUNGuild/master/FUNGuild_db.json"
"https://raw.githubusercontent.com/traitecoevo/fungaltraits/master/inst/extdata/fungal_traits.csv"
"https://static-content.springer.com/esm/art%3A10.1007%2Fs13225-020-00466-2/MediaObjects/13225_2020_466_MOESM4_ESM.xlsx"
"https://api.github.com/repos/UMNFuN/FUNGuild"
"https://api.github.com/repos/traitecoevo/fungaltraits"
)
for u in "${URLS[@]}"; do
  code=$(timeout 30 curl -sL -o /dev/null -w "%{http_code} %{size_download}" "$u" 2>/dev/null)
  echo "$code  $u"
done
echo "--- local copies already in project? ---"
find <KMD_ROOT> -iname "*funguild*" -o -iname "*fungaltrait*" 2>/dev/null | head -20
