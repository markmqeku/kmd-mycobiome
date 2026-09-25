import csv
G = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026/guilds"
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
yeast = {r[0] for r in csv.reader(open(f"{O}/nontarget_screen/query_abundance.tsv"), delimiter="\t")
         if r and r[0] != "ASV_ID_md5" and r[2] == "1"}
rows = list(csv.DictReader(open(f"{G}/guild_assignments_per_ASV.tsv"), delimiter="\t"))
grand = sum(int(r["total_reads"]) for r in rows)
yr = [r for r in rows if r["ASV_ID_md5"] in yeast]
yg = [r for r in yr if r["FUNGuild_guild"]]
print(f"yeast ASVs: {len(yr)}   FUNGuild-assigned: {len(yg)}   reads assigned: {sum(int(r['total_reads']) for r in yg):,}")
print("\nwhat rank did FUNGuild match them at, and to what?")
from collections import Counter
c = Counter((r["FG_match_rank"], r["FUNGuild_guild"][:60]) for r in yg)
for (rk, g), n in c.most_common(5):
    print(f"   rank={rk:8} n={n:3}  guild={g}")
fg_reads_all = sum(int(r["total_reads"]) for r in rows if r["FUNGuild_guild"])
fg_reads_yeast = sum(int(r["total_reads"]) for r in yg)
ft_reads_all = sum(int(r["total_reads"]) for r in rows if r["FungalTraits_primary_lifestyle"])
ft_reads_yeast = sum(int(r["total_reads"]) for r in yr if r["FungalTraits_primary_lifestyle"])
print(f"\nFUNGuild read assignment INCLUDING spurious yeast : {100*fg_reads_all/grand:.1f}%")
print(f"FUNGuild read assignment EXCLUDING spurious yeast : {100*(fg_reads_all-fg_reads_yeast)/grand:.1f}%")
print(f"  -> {100*fg_reads_yeast/grand:.1f} percentage points of FUNGuild's coverage is the misassigned lineage")
print(f"FungalTraits read assignment (yeast contributes {100*ft_reads_yeast/grand:.1f}%): {100*ft_reads_all/grand:.1f}%")
