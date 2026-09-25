#!/usr/bin/env python
"""12b figure and table manifest; 12c supplementary currency audit."""
import os, csv, glob
O = "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
A = f"{O}/assembly"; os.makedirs(A, exist_ok=True)
R = "__reanalysis_2026-06/phase3/outputs_27-07-2026"

FIGS = [
 ("Figure 1","Host DNA content by site",f"{R}/figures/Figure1_host_by_site.png","Results 3.1"),
 ("Figure 2","Bloemfontein host load against fungal ASV count",f"{R}/figures/Figure2_bloem_hostload_vs_ASVs.png","Results 3.2"),
 ("Figure 3","Phylogenetic placement of the 78 Tremellomycetes ASVs",f"{R}/figures/Figure3_tremellomycetes_placement.png","Results 3.3"),
 ("Figure 4","Taxonomic composition by site and condition",f"{R}/figures/Figure4_composition.png","Results 3.4"),
 ("Figure 5","Hill numbers q0, q1, q2 by site and condition",f"{R}/figures/Figure5_hill_numbers.png","Results 3.5"),
 ("Figure 6","Bray-Curtis and Jaccard ordination",f"{R}/figures/Figure6_ordination.png","Results 3.6"),
 ("Figure 7","Per-sample Curvibasidium relative abundance",f"{R}/figures/Figure7_curvibasidium.png","Results 3.7"),
 ("Figure 8","Placement of abundant unnamed ASVs by neighbourhood",f"{R}/figures/Figure8_placement_by_neighbourhood.png","Results 3.10"),
 ("Figure S1","Rarefaction curves, 16 core samples",f"{R}/figures/FigureS1_rarefaction_core16.png","Methods 2.7"),
 ("Figure S2","DADA2 read retention by stage",f"{R}/figures/FigureS2_dada2_retention.png","Methods 2.3"),
 ("Figure S3","Taxonomy assignment rate by rank",f"{R}/figures/FigureS3_assignment_by_rank.png","Methods 2.5"),
 ("Figure S4","Guild assignment coverage by tool",f"{R}/figures/FigureS4_guild_coverage.png","Results 3.9"),
]
TABS = [
 ("Table S1","ASV catalogue, 1,173 fungal ASVs, final taxonomy and representative sequences",f"{R}/asv_catalogue_final/ASV_catalogue.xlsx","Methods 2.8"),
 ("Table S2","Bloemfontein presence-absence inventory, 192 taxa, with per-sample host load",f"{R}/assembly/TableS2_bloemfontein_inventory.tsv","Results 3.8"),
 ("Table S3","Taxonomy threshold comparison, UNITE 8.2 / 10.0 permissive / 10.0 rank-threshold",f"{R}/assembly/TableS3_taxonomy_thresholds.tsv","Methods 2.5"),
 ("Table S4","DADA2 statistics, all 45 samples",f"{R}/assembly/TableS4_dada2_stats.tsv","Methods 2.3"),
 ("Table S5","Guild assignments per ASV with confidence rankings",f"{R}/guilds/guild_assignments_per_ASV.tsv","Results 3.9"),
 ("Table S6","Placement results for abundant unnamed ASVs",f"{R}/placement/placement_results.tsv","Results 3.10"),
 ("Table S7","Reference provenance for placement, accession, name and type status",f"{R}/placement/reference_provenance.tsv","Results 3.10"),
 ("Table S8","Hill numbers per sample at depth 57,141",f"{R}/gap_closure/hill_numbers_per_sample.tsv","Results 3.5"),
 ("Data","GenBank submission file, 78 Tremellomycetes ASVs, unsubmitted",f"{R}/genbank_submission/Tremellomycetes_ASVs.fasta","Methods 2.8"),
 ("Data","SRA metadata for raw read deposition, 45 samples, unsubmitted",f"{R}/assembly/SRA_metadata.tsv","Methods 2.8"),
 ("Data","Consolidated reference list with verification status",f"{R}/assembly/reference_list.tsv","References"),
]
base = "<KMD_ROOT>/"
with open(f"{A}/figure_table_manifest.tsv","w",newline="",encoding="utf-8") as fh:
    w=csv.writer(fh,delimiter="\t"); w.writerow(["number","title","file","cited_in","exists"])
    for n,t,p,s in FIGS+TABS:
        w.writerow([n,t,p,s,"yes" if os.path.exists(base+p) else "MISSING"])
print(f"{'item':11} {'exists':7} title")
missing=[]
for n,t,p,s in FIGS+TABS:
    ok = os.path.exists(base+p)
    if not ok: missing.append((n,p))
    print(f"{n:11} {'yes' if ok else 'MISSING':7} {t[:64]}")
print(f"\nmanifest -> {A}/figure_table_manifest.tsv")
if missing:
    print("\nSTILL TO BUILD:")
    for n,p in missing: print(f"   {n}: {p}")
