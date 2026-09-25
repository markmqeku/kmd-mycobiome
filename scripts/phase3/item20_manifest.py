#!/usr/bin/env python
"""Final manifest, pointing at the corrected artefacts in clean/ where one exists and at the
original otherwise. Every entry is existence-checked."""
import os, csv
K = "<KMD_ROOT>/"
R = "__reanalysis_2026-06/phase3/outputs_27-07-2026"
C = f"{R}/clean"
FIGS = [
 ("Figure 1","Host DNA content by site",f"{R}/figures/Figure1_host_by_site.png","Results 3.1","unchanged"),
 ("Figure 2","Bloemfontein host load against fungal ASV count",f"{C}/figures/Figure2_bloem_hostload_vs_ASVs.png","Results 3.2","corrected"),
 ("Figure 3","Phylogenetic placement of the Tremellomycetes ASVs",f"{R}/figures/Figure3_tremellomycetes_placement.png","Results 3.3","unchanged, but 7 tips now withdrawn"),
 ("Figure 4","Taxonomic composition by site and condition",f"{C}/figures/Figure4_composition.png","Results 3.4","corrected"),
 ("Figure 5","Hill numbers q0, q1, q2 by site and condition",f"{C}/figures/Figure5_hill_numbers.png","Results 3.5","corrected, depth 57,124"),
 ("Figure 6","Bray-Curtis and Jaccard ordination",f"{C}/figures/Figure6_ordination.png","Results 3.6","corrected, depth 57,124"),
 ("Figure 7","Per-sample Curvibasidium relative abundance",f"{C}/figures/Figure7_curvibasidium.png","Results 3.7","corrected"),
 ("Figure 8","Placement of abundant unnamed ASVs by neighbourhood",f"{R}/figures/Figure8_placement_by_neighbourhood.png","Results 3.10","unchanged"),
 ("Figure S1","Rarefaction curves, 16 core samples",f"{C}/figures/FigureS1_rarefaction_core16.png","Methods 2.7","corrected, depth 57,124"),
 ("Figure S2","DADA2 read retention by stage",f"{R}/figures/FigureS2_dada2_retention.png","Methods 2.3","unchanged"),
 ("Figure S3","Taxonomy assignment rate by rank",f"{R}/figures/FigureS3_assignment_by_rank.png","Methods 2.5","unchanged"),
 ("Figure S4","Guild assignment coverage by tool",f"{C}/figures/FigureS4_guild_coverage.png","Results 3.9","corrected"),
]
TABS = [
 ("Table S1","ASV catalogue, 1,142 fungal ASVs, final taxonomy and representative sequences",f"{C}/asv_catalogue/ASV_catalogue.xlsx","Methods 2.8","corrected"),
 ("Table S2","Bloemfontein presence-absence inventory, 199 taxa, distinct ASVs and occurrences",f"{C}/TableS2_bloemfontein_inventory.tsv","Results 3.8","corrected"),
 ("Table S3","Taxonomy threshold comparison, n = 1,142",f"{C}/TableS3_taxonomy_thresholds.tsv","Methods 2.5","corrected"),
 ("Table S4","DADA2 statistics, all 45 samples",f"{R}/assembly/TableS4_dada2_stats.tsv","Methods 2.3","unchanged"),
 ("Table S5","Guild assignments per ASV, 1,142 ASVs",f"{C}/guild_assignments_per_ASV.tsv","Results 3.9","corrected"),
 ("Table S6","Placement results for abundant unnamed ASVs",f"{R}/placement/placement_results.tsv","Results 3.10","unchanged"),
 ("Table S7","Reference provenance for placement",f"{R}/placement/reference_provenance.tsv","Results 3.10","unchanged"),
 ("Table S8","Hill numbers per sample at depth 57,124",f"{C}/TableS8_hill_numbers.tsv","Results 3.5","corrected"),
 ("Table S9","Non-fungal screen: call, evidence and reason for all 277 screened ASVs",f"{R}/full_screen/screen_calls.tsv","Results 3.11","new"),
 ("Table S10","Screen basis for all 1,173 ASVs, with reason each was or was not sent to NCBI",f"{R}/full_screen/screen_basis.tsv","Methods 2.4","new"),
 ("Table S11","ASVs retained as inconclusive",f"{R}/full_screen/inconclusive_asvs.tsv","Results 3.11","new"),
 ("Data","GenBank submission file, 71 Tremellomycetes ASVs, unsubmitted",f"{C}/genbank_submission/Tremellomycetes_ASVs.fasta","Methods 2.8","corrected"),
 ("Data","SRA metadata for raw read deposition, 45 samples, unsubmitted",f"{R}/assembly/SRA_metadata.tsv","Methods 2.8","unchanged"),
 ("Data","Consolidated reference list with verification status",f"{C}/reference_list.tsv","References","corrected"),
 ("Data","Hill numbers, superseded against corrected, side by side",f"{C}/hill_old_vs_new.tsv","Results 3.5","new"),
 ("Data","Screen validation: UNITE composition and the 149-host-ASV control",f"{R}/full_screen/screen_validation.txt","Methods 2.4","new"),
 ("Data","Outstanding INSERT checklist",f"{C}/INSERT_checklist.tsv","(editorial)","corrected"),
]
out = f"{K}{C}/figure_table_manifest.tsv"
missing = []
with open(out, "w", newline="", encoding="utf-8") as fh:
    w = csv.writer(fh, delimiter="\t")
    w.writerow(["number","title","file","cited_in","status","exists"])
    for n, t, p, s, st in FIGS + TABS:
        ok = os.path.exists(K + p)
        if not ok: missing.append((n, p))
        w.writerow([n, t, p, s, st, "yes" if ok else "MISSING"])
print(f"{'item':11} {'status':16} {'exists':8} title")
for n, t, p, s, st in FIGS + TABS:
    ok = os.path.exists(K + p)
    print(f"{n:11} {st:16} {'yes' if ok else 'MISSING':8} {t[:56]}")
print(f"\n{len(FIGS)} figures, {len(TABS)} tables and data files -> {out}")
print("ALL PRESENT" if not missing else f"MISSING: {missing}")
