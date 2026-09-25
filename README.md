# Fungal communities of karee malformation disease in *Searsia lancea*: analysis code and data

Code, parameters, phylogenetic placement files and summary data for:

> Mqeku, M., Cason, E., Kinge, R.T., Gryzenhout, M. Fungal communities of karee malformation disease in *Searsia lancea*: an abundant unresolved Tremellomycetes lineage and the limits of ITS2 metabarcoding. Submitted to Fungal Ecology.

Raw reads are in the NCBI Sequence Read Archive and the lineage sequences in GenBank (accessions given in the paper),
not here.

## Licences

- **Code** (`scripts/`): MIT licence, in `LICENSE`.
- **Data** (`data/`, `placement/`, `parameters/`): Creative Commons Attribution 4.0 International (CC BY 4.0), in `LICENSE-DATA`.

## Contents

| folder or file | contents |
|---|---|
| `scripts/phase0` to `scripts/phase3` | the analysis, in order: raw-read inventory, primers and DADA2 trials, denoising, then host and non-fungal screening, taxonomy, placement, guilds, rarefaction, diversity and figures |
| `parameters/manifests` | the three QIIME 2 import manifests: which read files entered denoising for which library |
| `placement/yeast_placement` | the first placement of the 78 ASVs first delimited as the lineage |
| `placement/yeast_placement_71` | the rerun on the 71 retained ASVs (Figure 3): input, MAFFT alignment (1,949 columns, 124 sequences), IQ-TREE tree, per-ASV placement |
| `placement/placement` | per-neighbourhood reference sets, alignments and trees for the abundant unnamed ASVs (Figure 8), with the NCBI organism name of every reference accession |
| `data/FINAL_NUMBERS_27-07-2026.md` | every value reported in the paper, with its source file and a change record |
| `data/metadata_v3_27-07-2026.tsv`, `data/TableS12_library_metadata.tsv` | library metadata: site, condition, tissue, age class, trees pooled, sequencing run, reads, inclusion |
| `data/TableS1_ASV_catalogue.tsv` | the 1,142 fungal ASVs with taxonomy, reads, prevalence and, for the lineage, the GenBank clone name |
| `data/KMD_ITS2_71_ASVs.fasta` | all 71 lineage ASVs as deposited in GenBank, labelled uncultured Tremellomycetes (67) or uncultured fungus (4) |
| `data/CLONE_MAP.tsv`, `data/KMD_ITS2_71_ASVs_source_modifiers.tsv` | clone name to ASV hash, the two tests and the organism; the GenBank source modifiers |
| `data/verified_ontology_and_taxonomy_terms.json` | the ENVO, Plant Ontology and NCBI Taxonomy terms used in the deposits, as verified |
| `data/r4_*.json`, `data/r6c_clade_summary.json`, `data/r6a_depth_reconciliation/` | Didymellaceae and Mycosphaerellaceae shares; placement-tree summary; ordination and sharing at the locked depth |

Scripts ran with the project folder as their data store. In these copies absolute paths are replaced by `<KMD_ROOT>`,
`<USER_HOME>`, `<HOME>` and `<TEMP>`, and e-mail addresses by `<email>`; set them to your own locations before running.
Script names carry the task number of the report or prompt item that used them.

## Pipeline, in order

1. **Libraries.** 45 usable sequencing libraries: Bloemfontein 29 (one tree per library), Christiana 9 and Pretoria 7 (each pooling tissue from four trees of the same site and condition). P6, P15 and P24 were truncated in delivery, YT-P was corrupted and 30-C was never used; OT-P was salvaged by removing one corrupted record, leaving 110,302 pairs.
2. **Primers.** ITS3 and ITS4 removed with cutadapt 4.9 in QIIME 2 2024.10.1. ITSxpress was not used, because it destroyed the 2 x 250 bp libraries.
3. **Denoising.** DADA2 1.30.0 run separately per sequencing run and the tables merged: run 3 (Bloemfontein, 2 x 250 bp) truncated at 228 and 200; run 13 (Christiana and Pretoria, 2 x 301 bp) at 280 and 230. ASV identifiers are MD5 hashes of the sequences.
4. **Host screen.** 149 ASVs matching *Searsia lancea* removed (1,487,006 reads); raw reads screened directly to confirm the host content of the delivered libraries. Host and lineage amplicon lengths rule out size selection as the cause of the host disparity (`scripts/phase3/r0c_size_selection.py`).
5. **Non-fungal screen.** 896 ASVs accepted as fungal by alignment to UNITE; 277 compared with NCBI nt by blastn; 31 non-fungal ASVs removed. Result: 1,142 fungal ASVs, 2,517,177 reads.
6. **Taxonomy.** VSEARCH consensus against the UNITE QIIME release for Fungi 10.0 of 4 April 2024 (doi:10.15156/BIO/2959336), truncated by rank-specific identity (species 97, genus 95, family 90, order 85 percent); UNITE 8.2 as a control.
7. **Placement.** The lineage (71 ASVs, `scripts/phase3/r6c_fig3_rerun.sh`) and 28 abundant unnamed ASVs placed with MAFFT and IQ-TREE 2 (GTR+G, 1,000 ultrafast bootstrap replicates). The lineage is assigned to Tremellomycetes only provisionally: the tree support is weak (44 for the clade holding the Tremellomycetes references).
8. **Guilds.** FungalTraits and FUNGuild; supplementary only.
9. **Diversity.** 16 core libraries rarefied to 57,124 reads (seed 42); Hill numbers with iNEXT and hillR in R; Bray-Curtis and Jaccard ordination. No inferential testing.
10. **Figures.** Drawn by the scripts named in `data/FINAL_NUMBERS_27-07-2026.md` and exported without title lines at 600 dpi by `scripts/phase3/p5_regen_figures.py`.
11. **Deposits.** BioSample terms and the GenBank set are built by `scripts/phase3/p5_biosample.py` and `p5_genbank.py`.

## Citation

See `CITATION.cff`. Cite the paper and this archive's DOI.
