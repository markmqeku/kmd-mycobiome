# FINAL_NUMBERS — KMD fungal ITS2 re-analysis
**Frozen 27-07-2026. This table is the SOLE source for the manuscript. ANALYSIS CLOSED.**
Every value traced to its source file. Paths relative to `KMD/__reanalysis_2026-06/`.
Root: `phase3/outputs_27-07-2026/` = `OUT`.

---

## 1. Pipeline and software
| item | value | source |
|---|---|---|
| Platform | Illumina MiSeq, instrument M03998 | raw FASTQ headers |
| Marker / primers | ITS2; ITS3 + ITS4 (White et al. 1990) | verified in reads, `phase0/PHASE0_INVENTORY_INTEGRITY_2026-06-17.md` |
| Primer detection rate | R1 99.8%, R2 99.1% (sample 25-C) | `phase1/PHASE1_PRIMER_DADA2_2026-06-17.md` |
| QIIME 2 | 2024.10.1 (q2cli 2024.10.1, Python 3.10.14, DADA2 1.30.0, cutadapt 4.9) | `phase2/RUN_LOG_methods.md` |
| Read lengths | Bloemfontein 2×250; Christiana + Pretoria 2×301 | `phase0/` inventory |
| DADA2 group A (Bloemfontein) | `--p-trunc-len-f 228 --p-trunc-len-r 200` | `phase2/run_P2c3.sh` |
| DADA2 group B (Christiana + Pretoria) | `--p-trunc-len-f 280 --p-trunc-len-r 230` | `phase2/run_P2c3.sh` |
| ITSxpress | **not used** (destroyed the 2×250 libraries: 25,000 → 51/798/98 pairs) | `phase1/PHASE1_PRIMER_DADA2_2026-06-17.md` |
| Taxonomy | UNITE **10.0**, fungi, dynamic (93,085 refs); `classify-consensus-vsearch` | `OUT/unite_current/`, `OUT/taxonomy_unite10/` |
| Rank cutoffs (fallback) | species ≥97, genus ≥95, family ≥90, order ≥85, floor 80 | `phase3/G4_G5_hybrid.py` |
| `classify-sklearn` | **not feasible** — training OOM-killed (7.3 GB resident, 29 GB virtual) | `OUT/taxonomy_sklearn/`, `phase3/G4.log` |

## 2. Reads and ASVs
| quantity | value | source |
|---|---|---|
| Total reads, all 45 samples (incl. host) | **4,010,617** | `OUT/table_exp/feature-table.tsv` |
| Total ASVs before host removal | **1,322** | `phase2/merged2/rep_exp/dna-sequences.fasta` |
| Host (*Searsia lancea*) ASVs removed | **149** | `OUT/nontarget_screen/host_asv_ids_ALL.txt` |
| Host reads removed | **1,487,006 (37.1%)** | `OUT/nontarget_screen/`, `phase3/host_load_by_sample.py` |
| **Fungal ASVs (final table)** | **1,142** | `OUT/clean/feature-table-clean.tsv` |
| **Fungal reads (final table)** | **2,517,177** | `OUT/clean/feature-table-clean.tsv` |
| *superseded* | 1,173 ASVs / 2,523,611 reads | after host removal only. **31 ASVs (6,434 reads, 0.255%) were later confirmed non-fungal by a whole-table NCBI screen and removed**: 15 land plants, 6 green algae, 10 contested green algae. 21 inconclusive ASVs (276 reads) retained. Source `OUT/full_screen/screen_calls.tsv` |
| Host identity evidence | 99.2–99.7% id, 98% qcov to **AY641514.1 *Searsia lancea*** | `OUT/nontarget_screen/blast_host.tsv` |
| Host load by site | Bloemfontein **84.3%**; Christiana **0.0%**; Pretoria **0.0%** | `phase3/host_load_by_sample.py` |

**Note:** the ~643,000 reads "recovered" by the truncation fix in the seven previously zeroed Bloemfontein samples were almost entirely host DNA (P2 100% host → 0 fungal reads; only P4 was a genuine fungal recovery). Do not present that recovery as a fungal-ecology result.

## 3. Sample set and per-group n
| set | n | source |
|---|---|---|
| Samples sequenced (design) | 49 | `phase0/metadata_consolidated.tsv` |
| Write-offs (corrupt/truncated raw files) | 4 — P6, P15, P24, YT-P | `phase0/PHASE0.5_SOURCE_SCAN_2026-06-17.md` |
| Excluded (unmappable/unrelated) | 6 — CR, PR, S1–S4 | same |
| Denoised successfully | 45 | `phase2/*/stats_exp/stats.tsv` |
| **Bloemfontein — EXCLUDED from quantitative analysis** | 29 (inventory only) | see §7 |
| **CORE ANALYSIS SET** | **16** = Christiana 9 + Pretoria 7 | `OUT/metadata_v2_27-07-2026.tsv` |

| group | n |
|---|---|
| Christiana asymptomatic | **3** (25-C, 26-C, 27-C) |
| Christiana symptomatic | **6** (OL-C, OM-C, OT-C, YL-C, YM-C, YT-C) |
| Pretoria asymptomatic | **2** (28-C, 29-C) |
| Pretoria symptomatic | **5** (OL-P, OM-P, OT-P, YL-P, YM-P) |

*condition_std*: `symptomatic` = malformed tissue; `asymptomatic` = normal tissue on a KMD-symptomatic tree; `reference_site` = Bloemfontein. **"Healthy" is not used.**

## 4. Rarefaction depth
| item | value | source |
|---|---|---|
| **Primary depth** | **57,124** | `phase3/item17b_hill.R` |
| Justification | **the deepest level retaining all 16 core samples** (set by the lowest sample, 27-C = 57,124 on the corrected table; the next-smallest core sample is 26-C at 64,674) | same |
| *superseded* | 57,141 | 27-C's total before 17 non-fungal reads were removed from it. At 57,141 the corrected table retains only 15/16. Re-locked per DEC-1 |
| **Not used as justification** | sample coverage. DADA2 emits no singletons, so Good's/iNEXT coverage = 1.0000 **by construction** and is uninformative | I-1 |
| **Not reported** | Chao1, ACE, iNEXT extrapolation — all depend on singleton/doubleton frequencies destroyed by DADA2 | I-1 |

**Methods limitation to state:** because DADA2 removes singletons, coverage-based estimators and asymptotic richness estimators are not interpretable for this dataset. Only observed richness at a fixed depth and Hill q1/q2 are reported. The iNEXT bootstrap intervals below are correspondingly **narrower than the true uncertainty** and should be read as within-sample resampling ranges, not as biological confidence intervals.

## 5. Diversity — descriptive only, no significance testing (every cell n = 2–6)
Observed richness at depth 57,124 (q0) and Hill q1/q2. Source: `phase3/item17b_hill.R`, `OUT/clean/hill_numbers_per_sample.tsv`.

| group | n | q0 (observed richness) | q1 exp(Shannon) | q2 inverse Simpson |
|---|---|---|---|---|
| Christiana asymptomatic | 3 | 78.4 (64.0–96.3) | 7.6 (5.6–9.1) | 4.4 (3.0–5.9) |
| Christiana symptomatic | 6 | 106.5 (80.8–130.5) | 14.2 (9.7–19.5) | 8.7 (5.8–11.8) |
| Pretoria asymptomatic | 2 | 143.5 (136.7–150.3) | 17.3 (14.7–19.8) | 8.0 (6.2–9.8) |
| Pretoria symptomatic | 5 | 116.0 (88.9–199.5) | 14.4 (11.1–16.9) | 8.3 (5.9–10.0) |

*Superseded* (depth 57,141, uncorrected table): q0 means 80.1 / 107.6 / 143.5 / 116.7. Only q0 moved materially; q1 and q2 group means changed by at most 0.01. Largest single-sample change −2.93 (q0). Side by side in `OUT/clean/hill_old_vs_new.tsv`.

Ranges are across samples within the group. **Direction is inconsistent between sites**: symptomatic higher at Christiana, asymptomatic higher at Pretoria (q0 and q1). Unchanged by the correction.

**Ordination** (`phase3/item17b_ord_clean.py` run at the locked depth 57,124, seed 42, corrected table; `phase3/r6a_ordination_57124.sh`, output `OUT/clean/r6a/ord_depth57124.txt`): PCoA axis 1 separates **site**, not condition. Bray-Curtis 63.1%, Jaccard **34.1%**. *Superseded: Jaccard 34.4%, which came from the same script with its `DEPTH` line still at the superseded 57,141 (Figure 6 was drawn at 57,124 and showed 34.1%); Bray-Curtis is 63.1% at both depths. Changed 25-09-2026, see change record.* all Christiana negative, all Pretoria positive, conditions interleaved within site. No PERMANOVA, no p-values.

> **Flagged, not silently corrected.** This table previously recorded Bray-Curtis 33.1%. That value **does not reproduce from its own cited source**: running the unmodified `H0_H4b.py` on the unmodified table gives **63.1%**, not 33.1%. The discrepancy therefore predates the non-fungal correction and is not caused by it. Jaccard reproduces (33.7% before, 34.4% after). The 33.1% figure needs an explanation before this row is relied on.

**Christiana position at 57,124** (`phase3/item23_position.py`, output `OUT/clean/r6a/item23_position_57124.txt`, and `OUT/clean/r6a/ord_depth57124.txt`): Bray-Curtis axis 1 for 25-C, 26-C, 27-C = **-0.428, -0.421, -0.318**; the six symptomatic libraries span **-0.381 to -0.117**. Largest difference in mean dissimilarity to the asymptomatic against the symptomatic group: **0.032** Bray-Curtis (0.0324, 26-C), **0.028** Jaccard (0.0276, 25-C); mean within-symptomatic Bray-Curtis dissimilarity **0.376** (0.3757). *Superseded (depth 57,141): -0.427, -0.420, -0.317; -0.381 to -0.121; 0.033. Changed 25-09-2026, see change record.*

## 6. Effort-corrected ASV sharing (replaces the naive 81% vs 43–48%)
Larger group subsampled without replacement to the smaller group's n, 1,000 iterations, on tables rarefied to equal depth. Source: `phase3/item17b_ord_clean.py` at depth 57,124 (corrected table; `OUT/clean/r6a/ord_depth57124.txt`).

| site | naive (unequal effort) | effort-corrected | difference |
|---|---|---|---|
| Christiana | asym 141 (49%) vs sym 240 (83%) | asym **141.0** vs sym **179.6** [155 to 204] | **+38.6** |
| Pretoria | asym 210 (60%) vs sym 261 (75%) | asym **210.0** vs sym **167.4** [115 to 230] | **-42.6** |

*Superseded* (uncorrected table): Christiana asym 143 vs sym 181.1, difference +38.1; Pretoria asym 210 vs sym 165.6, difference −44.4. *Superseded* (corrected table at the superseded depth 57,141, the values this table carried until 25-09-2026): Christiana 141.0 vs 180.4 [155 to 204], +39.4, naive 141 vs 241; Pretoria 209.0 vs 164.5 [115 to 222], -44.5, naive 209 vs 254. The direction and the conclusion are unchanged.

**The naive contrast was largely a sampling-effort artifact and the direction reverses at Pretoria.** The claim "malformed tissues harbour 81–82% of ASVs versus 43–48%" is not supportable.

## 7. Composition at the 95% genus threshold (threshold-robust: identical at 95/93/90)
Source: `phase3/H3_H5_H7.py`.

| Christiana asymptomatic (n=3) | % | Christiana symptomatic (n=6) | % |
|---|---|---|---|
| Tremellomycetes yeast lineage | 41.5 | Tremellomycetes yeast lineage | 39.6 |
| Dothideomycetes (class) | 25.2 | *Cladosporium* | 11.0 |
| Didymellaceae (family) | 18.8 | Dothideomycetes (class) | 10.0 |
| Saccotheciaceae (family) | 5.8 | Saccotheciaceae (family) | 8.7 |
| *Aureobasidium* | 1.8 | *Curvibasidium* | 5.8 |

| Pretoria asymptomatic (n=2) | % | Pretoria symptomatic (n=5) | % |
|---|---|---|---|
| Dothideaceae (family) | 26.7 | Dothideaceae (family) | 29.3 |
| Tremellomycetes yeast lineage | 14.6 | Tremellomycetes yeast lineage | 21.5 |
| Didymellaceae (family) | 12.5 | Didymellaceae (family) | 12.2 |
| Saccotheciaceae (family) | 10.2 | *Curvibasidium* | 11.3 |
| *Cladosporium* | 8.4 | *Cladosporium* | 10.3 |

### Historical claims — status
| taxon | result | verdict |
|---|---|---|
| ***Didymella*** | 0.00–0.01% in all four groups | **NOT REPRODUCIBLE.** Reads resolve only to Didymellaceae (family). The "64% → 10%" collapse cannot be recovered at any threshold. |
| ***Mycosphaerella*** | 0.00% in all four groups | **NOT REPRODUCIBLE.** The "3% → 22%" increase does not survive. |
| ***Curvibasidium*** | see §8 | **Partially supported — site-dependent.** |

### Didymellaceae and Mycosphaerellaceae, exact counts (R4, added 25-09-2026)
Source: `phase3/r4_didymella.py`, output `phase3/r4_didymella.json`. Fungi-only table, rank-specific thresholds, reads pooled within each group, not rarefied.

| group | n | fungal reads | Didymellaceae reads (share) | *Didymella* reads | Mycosphaerellaceae reads | *Mycosphaerella* reads |
|---|---|---|---|---|---|---|
| Christiana asymptomatic | 3 | 208,618 | 39,138 (18.76%) | 0 | 0 | 0 |
| Christiana symptomatic | 6 | 519,642 | 28,114 (5.41%) | 0 | 0 | 0 |
| Pretoria asymptomatic | 2 | 203,001 | 25,345 (12.49%) | 19 | 0 | 0 |
| Pretoria symptomatic | 5 | 1,315,622 | 160,528 (12.20%) | 0 | 9 | 0 |

*Didymella* reads: 19, from one ASV (0e16bd66), Pretoria asymptomatic only. Mycosphaerellaceae: 9 reads from two ASVs, Pretoria symptomatic only. Bloemfontein, presence only (28 libraries with fungal reads): Didymellaceae in 14, Mycosphaerellaceae in 1, *Didymella* and *Mycosphaerella* in 0. These counts replace the rounded "0.00 percent" wording.


## 8. *Curvibasidium* per-sample (I-2) — 10 ASVs at ≥95% identity
Source: `phase3/I2_curvibasidium.py`.

| site | condition | per-sample relative abundance | detected | ≥1% |
|---|---|---|---|---|
| Christiana | asymptomatic (n=3) | 0.33, 0.24, 0.10 | 3/3 | **0/3** |
| Christiana | symptomatic (n=6) | 20.23, 6.10, 5.25, 2.27, 1.00, 0.27 | 6/6 | **5/6** |
| Pretoria | asymptomatic (n=2) | 6.24, 3.03 | 2/2 | 2/2 |
| Pretoria | symptomatic (n=5) | 19.16, 12.19, 10.73, 2.99, 2.71 | 5/5 | 5/5 |

**Medians (corrected J-0, 28-07-2026).** For an even-numbered group the median is the mean of the two middle values.
- Christiana asymptomatic (n=3): 0.33, 0.24, 0.10 → **median 0.24%**
- Christiana symptomatic (n=6): 0.27, 1.00, 2.27, 5.25, 6.10, 20.23 → **median 3.76%** (= (2.27 + 5.25)/2). **Contrast ≈ 16-fold**, not 20-fold.
- Pretoria asymptomatic (n=2): 3.03, 6.24 → **median 4.63%**, from the unrounded values (2,603 of 85,895 and 7,304 of 117,106 reads; 4.634%). *Superseded 25-09-2026: 4.64%, the mean of the rounded values (D0d, change record).*
- Pretoria symptomatic (n=5): 2.71, 2.99, 10.73, 12.19, 19.16 → **median 10.73%**

> *Correction note: the earlier entry gave the Christiana symptomatic median as 5.25% and the contrast as ~20-fold. That arose from taking the upper-middle value of six rather than the mean of the two middle values. Text edit only; no recomputation of underlying data.*

**Verdict — must be reported with this caveat:**
**Enrichment in symptomatic tissue is consistent in DIRECTION at both sites, with clear separation at Christiana and overlapping ranges at Pretoria; magnitude at Christiana is influenced by one sample at 20.23%; n = 2 to 6; no test performed.**

Supporting detail: at Christiana no asymptomatic sample exceeds 0.35% and 5 of 6 symptomatic samples exceed 1%, though one symptomatic sample (OL-C, 0.27%) falls inside the asymptomatic range. At Pretoria two of five symptomatic samples (2.71%, 2.99%) fall below both asymptomatic values (3.03%, 6.24%).

## 9. Tremellomycetes yeast lineage — evidence
Source: `phase3/H2a_corrected.py`, `phase3/H2c_dates.py`, `OUT/yeast_placement/`.

| item | value |
|---|---|
| ASVs | **71** (retained as fungi; not host). *Superseded: 78; seven were confirmed green algae by the whole-table NCBI screen and removed, see section 11.* |
| Share of fungal reads | Christiana ~40%, Pretoria 14.6–21.5%. *The Bloemfontein share (15.2%) is not reported: Bloemfontein is presence and absence only (R5, 25-09-2026).* |
| Detection (R5, added 25-09-2026) | **24 of the 29** usable Bloemfontein libraries (absent from P1, P2, P12, P17, P27); 16 of 16 core libraries; 40 of 45 overall. Source `phase3/r2h_table_s12.py` |
| Nearest **named** Tremellomycetes reference | *Filobasidium* PP844358.1, **85.1–85.5% over 276 bp** (73% qcov) |
| Across all 74 with a hit | identity median **89.2%** (78.1–98.2); alignment median 121 bp; 17/74 align <100 bp |
| Nearest NCBI nt match | ***Cryptococcus* sp. isolate OTU796, MK019123.1 — 98.3–99.2% over 363 bp (98% qcov)** — an **unnamed environmental sequence** |
| MK019123.1 deposited | **2019-11-01** |
| UNITE 10.0 dataset DOI | **10.15156/BIO/2959336** ("UNITE QIIME release for Fungi", UNITE Community, 2024). **RESOLVED.** Identified from the RESCRIPt `UNITE_DOIS` table in `get_unite.py`, which maps version 10.0 + taxon_group fungi + singletons False to this DOI; those are exactly the parameters recorded in the artifact provenance. DataCite gives created **2024-04-04**, registered 2024-04-05. **Note:** the date 2024-04-04, withdrawn under DEC-19 as unsourced, is thereby confirmed and now has a source. Creators: Abarenkov, Zirk, Piirmann, Pohonen, Ivanov, Nilsson, Koljalg |
| Conclusion | the near-identical reference predates UNITE 10.0 by ~4.5 years → **a genuine UNITE coverage gap, not a lag** |
| UNITE's assignment | **wrong**: Agaricales (*Homophron*, a psathyrellaceous agaric) via a distant 87–88% match — rank truncation cannot detect a wrong-lineage match |
| Reported name | **"Tremellomycetes yeast lineage (ITS2-unresolved)"**. Genus NOT assigned; **LSU D1/D2 required** |
| Nomenclature checked | Index Fungorum: *Cryptococcus* has two homonyms (Kütz., Vuill.); *C. aciditolerans* → *Goffeauzyma aciditolerans*. NCBI's "*Cryptococcus* sp." is a legacy label, not used |
| 4 ASVs with no Tremellomycetes hit | 102, 9, 4, 4 reads; no confident NCBI hit either |
| **Tremellomycetes lineage placement, rerun tree (R6c, added 25-09-2026)** | **70 of 71 ASVs** fall within the clade containing all 48 Tremellomycetes references, in the tree rerun without the seven green-algal ASVs (`phase3/r6c_fig3_rerun.sh`; `OUT/yeast_placement_71/`), rooted on the four *Curvibasidium* references (Microbotryomycetes), none of which falls inside that clade. Ultrafast bootstrap support for the clade: **44**. Outside: ASV_ac534075 (2 reads). 32 ASVs carrying 646,258 of 647,397 lineage reads (**99.82%**) form a clade with MK019123.1 (ultrafast bootstrap **74**), including the most abundant ASV. Alignment: 1,949 columns, 124 sequences (71 ASVs, 48 Tremellomycetes and 4 *Curvibasidium* references, MK019123.1); IQ-TREE GTR+G, 1,000 ultrafast bootstrap replicates, seed 276583. The original 78-ASV tree under the same rooting and criterion: 70 of 71 retained ASVs and 7 of 7 green-algal ASVs inside (clade support 21). Source `phase3/r6c_clade_summary.py`, `OUT/yeast_placement_71/r6c_clade_summary.json`. |
| **Tremellomycetes reference hit (criterion clarified 25-09-2026)** | **68 of 71 ASVs have a blastn hit to the Tremellomycetes reference set.** *This row was titled "68 of 71 ASVs placed within Tremellomycetes"; the 68 is a reference-hit count, not a tree placement.* All 71 retained sequences match `OUT/yeast_placement/yeast_asvs.fasta` by exact sequence; 68 have a hit in `OUT/yeast_placement/yeast_vs_tremello.tsv` and 3 do not (`ASV_0505492f`, 102 reads; `ASV_5052240a`, 9 reads; `ASV_010fdafd`, 4 reads). The fourth no-hit ASV above, `ASV_e453cdd3` (4 reads), was among the seven green algae removed. Recomputed 23-09-2026 by an exact-sequence match of the 71 retained sequences against the placement input and hit table (see `KMD_CHECKPOINT.md` 2b). *Added 24-09-2026, see change record.* |
| Amplicon length, size-selection test (R0c, added 25-09-2026) | host ASVs **372 bp** (143 of 149 ASVs; 100% of host reads within 10 bp of 372); lineage ASVs **367 to 368 bp**, read-weighted median 368, **99.8%** of lineage reads within 10 bp of 372 (100% at Christiana and Pretoria). Source `phase3/r0c_size_selection.py` |
| **UNITE release** | **Fungal QIIME release 10.0, 4 April 2024, DOI 10.15156/BIO/2959336** ("UNITE QIIME release for Fungi", version 04.04.2024). Sources: the DOI landing page https://doi.org/10.15156/BIO/2959336 as recorded in `COWORK_FINDINGS/REFERENCES_COWORK.md` 2.1; the RESCRIPt 2024.10.0 source `rescript/get_unite.py`, `UNITE_DOIS["10.0"]["fungi"][False]`; the artifact provenance of `OUT/unite_current/unite10_seqs.qza` (version 10.0, taxon_group fungi, cluster_id dynamic, singletons false); and the UNITE repository row for this DOI, 18,895 RefS + 74,190 RepS = 93,085, the count in the artifact. *Added 24-09-2026, see change record.* |

**Withdrawn:** my earlier "100% identity to *Filobasidium*" (G-3) was an artifact of ranking hits by identity — it was 100% over **31 bp**. Superseded by the bitscore-ranked values above.

## 10. Bloemfontein — inventory only
Source: `phase3/H3_H5_H7.py`.

| item | value |
|---|---|
| Status | **EXCLUDED from all quantitative community analysis** |
| Reason | apparent richness tracks host load inversely: P10/P11/P12 ~100% host → 4–7 ASVs; P4/P5 <1% host → 313/314 ASVs. A technical gradient, not biological |
| Low-host samples used for inventory (host ≤10%) | **6** — P5 (0.10%), P32 (0.44%), P4 (0.92%), P18 (3.24%), P19 (3.33%), P25 (5.33%) |
| **Distinct taxa recorded** | **199** (presence/absence only), across **576 distinct ASVs** |
| Most frequent — **distinct ASVs** (occurrences) | Tremellomycetes yeast lineage 45 (70), unassigned 25 (31), Fungi (kingdom) 24 (28), Ascomycota 19 (26), Pleosporales 18 (27), *Aureobasidium* 14 (37), Sordariales 13 (17), Chaetomiaceae 12 (18), *Fusarium* 11 (17), Didymellaceae 10 (18) |

**Two corrections to this section, both recorded rather than silently applied:**

1. *Superseded taxon count:* **192**. The sample set is unchanged and provably identical (`P5, P32, P4, P18, P19, P25`); both counts were recomputed in one run to prove it (`phase3/item14_reconcile.py`). The 192 was computed on raw UNITE 10.0 with a flat 95/85 identity cut; the manuscript reports the rank-threshold final taxonomy throughout, which gives **199**. The cause is **not** the yeast relabelling, which is identical under both bases: 36 previously unassigned ASVs now resolve and 27 change rank (`phase3/item14b_cause.py`).
2. *Superseded "Most frequent" row:* the previous counts (unassigned 88, yeast lineage 78, *Aureobasidium* 37, …) were **(sample, ASV) occurrences labelled as ASVs**, inherited from H-5. Every entry overstated the ASV count. Both figures are now given.

Counts above are on the **corrected** table; the 31 non-fungal ASVs are removed. P4/P5 ASV counts fell from 329/326 to 313/314 as a result.
| Restriction | no abundances, no diversity metrics, no cross-site comparison (differing host load, run and chemistry) |

## 11. GenBank submission material (prepared, NOT submitted)
**Current (KMD PROMPT 5, D2, 25-09-2026):** `NCBI_SUBMISSION/GenBank/KMD_ITS2_71_ASVs.fasta`, **all 71 lineage ASVs**, rRNA-ITS Eukaryote workflow, clone names KMD_ITS2_ASV001 to KMD_ITS2_ASV071 by total read count (`NCBI_SUBMISSION/GenBank/CLONE_MAP.tsv`; GenBank help desk case CAS-1702031-X8R7H0: identifiers in "clone", not "isolate"). **67** pass both tests (blastn hit to the Tremellomycetes reference set, and inside the Tremellomycetes reference clade of the rerun tree): organism **uncultured Tremellomycetes**, taxid 543743. **4** pass one test only: organism **uncultured fungus**, taxid 175245 (KMD_ITS2_ASV030 = ASV_0505492f, 102 reads, tree only; KMD_ITS2_ASV054 = ASV_5052240a, 9 reads, tree only; KMD_ITS2_ASV060 = ASV_010fdafd, 4 reads, tree only; KMD_ITS2_ASV070 = ASV_ac534075, 2 reads, blastn only). Taxids verified in NCBI Taxonomy (`phase3/p5_verified_terms.json`). Source: `phase3/p5_genbank.py`.

*Superseded (Prompt 3 and 4 material, retired to `_retired/NCBI_GenBank_pre-prompt5/`):* the earlier GenBank material described below, including the 68-sequence file of Prompts 3 and 4.
`OUT/clean/genbank_submission/Tremellomycetes_ASVs.fasta` (**71 sequences**) and `Tremellomycetes_ASV_metadata.tsv`. Organism = **"Tremellomycetes sp."**, host *Searsia lancea*, ITS2, ITS3/ITS4, South Africa; identification note records the UNITE 10.0 gap and the Agaricales misassignment. **No genus assigned.**

*Superseded:* **78 sequences**. Seven were confirmed green algae by the whole-table NCBI screen and withdrawn: KMD_TremASV027, 037, 047, 048, 059, 066 and 069. All seven occur only in P4 and P5, and carry 257 reads combined. The clearest case is **KMD_TremASV066**, assigned by UNITE 10.0 to *Homophron spadiceum* at **81.0% identity with reported confidence 1.0** (its only hit in that database), whose closest NCBI match is an uncultured member of the Chlamydomonadales at 95.0% over 337 bp. **Nothing was ever submitted**, so no database correction is required. The yeast lineage as reported is therefore **71 ASVs, 647,397 reads, 40 of 45 samples** (superseded: 78 / 647,654 / 40).

---
**ANALYSIS CLOSED.** No further recomputation without instruction.

---
## Change record

| date | change | approved by | detail |
|---|---|---|---|
| 24-09-2026 | Section 9: row added, "Tremellomycetes lineage placement: 68 of 71 ASVs placed within Tremellomycetes" | Mark Mqeku, KMD PROMPT 3, P5 | value already used in Methods 2.6; recomputed from the placement files on 23-09-2026. No other value changed |
| 24-09-2026 | Section 9: row added, "UNITE release: fungal QIIME release 10.0, 4 April 2024, DOI 10.15156/BIO/2959336" | Mark Mqeku, KMD PROMPT 3, P5 | the DOI row already present is unchanged; the new row adds the release date and the evidence chain. No other value changed |

| 25-09-2026 | Section 5: Jaccard axis 1 34.4% -> **34.1%**; Christiana axis-1 coordinates and the 0.033 -> **0.032** Bray-Curtis difference added at 57,124 | Mark Mqeku, KMD PROMPT 4 (revision 4), R6a | the cited script still carried the superseded depth 57,141; rerun at 57,124 with the recorded seed; Figure 6 already showed 34.1% |
| 25-09-2026 | Section 6: effort-corrected sharing at 57,124: Christiana 179.6 (+38.6), Pretoria 210.0 vs 167.4 (-42.6) | Mark Mqeku, KMD PROMPT 4 (revision 4), R6a | same cause as above; superseded values kept in the section |
| 25-09-2026 | Section 7: Didymellaceae and Mycosphaerellaceae exact counts added | Mark Mqeku, KMD PROMPT 4 (revision 4), R4 | new computation, `phase3/r4_didymella.py` |
| 25-09-2026 | Section 9: Bloemfontein share withdrawn from reporting; detection 24 of 29 added | Mark Mqeku, KMD PROMPT 4 (revision 4), R5 | presence and absence only at Bloemfontein |
| 25-09-2026 | Section 9: rerun-tree placement **70 of 71** added; the 68-of-71 row retitled as a reference-hit count | Mark Mqeku, KMD PROMPT 4 (revision 4), R6c | **the placement count changed**: under the tree criterion the rerun (and the original tree) give 70 of 71, at clade support 44 (21 in the original tree) |
| 25-09-2026 | Section 9: amplicon lengths for the size-selection test added | Mark Mqeku, KMD PROMPT 4 (revision 4), R0c | values from `phase3/r0c_size_selection.py` |

| 25-09-2026 | Section 8: Pretoria asymptomatic *Curvibasidium* median 4.64% -> **4.63%** | Mark Mqeku, KMD PROMPT 5, D0d | computed from unrounded read counts (4.634%); 4.64 was the mean of the rounded values 3.03 and 6.24 |
| 25-09-2026 | Section 11: GenBank set is all 71 lineage ASVs, 67 uncultured Tremellomycetes and 4 uncultured fungus | Mark Mqeku, KMD PROMPT 5, D2 | two-test classification; replaces the 68-sequence set |

The analysis remains closed. The 24-09-2026 rows record values already established. The 25-09-2026 rows are the changes KMD PROMPT 4 (revision 4) authorised; every superseded value is kept beside its replacement.

*Resolved 25-09-2026 under KMD PROMPT 5 D0d (median set to 4.63%):* **Not changed, flagged 25-09-2026:** Section 8 gave the Pretoria asymptomatic *Curvibasidium* median as 4.64%, the mean of the rounded values 3.03 and 6.24. From the unrounded values (2,603 of 85,895 and 7,304 of 117,106 reads) it is 4.634%, which the manuscript reports as 4.63%. Mark to choose which to keep.
