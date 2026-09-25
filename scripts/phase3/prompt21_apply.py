#!/usr/bin/env python
"""PROMPT 21, corrections C-1 to C-15. Every substitution asserts its target is present exactly
once, so a silent miss is impossible. Before and after are printed for each."""
import io, sys
K = "<KMD_ROOT>"
log = []

def apply(fname, subs):
    p = f"{K}/{fname}"
    t = io.open(p, encoding="utf-8").read()
    for tag, a, b in subs:
        n = t.count(a)
        if n != 1:
            print(f"  !! {tag}: target appears {n} times in {fname}, SKIPPED"); continue
        t = t.replace(a, b)
        log.append((tag, fname, a, b))
    io.open(p, "w", encoding="utf-8").write(t)

# ---------------- C-1 ----------------
apply("DRAFT_RESULTS.md", [
 ("C-1a",
  "\\*Curvibasidium\\* was the only genus recovered at reliable thresholds whose abundance differed in the same direction between conditions at both sites (Figure 7). Averaged across each group, it rose from 0.24 to 5.82 percent at Christiana and from 4.88 to 11.28 percent at Pretoria.",
  "\\*Curvibasidium\\* relative abundance differed in the same direction between conditions at both sites (Figure 7). At Christiana the median was 0.24 percent in asymptomatic tissue (n = 3) and 3.76 percent in symptomatic tissue (n = 6). At Pretoria it was 4.63 percent (n = 2) and 10.73 percent (n = 5). It was the only genus recovered at reliable thresholds of which this was true."),
 ("C-1b", "Individual sample values show how strong this pattern is.",
  "Values for individual samples are given below."),
 # the third paragraph carries 4.64, which is a rounding error on the same median
 ("C-1c", "giving a median of 4.64 percent", "giving a median of 4.63 percent"),
 # ---------------- C-2 ----------------
 ("C-2", "Christiana and Pretoria contained no host sequence and retained all samples at every rarefaction depth tested.",
  "Christiana and Pretoria contained no host sequence and retained all 16 samples at the rarefaction depth used."),
 # ---------------- C-3 ----------------
 ("C-3a", "\\### 3.2 Host content controlled apparent fungal richness at Bloemfontein",
  "\\### 3.2 Apparent fungal richness varied inversely with host content at Bloemfontein"),
 ("C-3b", "the number of fungal ASVs recovered from a sample depended on how much host DNA that sample contained, not on its tissue type or the age of the tree",
  "the number of fungal ASVs recovered from a sample varied inversely with its host DNA content, and samples did not group by tissue type or age class"),
 # ---------------- C-4 ----------------
 ("C-4", "That sequence was deposited on 1 November 2019, before the release of the UNITE database used here. The lineage is therefore genuinely absent from the reference database rather than merely awaiting an update.",
  "That sequence was deposited on 1 November 2019 and was still absent from the UNITE release downloaded for this study in July 2026, nearly seven years later."),
 # ---------------- C-10 ----------------
 ("C-10", "This lineage accounted for 41.5 percent of reads in asymptomatic tissue at Christiana and 39.6 percent in symptomatic tissue, 14.6 and 21.5 percent respectively at Pretoria, and 15.2 percent at Bloemfontein.",
  "This lineage accounted for 41.5 percent of reads in asymptomatic tissue at Christiana (n = 3) and 39.6 percent in symptomatic tissue (n = 6), 14.6 percent (n = 2) and 21.5 percent (n = 5) respectively at Pretoria, and 15.2 percent at Bloemfontein (n = 29)."),
 # ---------------- C-11 ----------------
 ("C-11", "Removing the host ASVs left 1,173 ASVs. A second screen of that whole set against the NCBI nucleotide database, described in Section 2.4 and reported in Section 3.11, identified a further 31 ASVs as non-fungal. These were removed. The corrected fungal dataset used throughout the remainder of this paper comprises \\*\\*1,142 ASVs and 2,517,177 reads\\*\\*.",
  "Removing the 149 host ASVs left 1,173 ASVs and 2,523,611 reads. One sample, P2, contained no fungal reads at all, despite producing 57,959 reads in total. A second screen of that whole set against the NCBI nucleotide database, described in Section 2.4 and reported in Section 3.11, identified a further 31 ASVs as non-fungal, and these were removed. The corrected fungal dataset used throughout the remainder of this paper therefore comprises \\*\\*1,142 ASVs and 2,517,177 reads\\*\\*, and all values reported below are computed on it."),
 ("C-11b", "After the 149 host ASVs were removed, 1,173 fungal ASVs and 2,523,611 reads remained. One sample, P2, contained no fungal reads at all, despite producing 57,959 reads in total.\n", ""),
 # ---------------- C-12 ----------------
 ("C-12a", "A group of 71 ASVs, comprising 647,397 reads and found in 40 of the 45 samples, could not be assigned a name at reliable identity thresholds.",
  "A group of 71 ASVs, comprising 647,397 reads and found in 40 of the 45 samples, could not be assigned a name at reliable identity thresholds. Of these 71, 26 occur in the 16 core samples; the remaining 45 were recovered only at Bloemfontein."),
 ("C-12b", "Of its 25 ASVs present in the core samples, 19 never reached 1 percent in any sample while 6 did.",
  "Of its 25 ASVs present in the core samples after rarefaction to 57,124 reads, 19 never reached 1 percent in any sample while 6 did."),
])

# ---------------- Discussion ----------------
apply("DRAFT_DISCUSSION.md", [
 ("C-5", """That this signal survives reanalysis while two genera reported from the
same earlier sequence processing do not is a consequence of their differing evidential
basis: the yeast observation rests on isolation and morphology, independent of the sequence
pipeline, whereas the two genera rested on the pipeline itself at a permissive identity
threshold and did not survive reclassification.""",
  """The independent evidence is Swanepoel et al.'s culture-based recovery of yeast-like
morphospecies, which does not depend on any sequence pipeline. An earlier unpublished analysis
of the present material also reported *Filobasidiella*, though that observation derives from the
same superseded processing that produced *Didymella* and *Mycosphaerella* and is therefore
consistent with, rather than independent of, the present finding."""),
 ("C-6", """and sequenced together, so the separation between them is not attributable to batch. The
separation therefore reflects a difference between the two populations rather than a
technical artefact, and suggests that any search""",
  """and sequenced together. Sequencing batch can therefore be excluded as an explanation.
Collection dates and the extent of pruning also differed between the two sites, so the
separation cannot be attributed to population differences alone. It suggests that any search"""),
 ("C-7", "identical extraction,\nprimers, operator and thermocycler programme",
  "identical extraction,\nprimers and thermocycler programme"),
 ("C-8", "once sampling effort and reference\nquality are accounted for", "once sampling effort is accounted for"),
])

# ---------------- Introduction ----------------
apply("DRAFT_INTRODUCTION.md", [
 ("C-9a", "comparing symptomatic tissue with morphologically normal tissue from the same diseased\ntrees",
  "comparing symptomatic tissue with morphologically normal tissue collected at the same\ndiseased sites"),
 ("C-9b", """a survey for members of the
latent pathogen family Botryosphaeriaceae recovered only a single species from the host
(Jami et al., 2013)""",
  """a focused survey at one location
recovered only a single species of Botryosphaeriaceae from the host, although other hosts in
the same area carried more (Jami et al., 2014)"""),
 ("C-9c", """to characterise an abundant fungal lineage that current reference
databases place incorrectly""",
  """to characterise an abundant fungal lineage that could not be
resolved using current reference data"""),
])

print(f"applied {len(log)} substitutions\n")
for tag, f, a, b in log:
    print(f"--- {tag}  ({f}) ---")
    print(f"  BEFORE: {' '.join(a.split())[:150]}")
    print(f"  AFTER : {' '.join(b.split())[:150] if b.strip() else '(deleted)'}")
