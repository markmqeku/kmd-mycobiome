#!/usr/bin/env Rscript
# H-4: descriptive Hill numbers (q=0,1,2) with iNEXT confidence intervals + sample coverage.
# Core 16 only. NO significance testing: every cell is n=1-2.
suppressMessages({library(iNEXT)})
O <- "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
tab <- read.delim(file.path(O, "fungi_only/table_exp/feature-table.tsv"),
                  skip = 1, check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(O, "metadata_v2_27-07-2026.tsv"), check.names = FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]
rownames(meta) <- meta[[1]]
core <- rownames(meta)[meta$location %in% c("Christiana", "Pretoria")]
core <- intersect(core, colnames(tab))
cat("core samples:", length(core), "\n")

abun <- lapply(core, function(s) { v <- tab[[s]]; v <- v[v > 0]; sort(v, decreasing = TRUE) })
names(abun) <- core

DEPTH <- 57141   # H-1: highest depth retaining all 16
out <- iNEXT(abun, q = c(0, 1, 2), datatype = "abundance", size = c(DEPTH), nboot = 50)

cat("\n=== SAMPLE COVERAGE and observed richness (per sample) ===\n")
di <- out$DataInfo
di$site <- meta[di$Assemblage, "location"]
di$cond <- meta[di$Assemblage, "condition_std"]
di$tissue <- meta[di$Assemblage, "tissue"]
print(di[, c("Assemblage", "site", "cond", "tissue", "n", "S.obs", "SC")], row.names = FALSE)

cat("\n=== HILL NUMBERS at depth", DEPTH, "with 95% CI (per sample) ===\n")
iz <- out$iNextEst$size_based
iz <- iz[iz$m == DEPTH, c("Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC")]
iz$site <- meta[iz$Assemblage, "location"]
iz$cond <- meta[iz$Assemblage, "condition_std"]
iz$tissue <- meta[iz$Assemblage, "tissue"]
iz <- iz[order(iz$site, iz$cond, iz$tissue, iz$Order.q), ]
for (q in c(0, 1, 2)) {
  cat(sprintf("\n-- q = %d (%s) --\n", q,
      c("richness", "exp(Shannon)", "inverse Simpson")[q + 1]))
  s <- iz[iz$Order.q == q, ]
  for (i in seq_len(nrow(s))) {
    cat(sprintf("   %-6s %-11s %-13s %-14s  qD=%8.1f  [%8.1f, %8.1f]  SC=%.4f\n",
        s$Assemblage[i], s$site[i], s$cond[i], s$tissue[i],
        s$qD[i], s$qD.LCL[i], s$qD.UCL[i], s$SC[i]))
  }
}

cat("\n=== GROUP SUMMARY (site x condition): mean of per-sample estimates, n stated ===\n")
cat("NOTE: n = 1-2 per cell. These are descriptive point estimates, NOT tested.\n")
for (q in c(0, 1, 2)) {
  s <- iz[iz$Order.q == q, ]
  key <- paste(s$site, s$cond)
  cat(sprintf("\n-- q = %d --\n", q))
  for (k in unique(key)) {
    v <- s$qD[key == k]
    cat(sprintf("   %-28s n=%d  mean qD=%8.1f  range=[%.1f, %.1f]\n",
        k, length(v), mean(v), min(v), max(v)))
  }
}
cat("\nH4 INEXT DONE\n")
