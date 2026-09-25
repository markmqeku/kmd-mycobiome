#!/usr/bin/env Rscript
# Item 10b: export iNEXT per-sample Hill values at depth 57,141 as a machine-readable TSV,
# so Figure 5 has a proper source file instead of a parsed log.
suppressMessages({library(iNEXT)})
O <- "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
G <- file.path(O, "gap_closure"); dir.create(G, showWarnings = FALSE)
tab <- read.delim(file.path(O, "fungi_only/table_exp/feature-table.tsv"),
                  skip = 1, check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(O, "metadata_v2_27-07-2026.tsv"), check.names = FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]; rownames(meta) <- meta[[1]]
core <- intersect(rownames(meta)[meta$location %in% c("Christiana", "Pretoria")], colnames(tab))
abun <- lapply(core, function(s) { v <- tab[[s]]; v <- v[v > 0]; sort(v, decreasing = TRUE) })
names(abun) <- core
DEPTH <- 57141
out <- iNEXT(abun, q = c(0, 1, 2), datatype = "abundance", size = c(DEPTH), nboot = 50)
iz <- out$iNextEst$size_based
iz <- iz[iz$m == DEPTH, c("Assemblage", "Order.q", "qD", "qD.LCL", "qD.UCL", "SC")]
iz$site      <- meta[iz$Assemblage, "location"]
iz$condition <- meta[iz$Assemblage, "condition_std"]
iz$tissue    <- meta[iz$Assemblage, "tissue"]
iz$depth     <- DEPTH
names(iz)[names(iz) == "Assemblage"] <- "sample_id"
names(iz)[names(iz) == "Order.q"]    <- "hill_q"
iz <- iz[order(iz$site, iz$condition, iz$tissue, iz$hill_q),
         c("sample_id","site","condition","tissue","depth","hill_q","qD","qD.LCL","qD.UCL","SC")]
f <- file.path(G, "hill_numbers_per_sample.tsv")
write.table(iz, f, sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote", f, "with", nrow(iz), "rows\n")
for (q in c(0,1,2)) {
  s <- iz[iz$hill_q == q, ]
  cat(sprintf("\n-- q=%d group means --\n", q))
  for (k in unique(paste(s$site, s$condition))) {
    v <- s$qD[paste(s$site, s$condition) == k]
    cat(sprintf("   %-28s n=%d mean=%.1f\n", k, length(v), mean(v)))
  }
}
