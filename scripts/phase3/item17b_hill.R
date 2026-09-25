#!/usr/bin/env Rscript
# Item 17b: recompute per-sample Hill numbers on the CORRECTED table at the re-locked
# depth 57,124 (DEC-1). Same iNEXT call as item10b_hill_export.R so the two are comparable;
# only the input table and the depth differ. Old and new are written side by side.
suppressMessages({library(iNEXT)})
O <- "<KMD_ROOT>/__reanalysis_2026-06/phase3/outputs_27-07-2026"
C <- file.path(O, "clean"); dir.create(C, showWarnings = FALSE)
tab <- read.delim(file.path(C, "feature-table-clean.tsv"), check.names = FALSE, row.names = 1)
meta <- read.delim(file.path(O, "metadata_v2_27-07-2026.tsv"), check.names = FALSE)
meta <- meta[meta[[1]] != "#q2:types", ]; rownames(meta) <- meta[[1]]
core <- intersect(rownames(meta)[meta$location %in% c("Christiana", "Pretoria")], colnames(tab))
stopifnot(length(core) == 16)
abun <- lapply(core, function(s) { v <- tab[[s]]; v <- v[v > 0]; sort(v, decreasing = TRUE) })
names(abun) <- core
DEPTH <- 57124
cat("smallest core sample:", min(sapply(abun, sum)), "| depth:", DEPTH, "\n")
stopifnot(min(sapply(abun, sum)) >= DEPTH)
set.seed(1)
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
f <- file.path(C, "hill_numbers_per_sample.tsv")
write.table(iz, f, sep = "\t", row.names = FALSE, quote = FALSE)
cat("wrote", f, "with", nrow(iz), "rows\n")

# ---- side by side against the superseded values ----
old <- read.delim(file.path(O, "gap_closure/hill_numbers_per_sample.tsv"), check.names = FALSE)
cmp <- merge(old[, c("sample_id","hill_q","qD")], iz[, c("sample_id","hill_q","qD")],
             by = c("sample_id","hill_q"), suffixes = c("_old", "_new"))
cmp$delta <- cmp$qD_new - cmp$qD_old
cmp <- cmp[order(cmp$hill_q, cmp$sample_id), ]
write.table(cmp, file.path(C, "hill_old_vs_new.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
cat("\nper-sample qD, depth 57,141 -> 57,124\n")
cat(sprintf("%-8s %3s %12s %12s %10s\n", "sample", "q", "old", "new", "delta"))
for (i in seq_len(nrow(cmp))) cat(sprintf("%-8s %3d %12.4f %12.4f %+10.4f\n",
    cmp$sample_id[i], cmp$hill_q[i], cmp$qD_old[i], cmp$qD_new[i], cmp$delta[i]))
cat("\ngroup means, old -> new\n")
for (q in c(0,1,2)) {
  s <- iz[iz$hill_q == q, ]; so <- old[old$hill_q == q, ]
  for (k in unique(paste(s$site, s$condition))) {
    v  <- s$qD[paste(s$site, s$condition) == k]
    vo <- so$qD[paste(so$site, so$condition) == k]
    cat(sprintf("   q=%d %-26s n=%d  %8.2f -> %8.2f  (%+.2f)\n",
                q, k, length(v), mean(vo), mean(v), mean(v)-mean(vo)))
  }
}
cat(sprintf("\nlargest single-sample change: %+.4f\n", cmp$delta[which.max(abs(cmp$delta))]))
