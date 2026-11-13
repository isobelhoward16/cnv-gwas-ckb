#!/usr/bin/env Rscript
# Stage 2 check: rebuilds every expected matrix entry from the call set with a different method
# (data.table::foverlaps) and compares all three matrices entry by entry. Prints counts only.
#
# Expected coding:
#   deletion matrix   : 1 where a deletion covers the probe, NA where a duplication covers it, else 0
#   duplication matrix: 1 where a duplication covers the probe, NA where a deletion covers it, else 0
#   mirror matrix     : -1 where a deletion covers the probe, +1 where only a duplication does, else 0
source("R/config.R")
suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
md <- file.path(pa$results_root, "02_matrices", ver)
fin <- file.path(pa$results_root, "01_callset", ver, "final")
ok <- TRUE
chk <- function(pass, what, detail = "") {
  cat(sprintf("  [%s] %s %s\n", if (isTRUE(pass)) " OK " else "FAIL", what, detail))
  if (!isTRUE(pass)) ok <<- FALSE
}
fmt <- function(x) format(x, big.mark = ",")

cat("A. Samples (rows)\n")
cs <- readRDS(file.path(fin, "callset.rds"))
carriers <- readLines(file.path(fin, "samples", "carriers_ckids.txt"))
zero <- readLines(file.path(fin, "samples", "zero_cnv_ckids.txt"))
a1 <- fread(pa$array1_key)$IID; a2 <- fread(pa$array2_key)$IID
X <- list(deletion = readRDS(file.path(md, "combined_deletion_matrix_shared_probes_fixed.rds")),
          duplication = readRDS(file.path(md, "combined_duplication_matrix_shared_probes_fixed.rds")),
          mirror = readRDS(file.path(md, "combined_mirror_matrix_shared_probes_fixed.rds")))
rn <- rownames(X$mirror)
expected_rows <- intersect(c(carriers, zero), c(a1, a2))
chk(identical(rownames(X$deletion), rn) && identical(rownames(X$duplication), rn), "same rows in all three matrices")
chk(!anyDuplicated(rn), "no duplicated sample", fmt(length(rn)))
chk(setequal(rn, expected_rows), "rows = carriers + zero-CNV samples on either array",
    sprintf("(matrix %s | expected %s | carriers %s + zero %s; not on either array key: %s)",
            fmt(length(rn)), fmt(length(expected_rows)), fmt(length(carriers)), fmt(length(zero)),
            fmt(length(setdiff(c(carriers, zero), c(a1, a2))))))
ph <- readRDS(file.path(md, "combined_sample_phase_info_fixed.rds"))
chk(identical(ph$sample_id, rn), "phase info in matrix row order")
chk(all(ph[phase == "phase12"]$sample_id %in% a1) && all(ph[phase == "phase3"]$sample_id %in% a2), "phase matches array key")

cat("\nB. Probes (columns)\n")
m12 <- fread(pa$probe_map_phase12, col.names = c("chr", "probe_id", "genetic_dist", "pos"))
m3  <- fread(pa$probe_map_phase3,  col.names = c("chr", "probe_id", "genetic_dist", "pos"))
shared <- intersect(m12$probe_id, m3$probe_id)
pm <- readRDS(file.path(md, "combined_probe_map_shared_probes.rds"))
chk(setequal(colnames(X$mirror), shared) && ncol(X$mirror) == length(shared), "columns = probes on both arrays", fmt(length(shared)))
chk(identical(colnames(X$mirror), pm$probe_id), "column order = probe map order")
chk(identical(colnames(X$deletion), colnames(X$mirror)) && identical(colnames(X$duplication), colnames(X$mirror)), "same columns in all three matrices")

cat("\nC. Every entry rebuilt independently from the callset (foverlaps)\n")
cs[, phase := fifelse(CKID %in% a1, "phase12", fifelse(CKID %in% a2, "phase3", NA_character_))]
chk(!anyNA(cs$phase), "every call's sample is on an array key")
pairs_for <- function(map, ph_name) {
  p <- map[, .(chr = as.character(chr), probe_id, start = pos, end = pos)]
  setkey(p, chr, start, end)
  q <- cs[phase == ph_name, .(CKID, type = CN_type, chr = sub("^chr", "", Chromosome), start = StartPos, end = EndPos)]
  o <- foverlaps(q, p, by.x = c("chr", "start", "end"), type = "any", nomatch = NULL)
  unique(o[probe_id %in% shared, .(CKID, probe_id, type)])
}
pr <- rbind(pairs_for(m12, "phase12"), pairs_for(m3, "phase3"))
del <- unique(pr[type == "del", .(CKID, probe_id)]); dup <- unique(pr[type == "dup", .(CKID, probe_id)])
both <- fintersect(del, dup)
cat(sprintf("  sample-probe pairs: deletion %s | duplication %s | covered by both types %s\n", fmt(nrow(del)), fmt(nrow(dup)), fmt(nrow(both))))
expected <- list(
  deletion    = rbind(fsetdiff(del, dup)[, v := 1], dup[, .(CKID, probe_id, v = 99)]),
  duplication = rbind(fsetdiff(dup, del)[, v := 1], del[, .(CKID, probe_id, v = 99)]),
  mirror      = rbind(del[, .(CKID, probe_id, v = -1)], fsetdiff(dup, del)[, v := 1]))
for (m in names(X)) {
  s <- as.data.table(summary(X[[m]]))
  got <- data.table(CKID = rn[s$i], probe_id = colnames(X[[m]])[s$j], v = fifelse(is.na(s$x), 99, s$x))[v != 0]
  e <- expected[[m]]
  miss <- fsetdiff(e, got); extra <- fsetdiff(got, e)
  chk(nrow(miss) == 0 && nrow(extra) == 0, sprintf("%-11s all entries match", m),
      sprintf("(entries %s | expected %s | missing %s | unexpected %s)", fmt(nrow(got)), fmt(nrow(e)), fmt(nrow(miss)), fmt(nrow(extra))))
}
chk(all(!zero %in% rn[unique(as.data.table(summary(X$mirror))$i)]), "zero-CNV samples have no non-zero entry")
rm(pr, del, dup, both, expected); invisible(gc())

cat(if (ok) "\nAll checks passed.\n" else "\nSome checks FAILED.\n")
