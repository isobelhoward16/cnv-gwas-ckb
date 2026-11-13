#!/usr/bin/env Rscript
# Stage 10: CNV landscape numbers quoted in the Results and the BMI SD used in the Discussion, from the stage 1
# callset and stage 2 matrices, in participants passing CNV QC.
# Aggregate only. Output: results/10_paper/<version>/landscape_numbers.csv (item, value).
# Run from the repository root:  Rscript scripts/10_paper/06_landscape_numbers.R
suppressPackageStartupMessages(library(Matrix))
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version; R <- pa$results_root
od <- file.path(R, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
out <- list(); add <- function(item, value) { out[[length(out) + 1]] <<- data.table(item = item, value = as.character(value)); cat(sprintf("%-70s %s\n", item, value)) }
f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ","); fi <- function(x) formatC(x, format = "d", big.mark = ",")

qc <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
mdir <- file.path(R, "02_matrices", ver)
del <- readRDS(file.path(mdir, "combined_deletion_matrix_shared_probes_fixed.rds"))
ids <- intersect(rownames(del), qc); N <- length(ids)
add("Participants passing CNV QC (matrix rows in QC list)", fi(N))
cs <- readRDS(file.path(R, "01_callset", ver, "final", "callset.rds")); setDT(cs); cs <- cs[CKID %in% ids]
cs[, len := as.numeric(EndPos - StartPos)]   # as Figure 1 (02_figure1_callset.R: end - start)
add("CNV calls", fi(nrow(cs)))
nc <- uniqueN(cs$CKID); add("Participants with at least one CNV", fi(nc)); add("Participants with no CNV", fi(N - nc))
tab <- cs[, .N, by = CN][order(CN)]
for (i in seq_len(nrow(tab))) add(sprintf("Calls with CN = %s: n (%% of calls)", tab$CN[i]), sprintf("%s (%.1f%%)", fi(tab$N[i]), 100 * tab$N[i] / nrow(cs)))
add("Calls by CN_type", paste(cs[, .N, by = CN_type][, sprintf("%s %s (%.1f%%)", CN_type, fi(N), 100 * N / nrow(cs))], collapse = "; "))
per <- cs[, .N, by = CKID]; cnt <- c(per$N, rep(0L, N - nc))
add("CNVs per participant: mean | median (all QC-passing, zeros included)", sprintf("%.2f | %s", mean(cnt), median(cnt)))
add("Median length, deletions (kb)", f1(median(cs[CN_type == "del", len]) / 1e3))
add("Median length, duplications (kb)", f1(median(cs[CN_type == "dup", len]) / 1e3))
bt <- cs[, .(bp = sum(as.numeric(len))), by = CN_type]
add("Share of total CNV bp from duplications (%)", f1(100 * bt[CN_type == "dup", bp] / sum(bt$bp)))
bur <- cs[, .(kb = sum(as.numeric(len)) / 1e3), by = CKID]
add("Burden per carrier (kb): mean | median", sprintf("%s | %s", fi(round(mean(bur$kb))), fi(round(median(bur$kb)))))

# per-probe frequencies, 541,024 shared probes, all QC-passing participants and by array
dup <- readRDS(file.path(mdir, "combined_duplication_matrix_shared_probes_fixed.rds"))
stopifnot(identical(rownames(del), rownames(dup)))
pm <- readRDS(file.path(mdir, "combined_probe_map_shared_probes.rds")); setDT(pm)
spi <- readRDS(file.path(mdir, "combined_sample_phase_info_fixed.rds")); setDT(spi)
add("Shared probes (matrix columns)", fi(ncol(del)))
add("Probe map rows, CKB1 (phase12) full", fi(nrow(fread(pa$probe_map_phase12))))
add("Probe map rows, CKB2 (phase3) full", fi(nrow(fread(pa$probe_map_phase3))))
frq <- function(m, rows) { s <- m[rows, , drop = FALSE]; colSums(s == 1, na.rm = TRUE) / length(rows) }
grp <- list(all = ids, CKB1 = intersect(ids, spi[phase == "phase12", sample_id]), CKB2 = intersect(ids, spi[phase == "phase3", sample_id]))
add("Participants by array in this set: CKB1 | CKB2", sprintf("%s | %s", fi(length(grp$CKB1)), fi(length(grp$CKB2))))
for (g in names(grp)) for (ty in c("deletion", "duplication")) {
  f <- frq(if (ty == "deletion") del else dup, grp[[g]]); w <- which.max(f)
  add(sprintf("Per-probe %s frequency, %s: mean | median | max (probe, chr:pos)", ty, g),
      sprintf("%.5f | %.5f | %.4f (%s, chr%s:%s)", mean(f), median(f), max(f), colnames(del)[w],
              pm[probe_id == colnames(del)[w], chr], fi(pm[probe_id == colnames(del)[w], pos])))
  if (g == "all") add(sprintf("Probes with %s frequency > 1%%, all: n (%% of shared probes)", ty), sprintf("%s (%.2f%%)", fi(sum(f > 0.01)), 100 * mean(f > 0.01)))
}
# BMI SD in the BMI analysis sample (to convert the 16p11.2 BP2-BP3 duplication beta to kg/m2 in the Discussion)
ph <- readRDS(pa$phenotype_rds); setDT(ph)
k1 <- fread(pa$array1_key); k2 <- fread(pa$array2_key); setnames(k1, c("ccvid", "xIID_cnv")); setnames(k2, c("ccvid", "xIID_cnv"))
pa_all <- rbindlist(list(merge(ph[gwas_array_type == "Axiom_CKB1"], k1, by = "ccvid"), merge(ph[gwas_array_type == "Axiom_CKB2"], k2, by = "ccvid")), fill = TRUE)
b <- pa_all[xIID_cnv %in% ids & !is.na(bmi_calc), bmi_calc]
add("BMI (kg/m2) in QC-passing participants with BMI: n | mean | SD", sprintf("%s | %.4f | %.4f", fi(length(b)), mean(b), sd(b)))
fwrite(rbindlist(out), file.path(od, "landscape_numbers.csv")); cat("Saved:", file.path(od, "landscape_numbers.csv"), "\n")
