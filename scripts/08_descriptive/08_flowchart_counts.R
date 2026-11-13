#!/usr/bin/env Rscript
# Stage 8, step 8: the numbers in Figure S1 (study flowchart) and in the Methods sample and call QC sentences.
# Samples: PennCNV input (100,706) -> CNV sample QC exclusions (LRR SD, BAF drift, WF, more than 100 CNVs; one
#   sample can fail several) -> passed CNV QC (checked against the analysis set) -> phenotype available.
# Calls: raw calls in samples passing CNV QC -> call filters (per reason, and in filter order) -> final calls.
# Aggregate counts only. Output: results/08_descriptive/<version>/tables/flowchart_counts.txt (also printed)
# Run from the repository root:  Rscript scripts/08_descriptive/08_flowchart_counts.R
source("R/config.R"); source("R/callset_utils.R"); suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
q <- cfg$pipeline$sample_qc; fp <- cfg$pipeline$cnv_filter
fmt <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")
cd <- function(...) file.path(pa$results_root, "01_callset", ver, ...)
od <- file.path(pa$results_root, "08_descriptive", ver, "tables"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
out_f <- file.path(od, "flowchart_counts.txt"); sink(out_f, split = TRUE)
cat("Callset version:", ver, "\n")

# ---- samples ----
qc <- fread(pa$master_qc_summary, sep = "\t")
k1 <- fread(pa$array1_key); k2 <- fread(pa$array2_key); setnames(k1, c("ccvid", "id")); setnames(k2, c("ccvid", "id"))
qc[, array := fifelse(CKID %in% k1$id, "CKB1", fifelse(CKID %in% k2$id, "CKB2", NA_character_))]
qg <- fread(pa$qc_gwas_file)
qc[qg, NumCNV_raw := i.NumCNV, on = c(CKID = "GID")]
high <- readLines(cd("final", "samples", "high_cnv_excluded_ckids.txt"))
mat_ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))
final <- intersect(mat_ids, qg[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID])

by_arr <- function(v) sprintf("%s (CKB1 %s | CKB2 %s)", fmt(sum(v)), fmt(sum(v & qc$array == "CKB1", na.rm = TRUE)), fmt(sum(v & qc$array == "CKB2", na.rm = TRUE)))
cat("== Samples\n")
cat("Samples in the PennCNV QC summary (genotyped, PennCNV input):", by_arr(rep(TRUE, nrow(qc))), "\n")
cat("   with no array assignment from the keys:", sum(is.na(qc$array)), "\n")
r <- rep(TRUE, nrow(qc))
cat("Entering CNV sample QC:", by_arr(r), "\n")
f_lrr <- r & ((qc$LRR_SD > q$max_lrr_sd) %in% TRUE); f_baf <- r & ((abs(qc$BAF_drift) > q$max_baf_drift) %in% TRUE)
f_wf  <- r & ((abs(qc$WF) > q$max_wf) %in% TRUE);    f_na  <- r & (is.na(qc$LRR_SD) | is.na(qc$BAF_drift) | is.na(qc$WF))
f_raw <- r & ((qc$NumCNV_raw > 100) %in% TRUE);       f_hc  <- r & qc$CKID %in% high
f_num <- f_raw | f_hc
f_any <- f_lrr | f_baf | f_wf | f_na | f_num
cat("CNV sample QC, excluded (any criterion):", by_arr(f_any), "\n")
# each criterion by array (Figure S2 legend); one sample can fail several
cat("   LRR SD > 0.35:", by_arr(f_lrr), "\n")
cat("   BAF drift > 0.01:", by_arr(f_baf), "\n")
cat("   |WF| > 0.05:", by_arr(f_wf), "\n")
cat("   missing a metric:", by_arr(f_na), "\n")
cat("   More than 100 CNVs:", by_arr(f_num), "\n")
cat("      (raw PennCNV count > 100:", fmt(sum(f_raw)), "; > 100 after call filters:", fmt(sum(f_hc)),
    "; both:", fmt(sum(f_raw & f_hc)), ")\n")
p <- r & !f_any
cat("Passed CNV QC (computed here):", by_arr(p), "\n")
cat("Passed CNV QC (analysis set = matrix rows passing the GWAS QC file):", fmt(length(final)),
    "| computed-but-not-in-set:", sum(p & !qc$CKID %in% final), "| in-set-but-not-computed:", length(setdiff(final, qc$CKID[p])), "\n")
ph <- readRDS(pa$phenotype_rds); setDT(ph)
kk <- rbind(k1[, .(ccvid, id, gwas_array_type = "Axiom_CKB1")], k2[, .(ccvid, id, gwas_array_type = "Axiom_CKB2")])
ph <- merge(ph, kk, by = c("ccvid", "gwas_array_type"))[id %in% final]
tr <- grep("_rint$", names(ph), value = TRUE); tr <- intersect(tr, cfg$pipeline$analysis$phenotypes)
has <- rowSums(!is.na(as.matrix(ph[, ..tr]))) > 0
cat("With phenotype data for at least one of the", length(tr), "traits:", fmt(sum(has)), "\n")
nt <- sapply(tr, function(v) sum(!is.na(ph[[v]])))
cat("N per trait: range", fmt(min(nt)), "-", fmt(max(nt)), "\n")

# ---- calls ----
cat("\n== Calls in the", fmt(length(final)), "participants passing CNV QC\n")
x <- readRDS(cd("02_sample_qc", "calls_qc.rds")); setDT(x); x <- x[CKID %in% final]
cat("Raw calls (after merging adjacent fragments):", fmt(nrow(x)), "\n")
x[, `:=`(f_conf = !((Conf >= fp$min_conf) %in% TRUE), f_snp = !((NumSNP >= fp$min_snps) %in% TRUE),
         f_min = !((Length >= fp$min_length) %in% TRUE), f_max = !((Length <= fp$max_length) %in% TRUE))]
regc <- condense_regions(read_regions(file.path(cfg$root, cfg$pipeline$exclusion_regions_file)))
ov <- max_region_fraction(x$Chromosome, x$StartPos, x$EndPos, regc)
x[, f_reg := ov$frac >= fp$exclusion_min_fraction]
x[, f_any := f_conf | f_snp | f_min | f_max | f_reg]
cat("Calls excluded by the call filters:", fmt(sum(x$f_any)), "\n")
cat("   (one call can fail several) Conf <", fp$min_conf, ":", fmt(sum(x$f_conf)), "| < 10 probes:", fmt(sum(x$f_snp)),
    "| < 10 kb:", fmt(sum(x$f_min)), "| > 10 Mb:", fmt(sum(x$f_max)), "| >= 50% in an exclusion region:", fmt(sum(x$f_reg)), "\n")
cat("   in order (each call counted once): low confidence", fmt(sum(x$f_conf)),
    "; then < 10 probes or < 10 kb", fmt(sum(!x$f_conf & (x$f_snp | x$f_min))),
    "; then > 10 Mb", fmt(sum(!x$f_conf & !x$f_snp & !x$f_min & x$f_max)),
    "; then exclusion region", fmt(sum(!x$f_conf & !x$f_snp & !x$f_min & !x$f_max & x$f_reg)), "\n")
cs <- readRDS(cd("final", "callset.rds")); setDT(cs); cs <- cs[CKID %in% final]
cat("Calls kept by the filters:", fmt(sum(!x$f_any)), "| final callset in these participants:", fmt(nrow(cs)),
    "| difference (exact duplicate calls removed at freezing):", fmt(sum(!x$f_any) - nrow(cs)), "\n")
cat("Final calls: deletions", fmt(sum(cs$CN_type == "del")), "| duplications", fmt(sum(cs$CN_type == "dup")), "\n")
cat("Exclusion regions file:", cfg$pipeline$exclusion_regions_file, "| condensed regions:", nrow(regc), "\n")
sink()
cat("\nSaved:", out_f, "\n")
