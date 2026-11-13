#!/usr/bin/env Rscript
# Stage 10: the numbers quoted in the manuscript text, abstract and legends, computed from the pipeline outputs,
# with the file each comes from. A block whose inputs are missing prints NOT AVAILABLE and is skipped.
# Aggregate output only. Output: results/10_paper/<version>/manuscript_numbers.csv (item, value, source)
source("R/config.R"); suppressPackageStartupMessages(library(Matrix))   # needed for rownames() of the sparse matrix
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
mc <- cfg$pipeline$analysis$min_carriers
options(width = 200)
R <- pa$results_root
od <- file.path(R, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
out <- list()
add <- function(item, value, source) out[[length(out) + 1]] <<- data.table(item = item, value = as.character(value), source = source)
fmt <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
sci1 <- function(p) { if (is.na(p)) return("NA"); ex <- floor(log10(p)); m <- p / 10^ex; if (round(m, 2) >= 10) { m <- m / 10; ex <- ex + 1 }; sprintf("%.2f x 10^%d", m, ex) }
sci <- function(p) vapply(p, sci1, character(1))   # vectorised
block <- function(name, expr) {
  r <- tryCatch({ force(expr); TRUE }, error = function(e) { cat("NOT AVAILABLE:", name, "|", conditionMessage(e), "\n"); FALSE })
  if (r) cat("done:", name, "\n")
}
rel <- function(...) sub(paste0("^", R, "/?"), "", file.path(...))

block("sample and calls", {
  f <- file.path(R, "08_descriptive", ver, "data/derived/tables/ckb_denominator_full_qc.txt")
  add("Participants passing CNV QC", fmt(as.integer(readLines(f))), rel(f))
  qc <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
  ids <- intersect(rownames(readRDS(file.path(R, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds"))), qc)
  if (!length(ids)) stop("no participant IDs from the matrix")
  cs <- readRDS(file.path(R, "01_callset", ver, "final", "callset.rds")); setDT(cs); cs <- cs[CKID %in% ids]
  add("CNV calls in participants passing CNV QC", fmt(nrow(cs)), "01_callset/<ver>/final/callset.rds")
  add("Deletions", sprintf("%s (%.1f%%)", fmt(sum(cs$CN_type == "del")), 100 * mean(cs$CN_type == "del")), "callset.rds")
  add("Duplications", sprintf("%s (%.1f%%)", fmt(sum(cs$CN_type == "dup")), 100 * mean(cs$CN_type == "dup")), "callset.rds")
  add("Mean CNVs per participant (QC-passing)", sprintf("%.2f", nrow(cs) / length(ids)), "callset.rds / QC-passing N")
})
block("Table S2 (analysis population)", {
  f <- file.path(R, "08_descriptive", ver, "tables", "table_s2_population.csv"); s2 <- fread(f)
  for (k in c("Age (years)", "Sex (% female)", "Axiom_CKB1", "Axiom_CKB2")) add(paste("Table S2:", k), paste(s2[Characteristic == k, Statistic], "| N", s2[Characteristic == k, N]), rel(f))
})
block("thresholds and probes tested", {
  f <- file.path(R, "03_meff", ver, "meff_thresholds.csv"); th <- fread(f)
  for (i in seq_len(nrow(th))) add(sprintf("Meff, %s, %s", th$model[i], th$relatedness[i]),
    sprintf("%s | CNV-wide P %s | experiment-wide P %s", fmt(th$Meff[i]), sci(th$p_cnv_wide[i]), sci(th$p_experiment_wide[i])), rel(f))
  f2 <- file.path(R, "03_meff", ver, "meff_all_tasks.csv"); mt <- fread(f2)
  for (m in unique(mt$CNV_Type)) for (r in unique(mt$Relatedness)) { y <- mt[CNV_Type == m & Relatedness == r]
    add(sprintf("Probes tested, %s, %s (range over traits)", m, r), sprintf("%s-%s", fmt(min(y$N_Probes)), fmt(max(y$N_Probes))), rel(f2)) }
  add("Trait N range (all participants)", sprintf("%s-%s", fmt(min(mt[Relatedness == "all", N_Samples])), fmt(max(mt[Relatedness == "all", N_Samples]))), rel(f2))
})
block("signals", {
  f <- file.path(R, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"); x <- fread(f)
  add("Model-specific associations (all participants)", nrow(x), rel(f))
  add("  of which experiment-wide", sum(x$pass_gw_experiment), rel(f))
  add("  by model", paste(x[, .N, by = cnv_type][, sprintf("%s %d", cnv_type, N)], collapse = ", "), rel(f))
  add("Traits with at least one association", uniqueN(x$phenotype), rel(f))
  fu <- file.path(R, "05_signals", ver, "unrelated", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_unrelated.csv")
  if (file.exists(fu)) add("Model-specific associations (unrelated GWAS)", nrow(fread(fu)), rel(fu))
})
block("Table 1", {
  f <- file.path(R, "10_paper", ver, "table1_numeric.csv"); t1 <- fread(f)
  add("Table 1: locus-trait associations", nrow(t1), rel(f))
  add("Table 1: loci", uniqueN(sub("[ab]$", "", t1$Locus[t1$Locus != ""])), rel(f))
  add("Table 1: loci (lettered sub-loci counted separately)", uniqueN(t1$Locus), rel(f))
  for (i in seq_len(nrow(t1))) add(sprintf("Table 1 row: %s %s", t1$Locus[i], t1$Trait[i]),
    sprintf("%s model | beta %.2f (SE %.2f) | P %s | del/dup %d/%d | N %s | exp-wide %s | unrelated %s", t1$model_shown[i], t1$beta[i], t1$se[i],
            sci(t1$p[i]), t1$n_del[i], t1$n_dup[i], fmt(t1$N[i]), t1$exp_wide[i], t1$unrelated[i]), rel(f))
})
block("unrelated sensitivity (Table S9)", {
  f <- file.path(R, "06_followup", ver, "unrelated_leads", "table_s9_unrelated_sensitivity.csv"); s <- fread(f)
  add("Unrelated: leads with an estimate", sprintf("%d of %d", sum(!is.na(s$beta_unrel)), nrow(s)), rel(f))
  add("Unrelated: direction concordant", sprintf("%d of %d", sum(s$concordant, na.rm = TRUE), nrow(s)), rel(f))
  add("Unrelated: beta ratio median (range)", sprintf("%.2f (%.2f-%.2f)", median(s$ratio, na.rm = TRUE), min(s$ratio, na.rm = TRUE), max(s$ratio, na.rm = TRUE)), rel(f))
  add("Unrelated: CNV-wide significant", sprintf("%d of %d", sum(s$pass_unrel), nrow(s)), rel(f))
  add("Unrelated: N range at the leads", sprintf("%s-%s", fmt(min(s$N_unrel, na.rm = TRUE)), fmt(max(s$N_unrel, na.rm = TRUE))), rel(f))
})
block("fasting-adjusted random glucose", {
  th <- fread(file.path(R, "03_meff", ver, "meff_thresholds.csv"))[relatedness == "all"]
  d <- file.path(R, "04_gwas", ver, "all")
  for (m in c("deletion", "duplication", "mirror")) {
    f <- file.path(d, "sensitivity", sprintf("gwas_linear_unpruned_combined_%s_random_glucose_rint_mincar%s_all_fastadj.csv", m, mc))
    if (!file.exists(f)) { add(paste("Fasting-adjusted RPG,", m), "no file", rel(f)); next }
    a <- fread(file.path(d, sprintf("gwas_linear_unpruned_combined_%s_random_glucose_rint_mincar%s_all.csv", m, mc)), select = c("probe_id", "beta", "pval"))
    b <- fread(f, select = c("probe_id", "beta", "pval", "n_analyzed")); x <- merge(a, b, by = "probe_id", suffixes = c("", "_adj")); t <- th[model == m, p_cnv_wide]
    sig <- x[pval < t]
    add(paste("Fasting-adjusted RPG,", m), if (!nrow(sig)) "no main signal" else
      paste(sig[, sprintf("%s: P %s -> %s, beta %.3f -> %.3f, N %s", probe_id, sci(pval), sci(pval_adj), beta, beta_adj, fmt(n_analyzed))], collapse = "; "), rel(f))
  }
})
block("Table S4 (genomic disorders)", {
  f <- file.path(R, "08_descriptive", ver, "results/tables/main/Table_pathogenic_ckb_vs_ukb.tsv"); s <- fread(f)
  add("Table S4: loci", nrow(s), rel(f)); add("Table S4: loci with no CKB carrier", sum(s$ckb_carriers == 0), rel(f))
  add("Table S4: total CKB carriers", fmt(sum(s$ckb_carriers)), rel(f))
})
res <- rbindlist(out)
cat("\n== Manuscript numbers (", ver, ")\n"); print(res, nrows = 500)
fwrite(res, file.path(od, "manuscript_numbers.csv")); cat("\nSaved:", file.path(od, "manuscript_numbers.csv"), "\n")
