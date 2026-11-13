#!/usr/bin/env Rscript
# Stage 5, step 1: significant probes. From every association result for the sample set, keeps the
# probes with P <= 0.05 / Meff, where Meff is the maximum across the 13 traits for that model and
# sample set (results/03_meff/<version>/meff_thresholds.csv).
# Output: results/05_signals/<version>/<all|unrelated>/ (per-file summary, significant results and
#         significant_probes_unpruned_<sample set>.csv)
# Usage: Rscript scripts/05_signals/01_sig_extraction.R <all|unrelated>
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
})

source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
args <- commandArgs(trailingOnly = TRUE)
rel_arg <- if (length(args)) args[1] else "all"
stopifnot(rel_arg %in% c("all", "unrelated"))

cat("=== QC: UNPRUNED COMBINED PHASE SIGNIFICANCE + SIGNIFICANT PROBE EXTRACTION ===\n")
cat("Working directory:", getwd(), "\n\n")

# ============================================================================
# CONFIGURATION & PATHS
# ============================================================================

relatedness_group <- if (rel_arg == "all") "related" else "unrelated"   # file-name label for each sample set

gwas_base_dir <- file.path(pa$results_root, "04_gwas", ver)
meff_dir      <- file.path(pa$results_root, "03_meff", ver)

qc_dir <- file.path(pa$results_root, "05_signals", ver, rel_arg)
dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(qc_dir, "significance_results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(qc_dir, "inputs"), showWarnings = FALSE, recursive = TRUE)

cat("GWAS base dir:", gwas_base_dir, "\n")
cat("Relatedness group:", relatedness_group, "\n")
cat("M_eff dir:", meff_dir, "\n")
cat("QC output dir:", qc_dir, "\n\n")

# ============================================================================
# LOAD M_EFF THRESHOLDS
# ============================================================================

cat("=== LOADING UNPRUNED M_EFF THRESHOLDS ===\n")

tasks <- fread(file.path(meff_dir, "meff_all_tasks.csv"))
thr   <- fread(file.path(meff_dir, "meff_thresholds.csv"))
meff_data <- merge(
  tasks[, .(CNV_Type, Phenotype, Relatedness, N_Samples, N_Probes, Meff_trait = Meff)],
  thr[, .(CNV_Type = model, Relatedness = relatedness, Meff, Alpha_0.05 = p_cnv_wide)],
  by = c("CNV_Type", "Relatedness"))
cat("Thresholds (0.05 / maximum Meff across traits):\n"); print(thr[, .(model, relatedness, Meff, p_cnv_wide)])

cat("Loaded M_eff entries:", nrow(meff_data), "\n")
cat("CNV types:", paste(unique(meff_data$CNV_Type), collapse = ", "), "\n")
cat("Relatedness labels in M_eff:", paste(unique(meff_data$Relatedness), collapse = ", "), "\n\n")

# ============================================================================
# HELPER: Parse filename
# ============================================================================

# Expected pattern:
# gwas_linear_unpruned_combined_{cnv_type}_{phenotype}_mincar{N}_{sample set}.csv

parse_unpruned_gwas_filename <- function(filename) {
  if (!grepl("^gwas_linear_unpruned_combined_", filename)) return(NULL)

  cnv_type <- NA_character_
  if (grepl("_combined_deletion_",    filename)) cnv_type <- "deletion"
  if (grepl("_combined_duplication_", filename)) cnv_type <- "duplication"
  if (grepl("_combined_mirror_",      filename)) cnv_type <- "mirror"
  if (is.na(cnv_type)) return(NULL)

  relatedness_lookup <- if (grepl("_unrelated\\.csv$", filename)) {
    "unrelated"
  } else if (grepl("_all\\.csv$", filename)) {
    "all"
  } else {
    return(NULL)
  }

  # Extract phenotype: between _combined_{cnv_type}_ and _mincar15_
  pheno_pattern <- paste0("_combined_", cnv_type, "_(.+?)_mincar[0-9]+_")
  m <- regmatches(filename, regexec(pheno_pattern, filename))[[1]]
  if (length(m) < 2) return(NULL)

  phenotype <- m[2]  # e.g. bmi_at_25_rint

  list(
    cnv_type           = cnv_type,
    phenotype          = phenotype,
    relatedness_lookup = relatedness_lookup
  )
}

# ============================================================================
# ANALYZE GWAS FILES
# ============================================================================

cat("=== ANALYZING UNPRUNED GWAS FILES ===\n\n")

gwas_rel_dir <- file.path(gwas_base_dir, rel_arg)
if (!dir.exists(gwas_rel_dir)) stop("ERROR: GWAS directory not found: ", gwas_rel_dir)

gwas_files <- list.files(gwas_rel_dir,
                          pattern = "^gwas_linear_unpruned_combined_.*\\.csv$",
                          full.names = TRUE)
cat("Found", length(gwas_files), "GWAS files in:", gwas_rel_dir, "\n\n")
if (length(gwas_files) == 0) stop("ERROR: No unpruned GWAS files found")

all_results    <- list()
all_sig_probes <- list()

for (gwas_file in gwas_files) {
  filename <- basename(gwas_file)
  cat("Processing:", filename, "\n")

  parsed <- parse_unpruned_gwas_filename(filename)
  if (is.null(parsed)) {
    cat("  Skipping: does not match expected unpruned filename pattern\n\n")
    next
  }

  cnv_type           <- parsed$cnv_type
  phenotype          <- parsed$phenotype
  relatedness_lookup <- parsed$relatedness_lookup

  cat("  Parsed:", cnv_type, "|", phenotype, "| relatedness =", relatedness_lookup, "\n")

  meff_row <- meff_data[
    CNV_Type    == cnv_type &
    Phenotype   == phenotype &
    Relatedness == relatedness_lookup
  ]

  if (nrow(meff_row) == 0) {
    cat("  WARNING: No M_eff threshold found for:", cnv_type, phenotype, relatedness_lookup, "\n\n")
    next
  }

  threshold <- as.numeric(meff_row$Alpha_0.05[1])

  gwas_dt <- fread(gwas_file)
  if (!("pval" %in% names(gwas_dt))) {
    cat("  ERROR: No pval column found\n\n")
    next
  }

  n_total <- nrow(gwas_dt)
  n_sig   <- sum(gwas_dt$pval < threshold, na.rm = TRUE)
  min_p   <- min(gwas_dt$pval, na.rm = TRUE)

  cat("  Threshold:", sprintf("%.2e", threshold),
      "| significant probes:", n_sig, "/", n_total,
      "| min p:", sprintf("%.2e", min_p), "\n")

  all_results[[length(all_results) + 1]] <- data.table(
    analysis_type      = "unpruned",
    cnv_type           = cnv_type,
    phenotype          = phenotype,
    relatedness_group  = relatedness_group,
    relatedness_lookup = relatedness_lookup,
    meff_threshold     = threshold,
    n_samples          = meff_row$N_Samples[1],
    n_probes_meff      = meff_row$N_Probes[1],
    n_total_probes     = n_total,
    n_significant      = n_sig,
    min_pval           = min_p,
    is_significant     = (n_sig > 0),
    filename           = filename,
    gwas_file          = gwas_file
  )

  if (n_sig > 0) {
    sig_probes <- gwas_dt[pval < threshold]
    sig_probes[, `:=`(
      phenotype          = phenotype,
      cnv_type           = cnv_type,
      relatedness_group  = relatedness_group,
      relatedness_lookup = relatedness_lookup,
      analysis_type      = "unpruned",
      meff_threshold     = threshold,
      gwas_filename      = filename,
      gwas_file          = gwas_file
    )]

    for (col in c("probe_id", "chr", "pos", "pval", "beta", "se")) {
      if (!col %in% names(sig_probes)) sig_probes[, (col) := NA]
    }

    all_sig_probes[[length(all_sig_probes) + 1]] <- sig_probes
  }

  cat("\n")
}

if (length(all_results) == 0) stop("ERROR: No valid results found")
results_dt  <- rbindlist(all_results, use.names = TRUE, fill = TRUE)
sig_results <- results_dt[is_significant == TRUE]

# ============================================================================
# SAVE OUTPUTS
# ============================================================================

cat("=== SAVING QC OUTPUTS ===\n")

sig_summary_file <- file.path(qc_dir, "significance_results",
                               paste0("significance_summary_unpruned_", relatedness_group, ".csv"))
fwrite(results_dt, sig_summary_file)
cat("Saved:", sig_summary_file, "\n")

if (nrow(sig_results) > 0) {
  detailed_sig_file <- file.path(qc_dir, "significance_results",
                                  paste0("detailed_significant_results_unpruned_", relatedness_group, ".csv"))
  fwrite(sig_results, detailed_sig_file)
  cat("Saved:", detailed_sig_file, "\n")
}

sig_probes_file <- file.path(qc_dir,
                              paste0("significant_probes_unpruned_", relatedness_group, ".csv"))
if (length(all_sig_probes) > 0) {
  combined_sig_probes <- rbindlist(all_sig_probes, use.names = TRUE, fill = TRUE)
  fwrite(combined_sig_probes, sig_probes_file)
  cat("Saved significant probes:", sig_probes_file, "\n")
  cat("  Total significant probes:", nrow(combined_sig_probes), "\n")
  cat("  Phenotypes included:", paste(sort(unique(combined_sig_probes$phenotype)), collapse = ", "), "\n")
} else {
  cat("! No significant probes found; saved summaries only.\n")
}

pointer_file <- file.path(qc_dir, "inputs", "qc_input_paths.txt")
writeLines(c(
  paste0("qc_dir=", qc_dir),
  paste0("relatedness_group=", relatedness_group),
  paste0("sig_probes_file=", sig_probes_file),
  paste0("gwas_rel_dir=", gwas_rel_dir)
), pointer_file)
cat("Wrote pointer file:", pointer_file, "\n")

cat("\n=== QC SIGNIFICANCE EXTRACTION COMPLETE ===\n")
cat("QC dir:", qc_dir, "\n")
cat("Next step: run Script 2 using sig_probes_file above.\n")
