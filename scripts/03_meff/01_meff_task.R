#!/usr/bin/env Rscript
# Stage 3: effective number of independent tests (simpleM) for one model x phenotype x sample set.
# Uses the probes tested in the association analysis: samples passing QC with the phenotype, and
# probes with at least min_carriers carriers (any non-zero value, so deletions count in the mirror
# matrix). Missing values are imputed with the column mode before the eigen-decomposition.
# Output: results/03_meff/<version>/tasks/meff_summary_<model>_<phenotype>_<sample set>.csv
# Usage: Rscript scripts/03_meff/01_meff_task.R <deletion|duplication|mirror> <phenotype> <all|unrelated>
suppressMessages(library(data.table))
suppressMessages(library(Matrix))
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; an <- cfg$pipeline$analysis
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript scripts/03_meff/01_meff_task.R <model> <phenotype> <all|unrelated>")
CNV_TYPE <- args[1]; PHENO <- args[2]; TAG <- args[3]
stopifnot(CNV_TYPE %in% an$models, PHENO %in% an$phenotypes, TAG %in% c("all", "unrelated"))

ver <- cfg$pipeline$callset_version
MAT_DIR      <- file.path(pa$results_root, "02_matrices", ver)
MATRIX_FILE  <- file.path(MAT_DIR, paste0("combined_", CNV_TYPE, "_matrix_shared_probes_fixed.rds"))
MIN_CARRIERS <- as.integer(an$min_carriers)
TARGET_VAR   <- an$simplem_target_variance
OUTPUT_DIR   <- file.path(pa$results_root, "03_meff", ver, "tasks")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
PHENO_FILE   <- pa$phenotype_rds
QC_FILE      <- pa$qc_gwas_file
if (!file.exists(MATRIX_FILE)) stop("Matrix file not found: ", MATRIX_FILE)

matrix_basename <- basename(MATRIX_FILE)
matrix_info <- list(cnv_type = CNV_TYPE, phase_type = "combined", matrix_type = "unpruned")
task_id  <- sprintf("%s_%s_%s", CNV_TYPE, PHENO, TAG)
log_file <- file.path(OUTPUT_DIR, sprintf("meff_%s.log", task_id))
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(msg) {
  ts   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- sprintf("[%s] %s", ts, msg)
  cat(line, "\n")
  cat(line, "\n", file = log_file, append = TRUE)
}

log_msg(sprintf("Starting M_eff computation for matrix: %s | %s | %s", MATRIX_FILE, PHENO, TAG))
log_msg(sprintf("Min carriers threshold: %d (carrier = value != 0)", MIN_CARRIERS))

# =========================
# HELPERS
# =========================

norm_chr <- function(dt) {
  dt[, chr := as.character(chr)]
  dt[chr %in% c("X","x"), chr := "23"]
  dt[chr %in% c("Y","y"), chr := "24"]
  suppressWarnings(dt[, chr := as.integer(chr)])
  dt
}

meff_simpleM_chr <- function(X_chr, target = 0.995, cnv_type = "mirror") {
  sds  <- suppressWarnings(apply(X_chr, 2, sd, na.rm = TRUE))
  keep <- which(is.finite(sds) & sds > 0)
  if (length(keep) <= 1L) return(length(keep))
  X_chr <- as.matrix(X_chr[, keep, drop = FALSE])

  # impute NAs with column mode
  for (j in seq_len(ncol(X_chr))) {
    nas <- is.na(X_chr[, j])
    if (any(nas)) {
      vals <- X_chr[!nas, j]
      if (length(vals) > 0) {
        mode_val <- as.numeric(names(sort(table(vals), decreasing = TRUE))[1])
        X_chr[nas, j] <- mode_val
      } else {
        X_chr[nas, j] <- 0
      }
    }
  }

  C  <- suppressWarnings(cor(X_chr))
  C[!is.finite(C)] <- 0
  C  <- (C + t(C)) / 2
  ev <- tryCatch(eigen(C, symmetric = TRUE, only.values = TRUE)$values,
                 error = function(e) NULL)
  if (is.null(ev)) return(ncol(X_chr))
  ev  <- pmax(ev, 0)
  tot <- sum(ev)
  if (!is.finite(tot) || tot <= 0) return(0L)
  k <- which(cumsum(ev) / tot >= target)[1]
  if (is.na(k)) k <- ncol(X_chr)
  as.integer(k)
}

compute_meff <- function(X, pmap, target = 0.995, cnv_type = "mirror") {
  chrs      <- sort(unique(pmap$chr[pmap$chr %in% 1:22]))
  k_per_chr <- integer(length(chrs)); names(k_per_chr) <- chrs

  total_start <- Sys.time()
  for (cc in chrs) {
    ids  <- pmap[chr == cc, probe_id]
    cols <- intersect(ids, colnames(X))
    if (length(cols) == 0L) {
      k_per_chr[as.character(cc)] <- 0L
      log_msg(sprintf("Chr %d: 0 probes", cc))
      next
    }
    Xi <- X[, cols, drop = FALSE]
    if (inherits(Xi, "dgCMatrix")) Xi <- as.matrix(Xi)
    k <- meff_simpleM_chr(Xi, target, cnv_type)
    k_per_chr[as.character(cc)] <- k
    log_msg(sprintf("Chr %2d: %4d probes -> %d effective tests", cc, length(cols), k))
  }
  total_end <- Sys.time()
  list(
    Meff            = max(1L, sum(k_per_chr)),
    k_per_chr       = k_per_chr,
    total_time_mins = as.numeric(difftime(total_end, total_start, units = "mins"))
  )
}

get_sample_info <- function(matrix_info) {
  key1 <- fread(pa$array1_key)
  key2 <- fread(pa$array2_key)
  setnames(key1, c("ccvid","xIID_cnv"))
  setnames(key2, c("ccvid","xIID_cnv"))
  key1$gwas_array_type <- "Axiom_CKB1"
  key2$gwas_array_type <- "Axiom_CKB2"
  combined_key <- rbind(key1, key2)
  probe_map_file <- file.path(MAT_DIR, "combined_probe_map_shared_probes.rds")
  return(list(key = combined_key, probe_map_file = probe_map_file, is_combined = TRUE))
}

# =========================
# MAIN
# =========================

# Load matrix
log_msg("Loading matrix...")
X <- readRDS(MATRIX_FILE)
log_msg(sprintf("Matrix dimensions: %d samples x %d probes", nrow(X), ncol(X)))

# Get sample info based on matrix type
sample_info <- get_sample_info(matrix_info)
key       <- sample_info$key
probe_map <- norm_chr(as.data.table(readRDS(sample_info$probe_map_file)))
log_msg(sprintf("Loaded probe map: %d probes", nrow(probe_map)))

# Load phenotype data
pheno_data <- readRDS(PHENO_FILE)
setDT(pheno_data)

# ---- QC file: passing samples (same filter as scripts/04_gwas/01_gwas_task.R) ----
log_msg("Loading QC master file...")
qc_master <- fread(QC_FILE)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_samples <- qc_master[qc_pass == TRUE]$GID
log_msg(sprintf("QC-passing samples in master file: %d", length(qc_pass_samples)))
log_msg(sprintf("QC-failing samples excluded: %d", sum(!qc_master$qc_pass)))
# ---------------------------------------------------------------------

if (!PHENO %in% names(pheno_data)) stop("Phenotype not in phenotype file: ", PHENO)
RINT_PHENOS <- PHENO

all_results  <- list()
summary_data <- list()

for (pheno in RINT_PHENOS) {
  for (use_unrelated in (TAG == "unrelated")) {
    tag <- if (use_unrelated) "unrelated" else "all"
    log_msg(sprintf("Processing %s | %s", pheno, tag))

    # Filter phenotype data
    df <- pheno_data[!is.na(get(pheno))]
    if (use_unrelated) df <- df[unrelated == 1]

    # Merge with sample key(s)
    if (sample_info$is_combined) {
      df1      <- df[gwas_array_type == "Axiom_CKB1"]
      df2      <- df[gwas_array_type == "Axiom_CKB2"]
      key1     <- key[gwas_array_type == "Axiom_CKB1"]
      key2     <- key[gwas_array_type == "Axiom_CKB2"]
      merged1  <- merge(df1, key1, by = "ccvid", all = FALSE)
      merged2  <- merge(df2, key2, by = "ccvid", all = FALSE)
      merged   <- rbind(merged1, merged2)
    } else {
      merged <- merge(df, key, by = "ccvid", all = FALSE)
    }

    common <- intersect(merged$xIID_cnv, rownames(X))
    if (length(common) == 0) {
      log_msg(sprintf("Skipping %s (%s): no overlapping samples", pheno, tag))
      next
    }

    # ---- Apply QC filter from master QC file ----
    before_n <- length(common)
    common   <- intersect(common, qc_pass_samples)
    log_msg(sprintf("Removed %d QC-failing samples; %d remaining after QC filter",
                    before_n - length(common), length(common)))

    if (length(common) == 0) {
      log_msg(sprintf("Skipping %s (%s): no samples remaining after QC filter", pheno, tag))
      next
    }
    # ---------------------------------------------

    Xi <- X[common, , drop = FALSE]

    # Filter for minimum carriers
    carrier_counts <- colSums(Xi != 0, na.rm = TRUE)
    keep_probes    <- carrier_counts >= MIN_CARRIERS
    n_before       <- ncol(Xi)
    Xi             <- Xi[, keep_probes, drop = FALSE]
    n_after        <- ncol(Xi)
    log_msg(sprintf("Filtered to %d probes with >=%d carriers (from %d total)",
                    n_after, MIN_CARRIERS, n_before))

    if (n_after == 0) {
      log_msg(sprintf("Skipping %s (%s): no probes with >=%d carriers", pheno, tag, MIN_CARRIERS))
      next
    }

    res <- compute_meff(Xi, probe_map, target = TARGET_VAR, cnv_type = matrix_info$cnv_type)

    alpha_eff <- 0.05 / res$Meff
    row <- data.frame(
      Matrix_File      = matrix_basename,
      CNV_Type         = matrix_info$cnv_type,
      Phase_Type       = matrix_info$phase_type,
      Matrix_Type      = matrix_info$matrix_type,
      Phenotype        = pheno,
      Relatedness      = tag,
      Min_Carriers     = MIN_CARRIERS,
      N_Samples        = nrow(Xi),
      N_Probes         = ncol(Xi),
      Meff             = res$Meff,
      Alpha_0.05       = sprintf("%.3e", alpha_eff),
      Minus_log10_p05  = sprintf("%.2f", -log10(alpha_eff)),
      Time_mins        = sprintf("%.1f", res$total_time_mins),
      stringsAsFactors = FALSE
    )
    summary_data[[length(summary_data) + 1]] <- row
    all_results[[paste(pheno, tag, sep = "_")]] <-
      list(meta = row, k_per_chr = res$k_per_chr)
  }
}

# Save results (one row per task)
if (!length(summary_data)) stop("No result for ", task_id)
summary_table <- do.call(rbind, summary_data)
summary_csv   <- file.path(OUTPUT_DIR, sprintf("meff_summary_%s.csv", task_id))
results_rds   <- file.path(OUTPUT_DIR, sprintf("meff_results_%s.rds", task_id))

write.csv(summary_table, summary_csv, row.names = FALSE)
saveRDS(all_results, results_rds)

log_msg(sprintf("Saved summary table: %s", basename(summary_csv)))
log_msg(sprintf("Saved detailed results: %s", basename(results_rds)))
log_msg("=== M_EFF COMPUTATION COMPLETE ===")
