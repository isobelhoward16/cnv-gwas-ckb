#!/usr/bin/env Rscript

# Stage 6: array-stratified checks at the lead probe of each main signal. For each signal, fits the
# association model pooled across arrays, separately in each array (CKB1, CKB2), and with a CNV x array
# interaction, using the same covariates as the main analysis.
# Output: results/06_followup/<version>/array_checks/
# Usage: Rscript scripts/06_followup/04_array_validation.R
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
})
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
auwerx_csv     <- file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv")
matrix_type    <- "unpruned"
min_carriers   <- 1L
use_unrelated  <- FALSE
matrix_dir     <- file.path(pa$results_root, "02_matrices", ver)
out_prefix     <- "pooling_validation"

if (!file.exists(auwerx_csv)) stop("Input file missing: ", auwerx_csv)
if (!matrix_type %in% c("pruned", "unpruned")) stop("matrix_type must be pruned|unpruned")
if (is.na(min_carriers) || min_carriers < 1) stop("min_carriers must be >= 1")

log_msg <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), paste0(...)))
}

find_latest_file <- function(files) {
  files[which.max(file.info(files)$mtime)]
}

safe_coef <- function(fit, term) {
  s <- summary(fit)$coefficients
  if (!(term %in% rownames(s))) {
    return(c(Estimate = NA_real_, `Std. Error` = NA_real_, P = NA_real_))
  }
  if ("Pr(>|t|)" %in% colnames(s)) {
    pcol <- "Pr(>|t|)"
  } else if ("Pr(>|z|)" %in% colnames(s)) {
    pcol <- "Pr(>|z|)"
  } else {
    pcol <- NA_character_
  }

  c(
    Estimate = unname(s[term, "Estimate"]),
    `Std. Error` = unname(s[term, "Std. Error"]),
    P = if (!is.na(pcol)) unname(s[term, pcol]) else NA_real_
  )
}

carrier_fun <- function(x, cnv_type) {
  if (cnv_type %in% c("duplication", "deletion")) {
    x == 1
  } else {
    x != 0
  }
}

get_required_covars <- function(phenotype) {
  height_traits <- c("standing_height_mm", "sitting_height_mm")
  is_height <- phenotype %in% height_traits

  c("sex", "age", "age2", "region_code", paste0("PC", 1:11))
}

load_matrix <- function(cnv_type, matrix_type, matrix_dir) {
  if (matrix_type == "pruned") {
    mat_dir <- file.path(matrix_dir, "pruned")

    patt1 <- paste0("^combined_", cnv_type, "_matrix_.*shared_probes.*pruned.*intersection.*\\.rds$")
    files <- list.files(mat_dir, pattern = patt1, full.names = TRUE)

    if (length(files) == 0) {
      patt2 <- paste0("^combined_", cnv_type, "_matrix_.*shared_probes.*pruned.*\\.rds$")
      files <- list.files(mat_dir, pattern = patt2, full.names = TRUE)
    }
    if (length(files) == 0) stop("No pruned matrix found for ", cnv_type)

    fp <- find_latest_file(files)
    log_msg("Using pruned matrix for ", cnv_type, ": ", fp)
    return(readRDS(fp))
  } else {
    fp <- file.path(matrix_dir, paste0("combined_", cnv_type, "_matrix_shared_probes.rds"))
    if (!file.exists(fp)) {
      fp2 <- file.path(matrix_dir, paste0("combined_", cnv_type, "_matrix_shared_probes_fixed.rds"))
      if (!file.exists(fp2)) stop("No unpruned shared matrix found for ", cnv_type)
      fp <- fp2
    }
    log_msg("Using unpruned matrix for ", cnv_type, ": ", fp)
    return(readRDS(fp))
  }
}

# -------------------------------------------------
# 1) Lead probes of the stage 5 signals
# -------------------------------------------------
lead_dt <- fread(auwerx_csv)

required_cols <- c(
  "phenotype", "pheno_col_used", "cnv_type", "lead_probe",
  "lead_pos", "lead_beta", "lead_se", "pval",
  "cnvr_id", "locus_cytoband_final"
)
missing_cols <- setdiff(required_cols, names(lead_dt))
if (length(missing_cols) > 0) {
  stop("Input does not look like the stage 5 CNVR output. Missing columns: ",
       paste(missing_cols, collapse = ", "))
}

lead_dt <- unique(
  lead_dt[, .(
    phenotype,
    pheno_col_used,
    cnv_type,
    cnvr_id,
    chr,
    event_start,
    event_end,
    lead_probe,
    lead_pos,
    pval,
    lead_beta,
    lead_se,
    locus_cytoband_final,
    n_carrier,
    carrier_freq_pct
  )]
)

setorder(lead_dt, phenotype, cnv_type, chr, event_start, pval)

log_msg("Loaded CNVR rows: ", nrow(lead_dt))
log_msg("CNV type counts:\n", paste(capture.output(print(lead_dt[, .N, by = cnv_type])), collapse = "\n"))

n_expected <- 23L   # model-specific associations in stage 5 (all participants)
if (nrow(lead_dt) != n_expected) {
  stop("Expected the ", n_expected, " model-specific associations from the stage 5 output, found: ", nrow(lead_dt))
}

# -------------------------------------------------
# 2) Load phenotype data, keep default related/full sample
# -------------------------------------------------
pheno_path <- pa$phenotype_rds
pheno_data <- readRDS(pheno_path)
setDT(pheno_data)

if (use_unrelated) {
  if (!"unrelated" %in% names(pheno_data)) stop("Column 'unrelated' not found in phenotype data")
  pheno_data <- pheno_data[unrelated == 1]
}
log_msg("Phenotype rows after unrelated filter (use_unrelated=", use_unrelated, "): ", nrow(pheno_data))

key1 <- fread(pa$array1_key)
key2 <- fread(pa$array2_key)
setnames(key1, c("ccvid", "xIID_cnv"))
setnames(key2, c("ccvid", "xIID_cnv"))

ph1 <- pheno_data[gwas_array_type == "Axiom_CKB1"]
ph2 <- pheno_data[gwas_array_type == "Axiom_CKB2"]

ph1 <- merge(ph1, key1, by = "ccvid", all = FALSE)
ph2 <- merge(ph2, key2, by = "ccvid", all = FALSE)
pheno_with_cnv_id <- rbind(ph1, ph2, fill = TRUE)
setkey(pheno_with_cnv_id, xIID_cnv)

pheno_with_cnv_id[, phase := ifelse(gwas_array_type == "Axiom_CKB1", "phase12", "phase3")]
pheno_with_cnv_id[, phase := factor(phase, levels = c("phase12", "phase3"))]

if ("sex" %in% names(pheno_with_cnv_id)) pheno_with_cnv_id[, sex := as.factor(sex)]
if ("region_code" %in% names(pheno_with_cnv_id)) pheno_with_cnv_id[, region_code := as.factor(region_code)]
if (!"age2" %in% names(pheno_with_cnv_id) && "age" %in% names(pheno_with_cnv_id)) {
  pheno_with_cnv_id[, age2 := age^2]
}

log_msg("Phase counts before QC: ", paste(names(table(pheno_with_cnv_id$phase)), table(pheno_with_cnv_id$phase), collapse = "; "))

# -------------------------------------------------
# 3) Apply master QC filtering
# -------------------------------------------------
qc_master_path <- pa$qc_gwas_file
if (!file.exists(qc_master_path)) stop("Master QC file missing: ", qc_master_path)

qc_master <- fread(qc_master_path)
if (!all(c("GID", "LRR_SD", "WF", "BAF_drift", "NumCNV") %in% names(qc_master))) {
  stop("Master QC file missing one or more required columns: GID, LRR_SD, WF, BAF_drift, NumCNV")
}

qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_samples <- unique(qc_master[qc_pass == TRUE, GID])

before_qc <- nrow(pheno_with_cnv_id)
pheno_with_cnv_id <- pheno_with_cnv_id[xIID_cnv %in% qc_pass_samples]
after_qc <- nrow(pheno_with_cnv_id)

log_msg("Applied master QC filter: removed ", before_qc - after_qc, " phenotype rows")
log_msg("Phenotype rows after master QC: ", after_qc)
log_msg("Phase counts after master QC: ", paste(names(table(pheno_with_cnv_id$phase)), table(pheno_with_cnv_id$phase), collapse = "; "))

# -------------------------------------------------
# 4) Matrix cache
# -------------------------------------------------
mat_cache <- new.env(parent = emptyenv())

get_mat <- function(cnv_type) {
  if (!exists(cnv_type, envir = mat_cache, inherits = FALSE)) {
    assign(cnv_type, load_matrix(cnv_type, matrix_type, matrix_dir), envir = mat_cache)
  }
  get(cnv_type, envir = mat_cache, inherits = FALSE)
}

# -------------------------------------------------
# 5) Per-lead validation
# -------------------------------------------------
results <- vector("list", nrow(lead_dt))

for (i in seq_len(nrow(lead_dt))) {
  row <- lead_dt[i]
  phenotype  <- row$phenotype
  cnv_type   <- row$cnv_type
  lead_probe <- row$lead_probe

   pheno_col <- ifelse(grepl("_rint$", phenotype), phenotype, paste0(phenotype, "_rint"))

  if (!(pheno_col %in% names(pheno_with_cnv_id))) {
    results[[i]] <- data.table(
      phenotype = phenotype,
      pheno_col_used = pheno_col,
      cnv_type = cnv_type,
      cnvr_id = row$cnvr_id,
      lead_probe = lead_probe,
      lead_pos = row$lead_pos,
      locus_cytoband_final = row$locus_cytoband_final,
      input_pval = row$pval,
      input_beta = row$lead_beta,
      input_se = row$lead_se,
      matrix_found = TRUE,
      error = paste0("required _rint phenotype column not found: ", pheno_col)
    )
    next
  }

  log_msg("Processing ", i, "/", nrow(lead_dt), ": ", phenotype, " | ", cnv_type, " | ", lead_probe)

  mat <- get_mat(cnv_type)

  if (!(lead_probe %in% colnames(mat))) {
    results[[i]] <- data.table(
      phenotype = phenotype,
      pheno_col_used = pheno_col,
      cnv_type = cnv_type,
      cnvr_id = row$cnvr_id,
      lead_probe = lead_probe,
      lead_pos = row$lead_pos,
      locus_cytoband_final = row$locus_cytoband_final,
      input_pval = row$pval,
      input_beta = row$lead_beta,
      input_se = row$lead_se,
      matrix_found = FALSE,
      error = "lead_probe not found in matrix"
    )
    next
  }

  if (!(pheno_col %in% names(pheno_with_cnv_id))) {
    results[[i]] <- data.table(
      phenotype = phenotype,
      pheno_col_used = pheno_col,
      cnv_type = cnv_type,
      cnvr_id = row$cnvr_id,
      lead_probe = lead_probe,
      lead_pos = row$lead_pos,
      locus_cytoband_final = row$locus_cytoband_final,
      input_pval = row$pval,
      input_beta = row$lead_beta,
      input_se = row$lead_se,
      matrix_found = TRUE,
      error = paste0("phenotype column not found in phenotype data: ", pheno_col)
    )
    next
  }

  ph_sub <- pheno_with_cnv_id[!is.na(get(pheno_col))]
  common <- intersect(ph_sub$xIID_cnv, rownames(mat))

  if (length(common) == 0) {
    results[[i]] <- data.table(
      phenotype = phenotype,
      pheno_col_used = pheno_col,
      cnv_type = cnv_type,
      cnvr_id = row$cnvr_id,
      lead_probe = lead_probe,
      lead_pos = row$lead_pos,
      locus_cytoband_final = row$locus_cytoband_final,
      input_pval = row$pval,
      input_beta = row$lead_beta,
      input_se = row$lead_se,
      matrix_found = TRUE,
      error = "no overlapping samples after phenotype + QC + matrix alignment"
    )
    next
  }

  final_pheno <- ph_sub[J(common)]
  final_cnv   <- mat[common, lead_probe]
  final_cnv   <- as.numeric(final_cnv)

  n_car_total_pre <- sum(carrier_fun(final_cnv, cnv_type), na.rm = TRUE)
  if (n_car_total_pre < min_carriers) {
    results[[i]] <- data.table(
      phenotype = phenotype,
      pheno_col_used = pheno_col,
      cnv_type = cnv_type,
      cnvr_id = row$cnvr_id,
      lead_probe = lead_probe,
      lead_pos = row$lead_pos,
      locus_cytoband_final = row$locus_cytoband_final,
      input_pval = row$pval,
      input_beta = row$lead_beta,
      input_se = row$lead_se,
      matrix_found = TRUE,
      n_samples_total = length(common),
      carriers_total_pre_mask = n_car_total_pre,
      error = paste0("lead probe carrier count below min_carriers (", min_carriers, ")")
    )
    next
  }

  covars <- get_required_covars(phenotype)
  covars <- intersect(covars, names(final_pheno))

  mask_pool <- !is.na(final_cnv)
  df_pool <- copy(final_pheno[mask_pool])
  df_pool[, cnv := final_cnv[mask_pool]]
  df_pool[, y := get(pheno_col)]

  carriers_phase12_after_mask <- sum(carrier_fun(df_pool[phase == "phase12", cnv], cnv_type), na.rm = TRUE)
  carriers_phase3_after_mask  <- sum(carrier_fun(df_pool[phase == "phase3", cnv], cnv_type), na.rm = TRUE)

  excluded_phase12_by_mask <- sum(final_pheno$phase == "phase12") - sum(mask_pool & final_pheno$phase == "phase12")
  excluded_phase3_by_mask  <- sum(final_pheno$phase == "phase3") - sum(mask_pool & final_pheno$phase == "phase3")

  pooled_terms <- c("cnv", covars, "phase")
  pooled_terms <- unique(pooled_terms)
  pooled_terms <- pooled_terms[pooled_terms %in% c("cnv", names(df_pool))]
  pooled_formula <- as.formula(paste("y ~", paste(pooled_terms, collapse = " + ")))

  int_covars <- setdiff(covars, "phase")
  int_terms <- c("cnv * phase", int_covars)
  int_terms <- int_terms[nzchar(int_terms)]
  int_formula <- as.formula(paste("y ~", paste(int_terms, collapse = " + ")))

  strat_terms <- c("cnv", int_covars)
  strat_terms <- unique(strat_terms)
  strat_terms <- strat_terms[strat_terms %in% c("cnv", names(df_pool))]
  strat_formula <- as.formula(paste("y ~", paste(strat_terms, collapse = " + ")))

  pool_beta <- pool_se <- pool_p <- NA_real_
  p12_beta  <- p12_se  <- p12_p  <- NA_real_
  p3_beta   <- p3_se   <- p3_p   <- NA_real_
  int_beta  <- int_se  <- int_p  <- NA_real_
  model_error <- NA_character_

  tryCatch({
    fit_pool <- lm(pooled_formula, data = df_pool)
    s_pool <- safe_coef(fit_pool, "cnv")
    pool_beta <- s_pool["Estimate"]
    pool_se   <- s_pool["Std. Error"]
    pool_p    <- s_pool["P"]

    df_p12 <- df_pool[phase == "phase12"]
    df_p3  <- df_pool[phase == "phase3"]

    if (nrow(df_p12) > 0 && sum(carrier_fun(df_p12$cnv, cnv_type), na.rm = TRUE) > 0) {
      fit_p12 <- lm(strat_formula, data = df_p12)
      s_p12 <- safe_coef(fit_p12, "cnv")
      p12_beta <- s_p12["Estimate"]
      p12_se   <- s_p12["Std. Error"]
      p12_p    <- s_p12["P"]
    }

    if (nrow(df_p3) > 0 && sum(carrier_fun(df_p3$cnv, cnv_type), na.rm = TRUE) > 0) {
      fit_p3 <- lm(strat_formula, data = df_p3)
      s_p3 <- safe_coef(fit_p3, "cnv")
      p3_beta <- s_p3["Estimate"]
      p3_se   <- s_p3["Std. Error"]
      p3_p    <- s_p3["P"]
    }

    if (length(unique(df_pool$phase)) == 2) {
      fit_int <- lm(int_formula, data = df_pool)
      int_name <- grep("^cnv:phase|^phase.*:cnv", rownames(summary(fit_int)$coefficients), value = TRUE)
      if (length(int_name) == 1) {
        s_int <- safe_coef(fit_int, int_name)
        int_beta <- s_int["Estimate"]
        int_se   <- s_int["Std. Error"]
        int_p    <- s_int["P"]
      }
    }

  }, error = function(e) {
    model_error <<- e$message
  })

  results[[i]] <- data.table(
    phenotype = phenotype,
    pheno_col_used = pheno_col,
    cnv_type = cnv_type,
    cnvr_id = row$cnvr_id,
    chr = row$chr,
    event_start = row$event_start,
    event_end = row$event_end,
    lead_probe = lead_probe,
    lead_pos = row$lead_pos,
    locus_cytoband_final = row$locus_cytoband_final,

    input_pval = row$pval,
    input_beta = row$lead_beta,
    input_se = row$lead_se,

    matrix_found = TRUE,
    use_unrelated = use_unrelated,
    matrix_type = matrix_type,
    min_carriers = min_carriers,

    n_samples_total = length(common),
    n_samples_after_mask = nrow(df_pool),

    n_phase12_total = sum(final_pheno$phase == "phase12"),
    n_phase3_total  = sum(final_pheno$phase == "phase3"),

    n_phase12_after_mask = sum(df_pool$phase == "phase12"),
    n_phase3_after_mask  = sum(df_pool$phase == "phase3"),

    carriers_total_pre_mask = n_car_total_pre,
    carriers_phase12_after_mask = carriers_phase12_after_mask,
    carriers_phase3_after_mask  = carriers_phase3_after_mask,

    excluded_phase12_by_mask = excluded_phase12_by_mask,
    excluded_phase3_by_mask  = excluded_phase3_by_mask,

    pool_beta = pool_beta,
    pool_se   = pool_se,
    pool_p    = pool_p,

    p12_beta = p12_beta,
    p12_se   = p12_se,
    p12_p    = p12_p,

    p3_beta = p3_beta,
    p3_se   = p3_se,
    p3_p    = p3_p,

    int_beta = int_beta,
    int_se   = int_se,
    int_p    = int_p,

    error = model_error
  )
}

res_dt <- rbindlist(results, fill = TRUE)

log_msg("Finished validation rows: ", nrow(res_dt))
log_msg("Rows with errors: ", sum(!is.na(res_dt$error)))

# -------------------------------------------------
# 6) Save
# -------------------------------------------------
out_dir <- file.path(pa$results_root, "06_followup", ver, "array_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tag <- ifelse(use_unrelated, "unrelated", "related")

out_file <- file.path(
  out_dir,
  paste0(out_prefix, "_", matrix_type, "_mincar", min_carriers, "_", tag, ".csv")
)

fwrite(res_dt, out_file)
log_msg("Wrote validation results: ", out_file)

valid_pool <- sum(!is.na(res_dt$pool_p))
valid_int  <- sum(!is.na(res_dt$int_p))
log_msg("Valid pooled models: ", valid_pool)
log_msg("Valid interaction models: ", valid_int)
chk <- res_dt[!is.na(pool_beta), max(abs(pool_beta - input_beta))]
log_msg("Check: largest |pooled beta - main GWAS beta| = ", signif(chk, 3), " (expect about 0)")
log_msg("Interaction P < 0.05: ", sum(res_dt$int_p < 0.05, na.rm = TRUE), " of ", valid_int,
        " | Bonferroni (0.05/", valid_int, "): ", sum(res_dt$int_p < 0.05 / valid_int, na.rm = TRUE))
log_msg("Direction CKB1 vs CKB2 concordant: ", sum(sign(res_dt$p12_beta) == sign(res_dt$p3_beta), na.rm = TRUE), " of ",
        sum(!is.na(res_dt$p12_beta) & !is.na(res_dt$p3_beta)))
log_msg("Direction each array vs pooled concordant: CKB1 ", sum(sign(res_dt$p12_beta) == sign(res_dt$pool_beta), na.rm = TRUE),
        " | CKB2 ", sum(sign(res_dt$p3_beta) == sign(res_dt$pool_beta), na.rm = TRUE))
print(res_dt[, .(locus_cytoband_final, phenotype, cnv_type, carriers_phase12_after_mask, carriers_phase3_after_mask,
                 pool_beta = round(pool_beta, 3), p12_beta = round(p12_beta, 3), p3_beta = round(p3_beta, 3), int_p = signif(int_p, 3))])
log_msg("Done.")
