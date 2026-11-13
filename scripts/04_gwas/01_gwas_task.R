#!/usr/bin/env Rscript
# Stage 4: probe-level association test for one model x phenotype x sample set.
# For every probe with at least min_carriers carriers, the inverse rank-normal transformed trait is
# regressed on CNV dosage (deletion: 1/0; duplication: 1/0; mirror: -1/0/+1), adjusting for age, age^2,
# sex, recruitment region, array and 11 genetic principal components (section 8). Samples failing the
# CNV QC in paths$qc_gwas_file are removed first.
# Output: results/04_gwas/<version>/<all|unrelated>/gwas_linear_unpruned_combined_<model>_<phenotype>_mincar<N>_<sample set>.csv
# Usage: Rscript scripts/04_gwas/01_gwas_task.R <deletion|duplication|mirror> <phenotype> <all|unrelated>
# Sensitivity analyses (optional environment variables; the main analysis sets none of them):
#   EXTRA_COVARS  comma-separated extra covariates (e.g. hours_since_last_ate_x10); samples missing any are dropped
#   PROBE_LIST    file with one probe ID per line; only these probes are tested (after the carrier filter)
#   OUT_TAG       required with either of the above; output goes to <sample set>/sensitivity/ with _<OUT_TAG> appended
#   MIN_CARRIERS  sensitivity mode only: carrier filter for this run instead of the config value (e.g. 1 to test lead
#                 probes that have fewer than 15 carriers in the unrelated subset)
#   EXCLUDE_SAMPLES sensitivity mode only: file with one CNV sample ID per line, removed after the QC filter
library(data.table)
library(Matrix)
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; an <- cfg$pipeline$analysis; ver <- cfg$pipeline$callset_version
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("Usage: Rscript scripts/04_gwas/01_gwas_task.R <model> <phenotype> <all|unrelated>")
cnv_type      <- args[1]
phenotype     <- args[2]
relatedness   <- args[3]
analysis_type <- "linear"
additional_carrier_filter <- as.integer(an$min_carriers)
matrix_type   <- "unpruned"
mat_dir    <- file.path(pa$results_root, "02_matrices", ver)
matrix_path   <- file.path(mat_dir, paste0("combined_", cnv_type, "_matrix_shared_probes_fixed.rds"))

# ---- Validate ----
if (!cnv_type %in% c("duplication","deletion","mirror")) stop("cnv_type must be 'duplication','deletion','mirror'")
if (!relatedness %in% c("all","unrelated")) stop("relatedness must be 'all' or 'unrelated'")
if (!phenotype %in% an$phenotypes) stop("phenotype not in config/pipeline.yaml analysis$phenotypes")

analysis_suffix <- if (!is.null(additional_carrier_filter)) paste0("_mincar", additional_carrier_filter) else ""

# ---- Directories ----
extra_covs <- strsplit(Sys.getenv("EXTRA_COVARS"), ",")[[1]]; extra_covs <- trimws(extra_covs[nzchar(extra_covs)])
probe_list <- Sys.getenv("PROBE_LIST"); out_tag <- Sys.getenv("OUT_TAG")
is_sens    <- length(extra_covs) > 0 || nzchar(probe_list)
if (is_sens && !nzchar(out_tag)) stop("OUT_TAG must be set with EXTRA_COVARS or PROBE_LIST")
excl_file  <- Sys.getenv("EXCLUDE_SAMPLES")
if (nzchar(excl_file) && !is_sens) stop("EXCLUDE_SAMPLES can only be used in sensitivity mode (with PROBE_LIST or EXTRA_COVARS and OUT_TAG)")
output_dir <- file.path(pa$results_root, "04_gwas", ver, relatedness)
if (is_sens) output_dir <- file.path(output_dir, "sensitivity")
mc_override <- Sys.getenv("MIN_CARRIERS")
if (nzchar(mc_override)) {
  if (!is_sens) stop("MIN_CARRIERS can only be used in sensitivity mode (with PROBE_LIST or EXTRA_COVARS and OUT_TAG)")
  additional_carrier_filter <- as.integer(mc_override)
}
log_dir    <- file.path(output_dir, "logs")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir,    recursive = TRUE, showWarnings = FALSE)

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
log_file  <- file.path(log_dir, paste0("gwas_",analysis_type,"_",matrix_type,"_combined_",cnv_type,"_",phenotype,analysis_suffix,"_",relatedness,
                                      if (is_sens) paste0("_", out_tag) else "", ".log"))
if (file.exists(log_file)) file.remove(log_file)

log_msg <- function(msg) {
  ts   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[",ts,"] ",msg)
  cat(line,"\n")
  cat(line,"\n", file=log_file, append=TRUE)
}

log_msg(paste("Starting",analysis_type,"Combined CNV INTERSECTION ANALYSIS (",matrix_type,") for",cnv_type,"phenotype:",phenotype))
log_msg(paste("Output directory:", output_dir))
log_msg(paste("Log directory:", log_dir))

# ---- 1) Phenotype ----
pheno_data <- readRDS(pa$phenotype_rds)
setDT(pheno_data)

use_unrelated <- relatedness == "unrelated"
if (use_unrelated) {
  pheno_data <- pheno_data[unrelated == 1]
  log_msg(paste("Filtered to unrelated samples:", nrow(pheno_data), "remaining"))
  relatedness_suffix <- "_unrelated"
} else {
  log_msg("Using all samples (related and unrelated)")
  relatedness_suffix <- "_related"
}

if (!phenotype %in% names(pheno_data)) stop("Phenotype column missing")
if (analysis_type=="linear"   && !is.numeric(pheno_data[[phenotype]])) stop("Phenotype must be numeric for linear analysis")
if (analysis_type=="logistic" && !is.numeric(pheno_data[[phenotype]])) stop("Phenotype must be coded 0/1 for logistic analysis")

# ---- 2) Sample key for both phases ----
key1 <- fread(pa$array1_key)
key2 <- fread(pa$array2_key)
setnames(key1, c("ccvid","xIID_cnv"))
setnames(key2, c("ccvid","xIID_cnv"))

pheno_filtered1 <- pheno_data[gwas_array_type=="Axiom_CKB1" & !is.na(get(phenotype))]
pheno_filtered2 <- pheno_data[gwas_array_type=="Axiom_CKB2" & !is.na(get(phenotype))]

pheno_with_cnv_id1 <- merge(pheno_filtered1, key1, by="ccvid", all=FALSE)
pheno_with_cnv_id2 <- merge(pheno_filtered2, key2, by="ccvid", all=FALSE)

pheno_with_cnv_id <- rbind(pheno_with_cnv_id1, pheno_with_cnv_id2)
setkey(pheno_with_cnv_id, xIID_cnv)

log_msg(paste("Combined phenotype data from both phases:", nrow(pheno_with_cnv_id), "samples"))

# ---- 3) CNV matrix ----
cnv_file <- matrix_path
if (!file.exists(cnv_file)) stop("Matrix not found: ", cnv_file)

cnv_mat <- readRDS(cnv_file)
log_msg(paste("Loaded combined CNV matrix:", basename(cnv_file)))
log_msg(paste("Matrix dimensions:", nrow(cnv_mat), "samples x", ncol(cnv_mat), "probes"))

# ---- 4) Load combined probe map ----
probe_map_file <- file.path(mat_dir, "combined_probe_map_shared_probes.rds")
probe_map <- readRDS(probe_map_file)
log_msg(paste("Loaded combined probe map with", nrow(probe_map), "probes"))

# ---- 5) Align samples ----
common <- intersect(pheno_with_cnv_id$xIID_cnv, rownames(cnv_mat))
if (length(common)==0) stop("No overlapping samples")
log_msg(paste("Overlapping samples before QC filter:", length(common)))

# ---- 5b) Remove QC-failing samples from master QC file ----
qc_master <- fread(pa$qc_gwas_file)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_samples <- qc_master[qc_pass == TRUE]$GID

before_n <- length(common)
common   <- intersect(common, qc_pass_samples)
log_msg(paste("Removed", before_n - length(common), "QC-failing samples after cross-referencing QC master file"))
log_msg(paste("Samples remaining after QC filter:", length(common)))

if (length(common)==0) stop("No samples remaining after QC filter")
if (nzchar(excl_file)) {
  ex_ids <- trimws(readLines(excl_file, warn = FALSE)); before_n <- length(common)
  common <- setdiff(common, ex_ids)
  log_msg(paste("EXCLUDE_SAMPLES:", basename(excl_file), "| removed", before_n - length(common), "samples | remaining", length(common)))
}

final_pheno <- pheno_with_cnv_id[J(common)]
final_cnv   <- cnv_mat[common,,drop=FALSE]

# ---- 5a) Sensitivity only: extra covariates (complete cases) ----
if (length(extra_covs)) {
  miss <- setdiff(extra_covs, names(final_pheno)); if (length(miss)) stop("Extra covariates not in phenotype file: ", paste(miss, collapse = ", "))
  keep_cc <- complete.cases(final_pheno[, ..extra_covs])
  log_msg(paste("Extra covariates:", paste(extra_covs, collapse = ", "), "| dropped", sum(!keep_cc), "samples missing them"))
  final_pheno <- final_pheno[keep_cc]; final_cnv <- final_cnv[keep_cc, , drop = FALSE]
}

# ---- 5c) Optional: re-RINT the trait within this analysis sample (config analysis$rint) ----
if (identical(an$rint, "recompute")) {
  raw_col <- sub("_rint$", "", phenotype)
  if (!raw_col %in% names(final_pheno)) stop("rint = recompute needs the raw trait column ", raw_col)
  final_pheno[, (phenotype) := qnorm((frank(get(raw_col), ties.method = "average") - 0.5) / .N)]
  log_msg(paste("Re-RINT of", raw_col, "within the analysis sample (n =", nrow(final_pheno), ")"))
}

# ---- 6) Create phase info from gwas_array_type ----
final_pheno[, phase := ifelse(gwas_array_type == "Axiom_CKB1", "phase12", "phase3")]

log_msg("Created phase information directly from gwas_array_type")
log_msg("Phase mapping: Axiom_CKB1 -> phase12, Axiom_CKB2 -> phase3")

final_phase_counts <- table(final_pheno$phase, useNA="ifany")
log_msg("Final sample distribution by phase:")
for (phase_name in names(final_phase_counts)) {
  log_msg(paste("  ", phase_name, ":", final_phase_counts[phase_name], "samples"))
}

missing_phase <- sum(is.na(final_pheno$phase))
if (missing_phase > 0) {
  stop(paste("ERROR: Still have", missing_phase, "samples with missing phase information"))
} else {
  log_msg("All samples have phase information (no missing data)")
}

# ---- 7) Carrier count filter ----
if (!is.null(additional_carrier_filter)) {
  carriers  <- if (cnv_type %in% c("duplication","deletion")) colSums(final_cnv==1, na.rm=TRUE) else colSums(final_cnv!=0, na.rm=TRUE)
  keep      <- carriers >= additional_carrier_filter
  final_cnv <- final_cnv[,keep,drop=FALSE]
  log_msg(paste("Applied additional carrier filter >=", additional_carrier_filter, ": kept", sum(keep), "of", length(keep), "probes"))
}

# ---- 7b) Sensitivity only: restrict to a probe list ----
if (nzchar(probe_list)) {
  pl <- trimws(readLines(probe_list)); pl <- pl[nzchar(pl)]
  final_cnv <- final_cnv[, colnames(final_cnv) %in% pl, drop = FALSE]
  log_msg(paste("Probe list:", length(pl), "probes requested;", ncol(final_cnv), "pass the carrier filter and are tested"))
}

# ---- 8) Covariates ----
final_pheno[, age2 := age^2]
required_covs <- c("sex","age","age2","region_code","phase",paste0("PC",1:11), extra_covs)
log_msg("Including age^2 in covariates for all traits")

available_covs <- intersect(required_covs, names(final_pheno))
log_msg(paste("Available covariates:", paste(available_covs, collapse=", ")))

categorical_covs <- c("region_code","phase")
continuous_covs  <- setdiff(available_covs, categorical_covs)

covar_terms <- character(0)
if ("region_code" %in% available_covs) covar_terms <- c(covar_terms, "as.factor(region_code)")
if ("phase"       %in% available_covs) covar_terms <- c(covar_terms, "as.factor(phase)")
if (length(continuous_covs) > 0)       covar_terms <- c(covar_terms, continuous_covs)

covar_formula <- if (length(covar_terms) > 0) as.formula(paste("~", paste(covar_terms, collapse=" + "))) else as.formula("~ 1")

covs <- model.matrix(covar_formula, data=final_pheno)
log_msg(paste("Covariate formula:", deparse(covar_formula)))
log_msg(paste("Covariate matrix dimensions:", nrow(covs), "x", ncol(covs)))

# ---- 9) Association test ----
log_msg(paste("Running association tests on", ncol(final_cnv), "probes..."))

assoc_probe <- function(j) {
  tryCatch({
    cnv_vec <- final_cnv[,j]
    mask    <- !is.na(cnv_vec)
    y       <- final_pheno[[phenotype]][mask]
    Xcov    <- covs[mask,-1,drop=FALSE]
    df      <- data.frame(y=y, cnv_vec=cnv_vec[mask], Xcov)

    if (analysis_type=="linear") {
      fit <- lm(y ~ cnv_vec + ., data=df)
      s   <- summary(fit)$coefficients["cnv_vec",]
      pval <- s["Pr(>|t|)"]; beta <- s["Estimate"]; se <- s["Std. Error"]
    } else {
      fit  <- glm(y ~ cnv_vec + ., data=df, family=binomial)
      s    <- summary(fit)$coefficients["cnv_vec",]
      pval <- s["Pr(>|z|)"]; beta <- s["Estimate"]; se <- s["Std. Error"]
    }

    n_car <- if (cnv_type %in% c("duplication","deletion")) sum(cnv_vec==1, na.rm=TRUE) else sum(cnv_vec!=0, na.rm=TRUE)
    n_non <- sum(mask)

    data.table(probe_id=colnames(final_cnv)[j], beta=beta, se=se, pval=pval,
               af=n_car/n_non, n_carriers=n_car, n_analyzed=n_non, n_excluded=sum(!mask), error=NA)
  }, error=function(e) {
    data.table(probe_id=colnames(final_cnv)[j], beta=NA, se=NA, pval=NA,
               af=NA, n_carriers=NA, n_analyzed=NA, n_excluded=NA, error=e$message)
  })
}

n_probes  <- ncol(final_cnv)
batch_sz  <- 1000
res_list  <- vector("list", ceiling(n_probes/batch_sz))
for (b in seq_along(res_list)) {
  i1 <- (b-1)*batch_sz+1; i2 <- min(b*batch_sz, n_probes)
  res_list[[b]] <- rbindlist(lapply(i1:i2, assoc_probe))
  if (b %% 10 == 0) log_msg(paste("Processed", i2, "of", n_probes, "probes"))
}
res_dt <- rbindlist(res_list)
setDT(probe_map)
res_dt <- merge(res_dt, probe_map[,.(probe_id,chr,pos)], by="probe_id", all.x=TRUE)

# ---- 10) Summary ----
log_msg("=== RESULTS SUMMARY ===")
log_msg(paste("Total probes analyzed:", nrow(res_dt)))
log_msg(paste("Analysis errors:", sum(!is.na(res_dt$error))))
log_msg(paste("Successful tests:", sum(is.na(res_dt$error))))
valid_pvals <- res_dt[is.na(error) & !is.na(pval), pval]
if (length(valid_pvals) > 0) {
  log_msg(paste("P-value range:", sprintf("%.2e", min(valid_pvals)), "to", sprintf("%.2e", max(valid_pvals))))
  log_msg(paste("Significant at p<0.05:", sum(valid_pvals < 0.05)))
  log_msg(paste("Significant at p<0.01:", sum(valid_pvals < 0.01)))
  log_msg(paste("Significant at p<0.001:", sum(valid_pvals < 0.001)))
}
log_msg(paste("Sample composition:", ifelse(use_unrelated, "unrelated only", "all samples (related + unrelated)")))
log_msg("Array version included as covariate in all analyses")

# ---- 11) Save ----
out_csv <- file.path(output_dir, paste0("gwas_",analysis_type,"_",matrix_type,"_combined_",cnv_type,"_",phenotype,analysis_suffix,"_",relatedness,
                                         if (is_sens) paste0("_", out_tag) else "", ".csv"))
fwrite(res_dt, out_csv)
log_msg(paste("Results written to:", out_csv))
log_msg("=== ANALYSIS COMPLETE ===")
cat("=== ANALYSIS COMPLETE ===\n")
