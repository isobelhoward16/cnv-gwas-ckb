#!/usr/bin/env Rscript
# Stage 8, step 1: manifest of the participants passing CNV QC (the denominator for the descriptive
# analyses), with their CNV counts, array phase and region, and counts by phase and region.
# Run with DESCRIPTIVE_PATHS set (see slurm/08_descriptive.sbatch).
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(yaml)
})

# ---- Load config ----
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)
qc_thr <- yaml::read_yaml(file.path(dirname(paths_yaml), "qc_thresholds.yaml"))

gene_dir <- paths$gene_cleaned_merged_dir
qc_file  <- paths$qc_master_file
analysis_csv <- paths$analysis_dataframe_csv
zero_cnv_file <- paths$qc_passed_zero_cnv_file

out_root <- paths$out_root
out_manifest <- file.path(out_root, "data/derived/manifests/postQC_denominator_samples.tsv")
out_counts   <- file.path(out_root, "data/derived/tables/denominator_counts_by_phase_region.csv")

log_file <- file.path(out_root, "logs/01_make_postQC_denominator_manifest.log")
dir.create(dirname(out_manifest), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(out_counts), recursive=TRUE, showWarnings=FALSE)
dir.create(dirname(log_file), recursive=TRUE, showWarnings=FALSE)

sink(log_file, split=TRUE)
cat("=== Build post-QC denominator manifest ===\n")
cat("Started:", format(Sys.time()), "\n\n")

# ---- 1) CKIDs present in final CNV callset (gene_cleaned_merged) ----
cat("Reading CNV files from:", gene_dir, "\n")
cnv_files <- list.files(gene_dir, pattern="\\.txt$", full.names=TRUE)
stopifnot(length(cnv_files) > 0)
cat("Found", length(cnv_files), "CNV files\n")

ckids_with_cnvs <- unique(unlist(lapply(cnv_files, function(f) {
  # column 15 is CKID (confirmed by you)
  suppressWarnings(fread(f, header=FALSE, select=15)[[1]])
})))
ckids_with_cnvs <- ckids_with_cnvs[!is.na(ckids_with_cnvs) & ckids_with_cnvs != ""]
cat("Unique CKIDs with >=1 CNV record:", length(ckids_with_cnvs), "\n\n")

# ---- 2) QC-pass zero-CNV CKIDs (already computed) ----
cat("Reading QC-pass zero-CNV list:", zero_cnv_file, "\n")
zero_cnv <- fread(zero_cnv_file)
setnames(zero_cnv, 1, "ckid")
ckids_zero <- unique(zero_cnv$ckid)
cat("Unique QC-pass zero-CNV CKIDs:", length(ckids_zero), "\n\n")

# ---- 3) Union = post-QC denominator CKIDs ----
postqc_ckids <- unique(c(ckids_with_cnvs, ckids_zero))
cat("Union post-QC denominator CKIDs:", length(postqc_ckids), "\n\n")

# ---- 4) Master QC table: apply sample-QC thresholds + map CK<->G + parse phase/region ----
cat("Reading master QC file:", qc_file, "\n")
qc <- fread(qc_file)

need_cols <- c("GID","GSID","RC_name","File","LRR_SD","WF","BAF_drift","NumCNV")
missing <- setdiff(need_cols, names(qc))
if (length(missing)) stop("QC file missing required columns: ", paste(missing, collapse=", "))

qc_sub <- qc[, .(
  ckid = GID,
  gsid = GSID,
  region_name = RC_name,
  file_path = File,
  LRR_SD, WF, BAF_drift, NumCNV
)]

qc_sub[, qc_pass := (LRR_SD <= qc_thr$LRR_SD_max) &
                  (abs(WF) <= qc_thr$WF_abs_max) &
                  (BAF_drift <= qc_thr$BAF_drift_max)]

qc_pass <- qc_sub[qc_pass == TRUE]
cat("QC-pass samples in QC table:", nrow(qc_pass), "\n")

# Parse region_phase + phase from file path
qc_pass[, region_phase := str_match(file_path, "/temp_affyhmm/([^/]+)/")[,2]]
qc_pass[is.na(region_phase), region_phase := str_match(file_path, "/([^/]+_Phase(?:12|3))/")[,2]]
qc_pass[, phase := fifelse(grepl("Phase12", region_phase), "Phase12",
                    fifelse(grepl("Phase3", region_phase), "Phase3", NA_character_))]

# Keep only samples in our post-QC union set (carriers + zero-CNV)
denom <- qc_pass[ckid %in% postqc_ckids]
cat("QC-pass samples after restricting to (CNV carriers or zero-CNV list):", nrow(denom), "\n\n")

# ---- 5) Join analysis dataframe for sex/unrelated/array_type/etc via gsid ----
cat("Joining analysis dataframe:", analysis_csv, "\n")
adf <- fread(analysis_csv)

stopifnot("ccvid" %in% names(adf))
keep <- intersect(c("ccvid","sex","unrelated","gwas_array_type","region_code"), names(adf))
adf_sub <- adf[, ..keep]
setnames(adf_sub, "ccvid", "gsid")

denom <- merge(denom, adf_sub, by="gsid", all.x=TRUE)

# ---- 6) Final tidy + write ----
denom_out <- unique(denom[, .(
  ckid, gsid,
  phase, region_phase, region_name,
  sex = if ("sex" %in% names(denom)) sex else NA,
  unrelated = if ("unrelated" %in% names(denom)) unrelated else NA,
  gwas_array_type = if ("gwas_array_type" %in% names(denom)) gwas_array_type else NA,
  region_code = if ("region_code" %in% names(denom)) region_code else NA,
  LRR_SD, WF, BAF_drift, NumCNV
)])

fwrite(denom_out, out_manifest, sep="\t")
cat("Wrote manifest:", out_manifest, "\n")
cat("N samples:", nrow(denom_out), "\n\n")

# Counts table (for Table 1 scaffolding)
counts <- denom_out[, .N, by=.(phase, region_phase)][order(phase, region_phase)]
fwrite(counts, out_counts)
cat("Wrote counts table:", out_counts, "\n\n")

cat("Completed:", format(Sys.time()), "\n")
sink()
