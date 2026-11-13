#!/usr/bin/env Rscript
# Stage 10: deletion and duplication carriers at the lead probe of every stage 5 association, for the carrier
# column of Table 1. Stage 5 counts only the carriers the model tests, so for a single-type model the other type
# is counted here, in the same sample as stage 5 (participants with the trait, in that model's matrix and passing
# CNV QC). Counts come from the mirror matrix at the lead probe (-1 deletion, +1 duplication), or from the
# other-type matrix if the probe is not in the mirror matrix. Check: the mirror matrix must reproduce the stage 5
# count of the tested type at every lead probe.
# Aggregate output only. Output: results/10_paper/<version>/table1_carriers_both_types.csv
suppressPackageStartupMessages(library(Matrix))
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
od <- file.path(pa$results_root, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
mat_dir <- file.path(pa$results_root, "02_matrices", ver)
find_matrix_file <- function(ct) {                                  # as 02_cnvr_annotation.R
  f <- file.path(mat_dir, paste0("combined_", ct, "_matrix_shared_probes_fixed.rds")); if (file.exists(f)) return(f)
  fs <- list.files(mat_dir, pattern = paste0("^combined_", ct, "_matrix_shared_probes.*\\.rds$"), full.names = TRUE)
  fs <- fs[!grepl("pruned", fs)]; if (!length(fs)) stop("No unpruned matrix for ", ct); fs[which.max(file.info(fs)$mtime)]
}
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x <- x[, .(phenotype, cnv_type, lead_probe, chr, lead_pos, n_del, n_dup, total_n)]
cat("Stage 5 model-specific rows:", nrow(x), "\n")

# phenotypes, keys and master QC exactly as stage 5 (all participants)
ph <- readRDS(pa$phenotype_rds); setDT(ph)
k1 <- fread(pa$array1_key); k2 <- fread(pa$array2_key); setnames(k1, c("ccvid", "xIID_cnv")); setnames(k2, c("ccvid", "xIID_cnv"))
ph_all <- rbindlist(list(merge(ph[gwas_array_type == "Axiom_CKB1"], k1, by = "ccvid"), merge(ph[gwas_array_type == "Axiom_CKB2"], k2, by = "ccvid")),
                    use.names = TRUE, fill = TRUE)
qc <- fread(pa$qc_gwas_file); qc_pass <- qc[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]$GID
pheno_col <- function(p) { b <- sub("_rint$", "", p); if (b %in% names(ph_all)) b else if (paste0(b, "_rint") %in% names(ph_all)) paste0(b, "_rint") else stop("No column for ", p) }

mats <- list(); ids <- list()
for (ct in c("deletion", "duplication", "mirror")) { mats[[ct]] <- readRDS(find_matrix_file(ct)); ids[[ct]] <- rownames(mats[[ct]])
  if (!length(ids[[ct]])) stop("No row names in the ", ct, " matrix") }
other <- c(deletion = "duplication", duplication = "deletion")
out <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
  r <- x[i]; col <- pheno_col(r$phenotype)
  smp <- intersect(intersect(ph_all[!is.na(get(col)), xIID_cnv], ids[[r$cnv_type]]), qc_pass)     # stage 5 sample of this model
  src <- "mirror"; nd <- NA_integer_; nu <- NA_integer_
  if (r$lead_probe %in% colnames(mats$mirror)) {
    v <- as.numeric(mats$mirror[intersect(smp, ids$mirror), r$lead_probe]); nd <- sum(v < 0, na.rm = TRUE); nu <- sum(v > 0, na.rm = TRUE)
  } else if (r$cnv_type != "mirror" && r$lead_probe %in% colnames(mats[[other[r$cnv_type]]])) {
    o <- other[[r$cnv_type]]; v <- as.numeric(mats[[o]][intersect(smp, ids[[o]]), r$lead_probe]); src <- paste0(o, " matrix")
    if (o == "deletion") { nd <- sum(v > 0, na.rm = TRUE); nu <- r$n_dup } else { nu <- sum(v > 0, na.rm = TRUE); nd <- r$n_del }
  } else src <- "not in mirror or other-type matrix"
  data.table(phenotype = r$phenotype, cnv_type = r$cnv_type, lead_probe = r$lead_probe, n_sample = length(smp), stage5_n_del = r$n_del,
             stage5_n_dup = r$n_dup, n_del_at_lead = nd, n_dup_at_lead = nu, source = src)
}))
# checks: the tested type must match stage 5 exactly
out[, check_tested := fifelse(cnv_type == "deletion", n_del_at_lead == stage5_n_del, fifelse(cnv_type == "duplication", n_dup_at_lead == stage5_n_dup,
                                                                                             n_del_at_lead == stage5_n_del & n_dup_at_lead == stage5_n_dup))]
print(out[, .(phenotype, cnv_type, lead_probe, n_sample, stage5_n_del, stage5_n_dup, n_del_at_lead, n_dup_at_lead, source, check_tested)])
nbad <- out[is.na(check_tested) | check_tested == FALSE, .N]
cat("\nCHECK tested-type counts reproduce stage 5:", nrow(out) - nbad, "of", nrow(out), if (nbad) "-> FAIL" else "-> PASS", "\n")
fwrite(out, file.path(od, "table1_carriers_both_types.csv"))
cat("Saved:", file.path(od, "table1_carriers_both_types.csv"), "\n")
if (nbad) quit(status = 1)
