#!/usr/bin/env Rscript
# Stage 7, step 3: look-up of each CKB association in the UK Biobank CNV-GWAS summary statistics of
# Auwerx et al. 2022 (SD3 mirror, SD4 duplication-only, SD5 deletion-only; hg19) (Table S8).
# Same model, matched or proxy trait, UK Biobank probes tested inside the CKB CNVR. Reports the UK Biobank probe
# nearest the CKB lead probe (the primary comparison) and the smallest-P UK Biobank probe in the CNVR.
# Output: results/07_external/<version>/ukb/table_s8_ukb_lookup.csv
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
if (is.null(pa$ukb_auwerx_dir)) stop("Set ukb_auwerx_dir in config/paths.yaml")
od <- file.path(pa$results_root, "07_external", ver, "ukb"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
files <- c(mirror = "SD3_Auwerx_2021_CNV_GWAS_SummaryStats_Mirror.txt", duplication = "SD4_Auwerx_2021_CNV_GWAS_SummaryStats_DuplicationOnly.txt",
           deletion = "SD5_Auwerx_2021_CNV_GWAS_SummaryStats_DeletionOnly.txt")
trait_map <- data.table(trait = c("bmi_calc", "weight_kg", "standing_height_mm", "sitting_height_mm", "dbp_mean", "sbp_mean", "heart_rate_mean", "random_glucose", "fat_percent", "bmi_at_25"),
                        ukb_trait = c("BMI", "weight", "height", "height", "diastolic_BP", "systolic_BP", "heart_rate", "glucose", "body_fat_mass", "BMI"),
                        match = c("exact", "exact", "exact", "proxy", "exact", "exact", "exact", "exact", "proxy", "proxy"))
ss <- rbindlist(lapply(names(files), function(m) {
  y <- fread(file.path(pa$ukb_auwerx_dir, files[[m]]), select = c("PHENO", "CHR", "POS", "ID", "OBS_CT", "BETA", "SE", "P"))
  y[PHENO %in% unique(trait_map$ukb_trait)][, model := m] }))
ss[, CHR := as.integer(CHR)]
cat("UKB rows loaded:", nrow(ss), "\n")
chk <- ss[PHENO == "BMI" & CHR == 16 & POS >= 29650000 & POS <= 30200000, .(n_probes = .N, median_beta = round(median(BETA), 3), min_P = signif(min(P), 3)), by = model]
cat("SIGN CHECK, UKB BMI at 16p11.2 BP4-BP5 (expect deletion > 0, duplication < 0, mirror < 0):\n"); print(chk)

x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x[, `:=`(start = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), end = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end), trait = sub("_rint$", "", phenotype))]
sig <- merge(x[, .(cytoband = locus_cytoband_final, trait, model = cnv_type, chr, start, end, lead = lead_pos, ckb_beta = lead_beta, ckb_P = pval)], trait_map, by = "trait", all.x = TRUE)
cat("\nSignals:", nrow(sig), "| without a UKB trait (not looked up):", sig[is.na(ukb_trait), paste(cytoband, trait, model, collapse = "; ")], "\n")
sig <- sig[!is.na(ukb_trait)]
res <- rbindlist(lapply(seq_len(nrow(sig)), function(i) {
  s <- sig[i]; y <- ss[model == s$model & PHENO == s$ukb_trait & CHR == s$chr & POS >= s$start & POS <= s$end]
  base <- s[, .(cytoband, trait, ukb_trait, match, model, ckb_beta = round(ckb_beta, 3), ckb_P = signif(ckb_P, 3))]
  if (!nrow(y)) return(cbind(base, data.table(n_ukb_probes = 0L)))
  near <- y[which.min(abs(POS - s$lead))]; best <- y[which.min(P)]
  cbind(base, data.table(n_ukb_probes = nrow(y), near_dist_kb = round(abs(near$POS - s$lead) / 1e3, 1),
    near_beta = round(near$BETA, 3), near_se = round(near$SE, 3), near_P = signif(near$P, 3), near_same_dir = sign(near$BETA) == sign(s$ckb_beta),
    best_pos = best$POS, best_beta = round(best$BETA, 3), best_se = round(best$SE, 3), best_P = signif(best$P, 3),
    best_P_x_nprobes = signif(pmin(1, best$P * nrow(y)), 3), best_same_dir = sign(best$BETA) == sign(s$ckb_beta), ukb_N = near$OBS_CT))
}), fill = TRUE)
setorder(res, cytoband, trait, model)
cat("\n== Table S8 (near = UKB probe nearest the CKB lead; best = smallest P in the CNVR)\n")
print(res[, .(cytoband, trait, ukb_trait, match, model, ckb_beta, ckb_P, n_ukb_probes, near_dist_kb, near_beta, near_se, near_P, near_same_dir, best_beta, best_P, best_P_x_nprobes, best_same_dir)])
fwrite(res, file.path(od, "table_s8_ukb_lookup.csv")); cat("\nSaved:", file.path(od, "table_s8_ukb_lookup.csv"), "\n")
