#!/usr/bin/env Rscript
# Stage 10: Tables S1, S5 and S6 as plain CSV (the formatted versions are built from these).
#   S1: exclusion regions (config/exclusion_regions_hg19.txt).
#   S5: CNV-wide significant associations, single-type models (deletion, duplication); S6: mirror model.
#       Phenotype | Cytoband | Chr | Sig probe start | Sig probe end | N sig probes | CNV region start (hg19) |
#       CNV region end (hg19) | Lead probe position | P-value | beta | SE | N del | N dup | N carriers |
#       Carrier freq (%) | Total N | CNV type (S5) or Also significant (DEL/DUP) (S6) | GWAS Catalog hit
#   - sig probe start/end and N = the significant probes tagged to the CNVR (stage 5 event_start/end, n_sig_probes_in_event);
#   - CNV region = the r2 >= 0.5 extension of Auwerx et al. 2022 (stage 5 auwerx_event_start/end; the stage 5 CNVR
#     where it is missing), as in Table 1;
#   - P-value: "x10^" with superscript digits, dagger when the experiment-wide threshold is passed;
#   - N del / N dup: "-" for the type a single-type model does not test;
#   - GWAS Catalog hit: "Yes" if >= 1 trait-matched (Table S3 names) Catalog association lies in the CNVR
#     (stage 7 catalog_counts_summary.csv, matched by trait, chromosome and overlapping CNVR); blank if stage 7 has not run;
#   - Cytoband: stage 5 cytoband, "1q21.2" read as 1q21.1 (as in Table 1).
# Also writes the thresholds quoted in the captions. Aggregate output only. Output: results/10_paper/<version>/
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
od <- file.path(pa$results_root, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- S1 ----
reg <- trimws(readLines(file.path(cfg$root, cfg$pipeline$exclusion_regions_file), warn = FALSE)); reg <- reg[reg != "" & !startsWith(reg, "#")]
reg <- sub("[[:space:]].*$", "", reg)
# immune-receptor loci first (IGL, IGH, IGK, TRA), then the rest in config order
imm <- c(IGL = "chr22:22380474-23265085", IGH = "chr14:106052774-107288051", IGK = "chr2:89156874-90274235", TRA = "chr14:22090057-23021075")
if (!all(imm %in% reg)) stop("Immune-receptor regions in the config differ from the expected hg19 coordinates")
if (anyDuplicated(reg)) stop("Duplicate regions in the config")
reg <- c(unname(imm), setdiff(reg, imm))
fwrite(data.table(region_hg19 = reg), file.path(od, "table_s1_exclusion_regions.csv"))
cat("Table S1:", length(reg), "regions (", length(imm), "immune-receptor loci first, then", length(reg) - length(imm), "centromere/telomere)\n")

# ---- S5/S6 need stage 5; stop cleanly if it has not run yet ----
x_file <- file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv")
if (!file.exists(x_file)) { cat("Stage 5 output not found yet: Table S1 written, Tables S5/S6 skipped. Rerun after stage 5.\n"); quit(save = "no", status = 0) }

# ---- S5/S6 ----
x <- fread(x_file)
thr <- fread(file.path(pa$results_root, "03_meff", ver, "meff_thresholds.csv"))[relatedness == "all"]
names_s <- c(bmi_calc = "BMI", bmi_at_25 = "BMI at age 25", standing_height_mm = "Standing height", sitting_height_mm = "Sitting height",
             weight_kg = "Weight", waist_mm = "Waist circumference", hip_mm = "Hip circumference", waist_hip_ratio = "Waist-hip ratio",
             fat_percent = "Body fat percentage", sbp_mean = "Systolic blood pressure", dbp_mean = "Diastolic blood pressure",
             heart_rate_mean = "Resting heart rate", random_glucose = "Random plasma glucose")
x[, ph := sub("_rint$", "", phenotype)]
miss <- setdiff(unique(x$ph), names(names_s)); if (length(miss)) stop("No display name for: ", paste(miss, collapse = ", "))
x[, `:=`(rs = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), re = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end))]
sup <- c("0" = "⁰", "1" = "¹", "2" = "²", "3" = "³", "4" = "⁴", "5" = "⁵", "6" = "⁶", "7" = "⁷", "8" = "⁸", "9" = "⁹", "-" = "⁻")
fmt_p <- function(p, dag = FALSE) {
  ex <- floor(log10(p)); m <- p / 10^ex; if (round(m, 2) >= 10) { m <- m / 10; ex <- ex + 1 }
  e <- paste(sup[strsplit(as.character(ex), "")[[1]]], collapse = "")
  paste0(sprintf("%.2f×10", m), e, if (dag) "†" else "")
}
# GWAS Catalog hits (stage 7)
cf <- file.path(pa$results_root, "07_external", ver, "catalog", "catalog_counts_summary.csv")
cat_tab <- if (file.exists(cf)) fread(cf) else { cat("NOTE: stage 7 catalog counts not found; GWAS Catalog hit left blank\n"); NULL }
cat_hit <- function(p, ch, s, e) {
  if (is.null(cat_tab)) return("")
  y <- cat_tab[phenotype == p & chr == ch & cnvr_start <= e & cnvr_end >= s]
  if (!nrow(y)) return("No (not in catalog table)")
  if (any(y$n_rsid_CNVR > 0)) "Yes" else "No"
}
build <- function(d, single) {
  d <- copy(d)[order(ph, chr, rs)]
  out <- d[, .(Phenotype = names_s[ph], Cytoband = sub("^1q21\\.2$", "1q21.1", locus_cytoband_final), Chr = chr,
               `Sig probe start` = event_start, `Sig probe end` = event_end, `N sig probes` = n_sig_probes_in_event,
               `CNV region start (hg19)` = rs, `CNV region end (hg19)` = re, `Lead probe position` = lead_pos,
               `P-value` = mapply(fmt_p, pval, pass_gw_experiment), beta = round(lead_beta, 3), SE = round(lead_se, 3),
               `N del` = if (single) fifelse(cnv_type == "deletion", as.character(n_del), "-") else as.character(n_del),
               `N dup` = if (single) fifelse(cnv_type == "duplication", as.character(n_dup), "-") else as.character(n_dup),
               `N carriers` = n_carrier, `Carrier freq (%)` = round(carrier_freq_pct, 3), `Total N` = total_n)]
  if (single) out[, `CNV type` := d$cnv_type] else
    out[, `Also significant (DEL/DUP)` := mapply(function(p, ch, s, e) {
      y <- x[cnv_type != "mirror" & ph == p & chr == ch & rs <= e & re >= s]
      if (!nrow(y)) "" else paste(sprintf("%s (p=%s)", y$cnv_type, sapply(y$pval, fmt_p)), collapse = "; ")
    }, d$ph, d$chr, d$rs, d$re)]
  out[, `GWAS Catalog hit` := mapply(cat_hit, d$ph, d$chr, d$rs, d$re)]
  out
}
s5 <- build(x[cnv_type %in% c("deletion", "duplication")], TRUE)
s6 <- build(x[cnv_type == "mirror"], FALSE)
cat("\n== Table S5 (single-type models):", nrow(s5), "rows\n"); print(s5)
cat("\n== Table S6 (mirror model):", nrow(s6), "rows\n"); print(s6)
fwrite(s5, file.path(od, "table_s5_single_type.csv"), bom = TRUE); fwrite(s6, file.path(od, "table_s6_mirror.csv"), bom = TRUE)
cap <- thr[, .(model, Meff, p_cnv_wide = signif(p_cnv_wide, 3), p_experiment_wide = signif(p_experiment_wide, 3))]
cat("\n== Thresholds for the captions (all participants)\n"); print(cap)
fwrite(cap, file.path(od, "table_s5_s6_thresholds.csv"))
cat("\nBoundary check (CNV region used vs the stage 5 CNVR and the call extent):\n")
print(x[, .(phenotype = ph, cnv_type, chr, used = sprintf("%d-%d", rs, re), stage5_cnvr = sprintf("%d-%d", event_start, event_end),
            call_extent = sprintf("%s-%s", cnv_region_start, cnv_region_end))][order(chr, phenotype)])
cat("\nSaved: table_s1_exclusion_regions.csv, table_s5_single_type.csv, table_s6_mirror.csv, table_s5_s6_thresholds.csv in", od, "\n")
