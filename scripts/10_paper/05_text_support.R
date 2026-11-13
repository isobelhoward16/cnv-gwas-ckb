#!/usr/bin/env Rscript
# Stage 10: numbers quoted in the Results text that are not in a table (aggregate only):
#  (1) text_support_loci.csv: one row per stage 5 association: locus label (as in Table 1), lead probe and position,
#      CNVR, genes (lead probe, significant probes, event window, CNV region), and the untransformed trait mean,
#      median and SD in deletion carriers, duplication carriers and non-carriers (from stage 5);
#  (2) text_support_lookups.csv: glucose, SBP and DBP at the 16p11.2 BP2-BP3 lead probe, all three models, from
#      slurm/10_text_support_lookups.sbatch (same code and model as the main analysis, all participants, this probe
#      only, no minimum-carrier filter).
# Run from the repository root:  Rscript scripts/10_paper/05_text_support.R
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
od <- file.path(pa$results_root, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x[, `:=`(ph = sub("_rint$", "", phenotype), s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start),
         e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end))]
x[, label := sub("^1q21\\.2$", "1q21.1", locus_cytoband_final)]
x[label == "16p11.2", label := fifelse(s < 29.3e6, "16p11.2a (BP2-BP3)", "16p11.2b (BP4-BP5)")]
keep <- intersect(c("label", "ph", "cnv_type", "lead_probe", "chr", "lead_pos", "s", "e", "event_start", "event_end", "cnv_region_start",
                    "cnv_region_end", "pval", "lead_beta", "lead_se", "n_del", "n_dup", "total_n", "carrier_freq_pct", "lead_genes",
                    "event_span_genes", "auwerx_genes", "cnv_region_genes", "mean_all", "median_all", "mean_del", "median_del", "sd_del",
                    "mean_dup", "median_dup", "sd_dup", "mean_noncarrier", "median_noncarrier", "sd_noncarrier"), names(x))
loci <- x[, ..keep]; setnames(loci, c("s", "e"), c("cnvr_start", "cnvr_end")); setorder(loci, chr, cnvr_start, ph, cnv_type)
fwrite(loci, file.path(od, "text_support_loci.csv"))
cat("text_support_loci.csv:", nrow(loci), "rows (expected 23)\n"); print(loci[, .(label, ph, cnv_type, lead_pos, lead_genes, mean_del, mean_dup, mean_noncarrier)])

gdir <- file.path(pa$results_root, "04_gwas", ver, "all", "sensitivity")
gfile <- function(model, ph) { f <- list.files(gdir, pattern = paste0("^gwas_linear_unpruned_combined_", model, "_", ph, "(_rint)?(_mincar[0-9]+)?_.*_lead16a\\.csv$"), full.names = TRUE)
  if (!length(f)) return(NA_character_); f[which.max(file.mtime(f))] }
lead_16a <- x[label == "16p11.2a (BP2-BP3)" & ph == "bmi_calc" & cnv_type == "mirror", lead_probe]
if (length(lead_16a) != 1) stop("16p11.2 BP2-BP3 BMI mirror lead not found")
if (lead_16a != "AX-31870881") stop("16p11.2 BP2-BP3 lead is ", lead_16a, ", but the lookup GWAS was run for AX-31870881")
look <- CJ(probe = lead_16a, ph = c("random_glucose", "sbp_mean", "dbp_mean"), model = c("deletion", "duplication", "mirror"))
lk <- rbindlist(lapply(seq_len(nrow(look)), function(i) { f <- gfile(look$model[i], look$ph[i])
  if (is.na(f)) return(data.table(probe = look$probe[i], ph = look$ph[i], model = look$model[i], note = "no GWAS file"))
  g <- fread(f, select = c("probe_id", "beta", "se", "pval", "n_carriers", "n_analyzed"))[probe_id == look$probe[i]]
  if (!nrow(g)) return(data.table(probe = look$probe[i], ph = look$ph[i], model = look$model[i], note = "probe not in the output"))
  data.table(probe = look$probe[i], ph = look$ph[i], model = look$model[i], beta = g$beta, se = g$se, pval = g$pval,
             n_carriers = g$n_carriers, n_analyzed = g$n_analyzed, note = "", file = basename(f)) }), fill = TRUE)
fwrite(lk, file.path(od, "text_support_lookups.csv")); cat("\ntext_support_lookups.csv (16p11.2 BP2-BP3 lead", lead_16a, ")\n"); print(lk[, !"file"])
