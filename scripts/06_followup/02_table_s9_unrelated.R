#!/usr/bin/env Rscript
# Stage 6, step 2: each main signal at its lead probe in all participants and in unrelated participants.
# All-participant results come from the main association analysis; unrelated results from the lead-probe
# run (same model, unrelated participants, no carrier filter). Significance in unrelated participants
# uses the CNV-wide threshold re-estimated in that subset (0.05 / maximum unrelated Meff).
# Output: results/06_followup/<version>/unrelated_leads/table_s9_unrelated_sensitivity.csv
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version; mc <- cfg$pipeline$analysis$min_carriers
options(width = 250)
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
thr <- fread(file.path(pa$results_root, "03_meff", ver, "meff_thresholds.csv"))[relatedness == "unrelated", .(model, p_thr_unrel = p_cnv_wide)]
g <- function(set, m, ph, tag = "") file.path(pa$results_root, "04_gwas", ver, set, if (nzchar(tag)) "sensitivity" else "",
  sprintf("gwas_linear_unpruned_combined_%s_%s_mincar%s_%s%s.csv", m, ph, mc, set, if (nzchar(tag)) paste0("_", tag) else ""))
cols <- c("probe_id", "beta", "se", "pval", "n_carriers", "n_analyzed")
res <- rbindlist(lapply(seq_len(nrow(x)), function(i) {
  m <- x$cnv_type[i]; ph <- x$phenotype[i]; pr <- x$lead_probe[i]
  a <- fread(g("all", m, ph), select = cols)[probe_id == pr]
  fu <- g("unrelated", m, ph, "leads"); u <- if (file.exists(fu)) fread(fu, select = cols)[probe_id == pr] else data.table()
  gu <- fread(g("unrelated", m, ph), select = c("probe_id", "pval"))
  rank <- if (nrow(u) && pr %in% gu$probe_id) round(100 * mean(gu$pval <= u$pval, na.rm = TRUE), 2) else NA_real_
  data.table(cytoband = x$locus_cytoband_final[i], phenotype = sub("_rint$", "", ph), model = m, lead_probe = pr,
             beta_all = a$beta, se_all = a$se, p_all = a$pval, carriers_all = a$n_carriers, N_all = a$n_analyzed,
             beta_unrel = if (nrow(u)) u$beta else NA_real_, se_unrel = if (nrow(u)) u$se else NA_real_,
             p_unrel = if (nrow(u)) u$pval else NA_real_, carriers_unrel = if (nrow(u)) u$n_carriers else NA_integer_,
             N_unrel = if (nrow(u)) u$n_analyzed else NA_integer_, rank_pct_unrel = rank)
}))
res <- merge(res, thr, by = "model")
res[, `:=`(ratio = beta_unrel / beta_all, concordant = sign(beta_unrel) == sign(beta_all), pass_unrel = !is.na(p_unrel) & p_unrel <= p_thr_unrel)]
setorder(res, cytoband, phenotype, model)
cat("== Table S9: unrelated sensitivity at the", nrow(res), "lead probes (all-participant signals)\n")
print(res[, .(cytoband, phenotype, model, lead_probe, beta_all = round(beta_all, 3), p_all = signif(p_all, 3), carriers_all,
              beta_unrel = round(beta_unrel, 3), se_unrel = round(se_unrel, 3), p_unrel = signif(p_unrel, 3), carriers_unrel, N_unrel,
              ratio = round(ratio, 2), concordant, p_thr_unrel = signif(p_thr_unrel, 3), pass_unrel, rank_pct_unrel)])
cat("\n== Summary for the text\n")
cat("Unrelated N range (at the leads):", min(res$N_unrel, na.rm = TRUE), "-", max(res$N_unrel, na.rm = TRUE), "\n")
cat("Leads with an unrelated estimate:", sum(!is.na(res$beta_unrel)), "of", nrow(res), "\n")
cat("Directionally concordant:", sum(res$concordant, na.rm = TRUE), "of", nrow(res), "\n")
cat("Ratio unrelated/full beta: median", round(median(res$ratio, na.rm = TRUE), 2), "| range", round(min(res$ratio, na.rm = TRUE), 2), "-", round(max(res$ratio, na.rm = TRUE), 2), "\n")
cat("CNV-wide significant in unrelated:", sum(res$pass_unrel), "of", nrow(res), "\n"); print(res[pass_unrel == TRUE, .(cytoband, phenotype, model, p_unrel = signif(p_unrel, 3))])
cat("Lead in top 1% / top 5% of the unrelated GWAS:", sum(res$rank_pct_unrel <= 1, na.rm = TRUE), "/", sum(res$rank_pct_unrel <= 5, na.rm = TRUE), "of", sum(!is.na(res$rank_pct_unrel)), "ranked\n")
od <- file.path(pa$results_root, "06_followup", ver, "unrelated_leads"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
fwrite(res, file.path(od, "table_s9_unrelated_sensitivity.csv")); cat("\nSaved:", file.path(od, "table_s9_unrelated_sensitivity.csv"), "\n")
