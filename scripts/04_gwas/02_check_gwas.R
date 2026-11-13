#!/usr/bin/env Rscript
# Stage 4 check: every model x phenotype x sample-set result exists and is complete, the number of
# probes tested equals the Meff probe count for the same task, and the number of probes passing the
# CNV-wide threshold.
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; an <- cfg$pipeline$analysis; ver <- cfg$pipeline$callset_version
g <- file.path(pa$results_root, "04_gwas", ver)
tk <- CJ(model = an$models, phenotype = an$phenotypes, relatedness = an$relatedness)
tk[, file := file.path(g, relatedness, sprintf("gwas_linear_unpruned_combined_%s_%s_mincar%s_%s.csv", model, phenotype, an$min_carriers, relatedness))]
tk[, exists := file.exists(file)]
cat("GWAS results expected:", nrow(tk), "| found:", sum(tk$exists), "\n")
if (any(!tk$exists)) { cat("Missing (rerun these array tasks):\n"); print(tk[exists == FALSE, .(model, phenotype, relatedness, task = match(paste(model, phenotype, relatedness), system("Rscript R/tasks.R", intern = TRUE)))]) }
res <- tk[exists == TRUE, {
  r <- fread(file)
  .(n_probes = nrow(r), n_errors = sum(!is.na(r$error) & r$error != ""), min_p = min(r$pval, na.rm = TRUE),
    n_samples = max(r$n_analyzed, na.rm = TRUE))
}, by = .(model, phenotype, relatedness)]
mt <- file.path(pa$results_root, "03_meff", ver, "meff_all_tasks.csv")
if (file.exists(mt)) {
  m <- fread(mt)[, .(model = CNV_Type, phenotype = Phenotype, relatedness = Relatedness, meff_probes = N_Probes)]
  res <- merge(res, m, by = c("model", "phenotype", "relatedness"), all.x = TRUE)
  cat("Probe count equals Meff probe count:", sum(res$n_probes == res$meff_probes, na.rm = TRUE), "of", nrow(res), "\n")
  th <- fread(file.path(pa$results_root, "03_meff", ver, "meff_thresholds.csv"))
  res <- merge(res, th[, .(model, relatedness, p_cnv_wide)], by = c("model", "relatedness"), all.x = TRUE)
  res[, n_sig := tk[exists == TRUE][.SD, on = .(model, phenotype, relatedness), mapply(function(f, p) sum(fread(f, select = "pval")$pval <= p, na.rm = TRUE), x.file, i.p_cnv_wide)]]
}
cat("Tasks with errors:", sum(res$n_errors > 0), "\n\n")
print(res[order(relatedness, model, phenotype)])
fwrite(res, file.path(g, "gwas_check.csv"))
