#!/usr/bin/env Rscript
# Stage 3 summary: combines the Meff tasks, takes the maximum Meff across phenotypes for each model
# and sample set, and writes the CNV-wide (0.05 / Meff) and experiment-wide (0.05 / (Meff x 11))
# thresholds used by the later stages.
# Output: results/03_meff/<version>/meff_all_tasks.csv, meff_thresholds.csv
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; an <- cfg$pipeline$analysis; ver <- cfg$pipeline$callset_version
d <- file.path(pa$results_root, "03_meff", cfg$pipeline$callset_version)
fs <- list.files(file.path(d, "tasks"), pattern = "^meff_summary_.*\\.csv$", full.names = TRUE)
x <- rbindlist(lapply(fs, fread))
expected <- CJ(CNV_Type = an$models, Phenotype = an$phenotypes, Relatedness = an$relatedness)
miss <- expected[!x, on = .(CNV_Type, Phenotype, Relatedness)]
cat("Tasks expected:", nrow(expected), "| found:", nrow(x), "| missing:", nrow(miss), "\n")
if (nrow(miss)) { print(miss); stop("Some Meff tasks are missing: rerun them before using the thresholds") }
fwrite(x[order(CNV_Type, Relatedness, Phenotype)], file.path(d, "meff_all_tasks.csv"))

thr <- x[, .(Meff = max(Meff), phenotype_at_max = Phenotype[which.max(Meff)],
             N_Probes_at_max = N_Probes[which.max(Meff)], N_Samples_at_max = N_Samples[which.max(Meff)]),
         by = .(model = CNV_Type, relatedness = Relatedness)]
thr[, `:=`(p_cnv_wide = 0.05 / Meff, p_experiment_wide = 0.05 / (Meff * an$n_effective_traits))]
setorder(thr, relatedness, model)
fwrite(thr, file.path(d, "meff_thresholds.csv"))
cat("\nThresholds (maximum Meff across the", length(an$phenotypes), "phenotypes):\n")
print(thr[, .(model, relatedness, Meff, N_Probes_at_max, N_Samples_at_max,
              p_cnv_wide = signif(p_cnv_wide, 3), p_experiment_wide = signif(p_experiment_wide, 3))])

cat("\nPer phenotype, all participants (N_Probes should equal the probes tested in the GWAS):\n")
print(dcast(x[Relatedness == "all"], Phenotype ~ CNV_Type, value.var = c("N_Probes", "Meff")))

