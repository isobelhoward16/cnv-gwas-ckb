#!/usr/bin/env Rscript
# Stage 6, step 1: lead probes of the main signals (all participants), one list per model x trait,
# for the analysis in unrelated participants.
# Output: results/06_followup/<version>/unrelated_leads/lists/*.txt and tasks.txt
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
d <- file.path(pa$results_root, "06_followup", ver, "unrelated_leads"); dir.create(file.path(d, "lists"), recursive = TRUE, showWarnings = FALSE)
tasks <- unique(x[, .(model = cnv_type, phenotype)])
for (i in seq_len(nrow(tasks))) {
  f <- file.path(d, "lists", sprintf("leads_%s_%s.txt", tasks$model[i], tasks$phenotype[i]))
  writeLines(unique(x[cnv_type == tasks$model[i] & phenotype == tasks$phenotype[i], lead_probe]), f)
  tasks[i, list_file := f]
}
fwrite(tasks, file.path(d, "tasks.txt"), sep = " ", col.names = FALSE)
cat("Signals:", nrow(x), "| model x trait tasks:", nrow(tasks), "\nWritten:", file.path(d, "tasks.txt"), "\n")
