#!/bin/bash
#SBATCH --job-name=rpg_fastadj
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 6: random plasma glucose adjusted for time since last meal (hours_since_last_ate_x10).
# Same GWAS code and sample as the main analysis (scripts/04_gwas/01_gwas_task.R, all participants),
# tested at every probe with main-analysis P < 1e-4 for RPG in that model (covers every RPG signal
# and near-miss).
# Submit from the repository root:  sbatch slurm/06_fasting_rpg.sbatch
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
OUT=results/04_gwas/$VER/all/sensitivity; mkdir -p "$OUT/probe_lists"
echo "Started $(date)"
Rscript -e '
source("R/config.R"); cfg <- load_config(); pa <- cfg$paths; v <- cfg$pipeline$callset_version; mc <- cfg$pipeline$analysis$min_carriers
for (m in c("deletion", "duplication", "mirror")) {
  g <- fread(file.path(pa$results_root, "04_gwas", v, "all", sprintf("gwas_linear_unpruned_combined_%s_random_glucose_rint_mincar%s_all.csv", m, mc)), select = c("probe_id", "pval"))
  p <- g[!is.na(pval) & pval < 1e-4, probe_id]
  writeLines(p, file.path(pa$results_root, "04_gwas", v, "all", "sensitivity", "probe_lists", paste0("rpg_", m, ".txt")))
  cat(m, ": probes with main P < 1e-4:", length(p), "\n") }'
for M in deletion duplication mirror; do
  echo "== $M $(date)"
  EXTRA_COVARS=hours_since_last_ate_x10 PROBE_LIST="$OUT/probe_lists/rpg_${M}.txt" OUT_TAG=fastadj \
    Rscript scripts/04_gwas/01_gwas_task.R "$M" random_glucose_rint all
done
Rscript -e '
source("R/config.R"); cfg <- load_config(); pa <- cfg$paths; v <- cfg$pipeline$callset_version; mc <- cfg$pipeline$analysis$min_carriers
thr <- fread(file.path(pa$results_root, "03_meff", v, "meff_thresholds.csv"))[relatedness == "all"]
d <- file.path(pa$results_root, "04_gwas", v, "all")
for (m in c("deletion", "duplication", "mirror")) {
  a <- fread(file.path(d, sprintf("gwas_linear_unpruned_combined_%s_random_glucose_rint_mincar%s_all.csv", m, mc)), select = c("probe_id", "beta", "pval", "n_analyzed"))
  f <- file.path(d, "sensitivity", sprintf("gwas_linear_unpruned_combined_%s_random_glucose_rint_mincar%s_all_fastadj.csv", m, mc))
  if (!file.exists(f)) { cat(m, ": no fasting-adjusted file\n"); next }
  b <- fread(f, select = c("probe_id", "beta", "pval", "n_analyzed"))
  x <- merge(a, b, by = "probe_id", suffixes = c("_main", "_fastadj"))[order(pval_main)]
  t <- thr[model == m, p_cnv_wide]
  cat("\n==", m, "| CNV-wide threshold", signif(t, 3), "| probes", nrow(x), "| significant main:", sum(x$pval_main < t),
      "| still significant after adjustment:", sum(x$pval_main < t & x$pval_fastadj < t), "\n")
  print(head(x[, .(probe_id, beta_main = round(beta_main, 3), P_main = signif(pval_main, 3), beta_fastadj = round(beta_fastadj, 3),
                   P_fastadj = signif(pval_fastadj, 3), N_main = n_analyzed_main, N_fastadj = n_analyzed_fastadj)], 15)) }'
echo "Finished $(date)"
