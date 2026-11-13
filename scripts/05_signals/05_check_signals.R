#!/usr/bin/env Rscript
# Stage 5, step 4: consistency checks across Meff (stage 3), GWAS (stage 4) and signals (stage 5).
# Every check prints PASS or FAIL with counts.
#   1. Population: the regression sample for each trait (n_analyzed + n_excluded, constant over probes)
#      equals the Meff sample (N_Samples) for the same trait, model and sample set; and is the same
#      across the deletion, duplication and mirror models.
#   2. Probes: probes tested in each GWAS = probes counted by Meff.
#   3. Threshold: one CNV-wide threshold per model and sample set (0.05 / max Meff); the significant-probe
#      file holds exactly the GWAS probes with P below it.
#   4. CNVR membership: every significant probe sits in exactly one CNVR of its trait and model.
#   5. Leads: each CNVR lead has the smallest P among its members, and its P and beta equal the GWAS.
#   6. Carriers and N: carriers at the lead equal the GWAS carrier count at that probe; total N equals
#      the GWAS n_analyzed at that probe.
#   7. Flags: experiment-wide flag = P <= 0.05 / (Meff x n_effective_traits).
# Then prints the signal table (trait, model, cytoband, lead P, carriers, experiment-wide).
# Usage: Rscript scripts/05_signals/05_check_signals.R [all|unrelated]
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version; an <- cfg$pipeline$analysis
args <- commandArgs(trailingOnly = TRUE)
rel <- if (length(args)) args[1] else "all"
stopifnot(rel %in% c("all", "unrelated"))
rg <- if (rel == "all") "related" else "unrelated"
mc <- an$min_carriers; ntr <- an$n_effective_traits
res_all <- logical(0)
report <- function(name, ok, detail = "") {
  cat(sprintf("[%s] %s%s\n", if (isTRUE(ok)) "PASS" else "FAIL", name, if (nchar(detail)) paste0(" | ", detail) else ""))
  res_all <<- c(res_all, setNames(isTRUE(ok), name))
}
cat("== Stage 5 checks |", rel, "| callset", ver, "\n\n")

# ---- inputs ----
meff <- rbindlist(lapply(list.files(file.path(pa$results_root, "03_meff", ver, "tasks"), "^meff_summary_.*\\.csv$", full.names = TRUE), fread))
meff <- meff[Relatedness == rel, .(model = CNV_Type, phenotype = Phenotype, N_meff = N_Samples, probes_meff = N_Probes, Meff)]
thr <- fread(file.path(pa$results_root, "03_meff", ver, "meff_thresholds.csv"))[relatedness == rel]
gdir <- file.path(pa$results_root, "04_gwas", ver, rel)
gfile <- function(ct, ph) file.path(gdir, sprintf("gwas_linear_unpruned_combined_%s_%s_mincar%s_%s.csv", ct, ph, mc, rel))
gw <- list()
for (i in seq_len(nrow(meff))) {
  f <- gfile(meff$model[i], meff$phenotype[i])
  if (!file.exists(f)) stop("Missing GWAS file: ", f)
  gw[[paste(meff$model[i], meff$phenotype[i])]] <- fread(f, select = c("probe_id", "beta", "pval", "n_carriers", "n_analyzed", "n_excluded"))
}
sdir <- file.path(pa$results_root, "05_signals", ver, rel)
sig  <- fread(file.path(sdir, paste0("significant_probes_unpruned_", rg, ".csv")))
cnvr <- fread(file.path(sdir, "cnvr_outputs", paste0("cnvr_by_model_supplemental_unpruned_", rg, ".csv")))
memb <- fread(file.path(sdir, "cnvr_outputs", paste0("cnvr_probe_membership_unpruned_", rg, ".csv")))

# ---- 1. population ----
pop <- rbindlist(lapply(seq_len(nrow(meff)), function(i) {
  g <- gw[[paste(meff$model[i], meff$phenotype[i])]]
  tot <- unique(g[!is.na(n_analyzed), n_analyzed + n_excluded])
  data.table(model = meff$model[i], phenotype = meff$phenotype[i], N_meff = meff$N_meff[i],
             N_regression = if (length(tot) == 1) tot else NA_integer_, n_distinct_totals = length(tot))
}))
report("1a. regression sample is the same for every probe of a GWAS", all(pop$n_distinct_totals == 1),
       sprintf("%d of %d GWAS", sum(pop$n_distinct_totals == 1), nrow(pop)))
report("1b. regression N = Meff N (same trait, model, sample set)", all(pop$N_regression == pop$N_meff, na.rm = FALSE),
       sprintf("%d of %d match", sum(pop$N_regression == pop$N_meff, na.rm = TRUE), nrow(pop)))
wide <- dcast(pop, phenotype ~ model, value.var = "N_regression")
same <- wide[, apply(.SD, 1, function(x) length(unique(x)) == 1), .SDcols = setdiff(names(wide), "phenotype")]
report("1c. N is the same in the deletion, duplication and mirror models", all(same), sprintf("%d of %d traits", sum(same), length(same)))
cat("    N per trait (", rel, "):\n", sep = ""); print(wide)
if (any(pop$N_regression != pop$N_meff, na.rm = TRUE) || any(is.na(pop$N_regression))) print(pop[N_regression != N_meff | is.na(N_regression)])

# ---- 2. probes ----
pr <- rbindlist(lapply(seq_len(nrow(meff)), function(i) data.table(model = meff$model[i], phenotype = meff$phenotype[i],
        probes_meff = meff$probes_meff[i], probes_gwas = nrow(gw[[paste(meff$model[i], meff$phenotype[i])]]))))
report("2. probes tested = probes in Meff", all(pr$probes_meff == pr$probes_gwas), sprintf("%d of %d", sum(pr$probes_meff == pr$probes_gwas), nrow(pr)))

# ---- 3. threshold and significant probes ----
report("3a. one threshold per model (0.05 / max Meff)", nrow(thr) == 3 && all(abs(thr$p_cnv_wide - 0.05 / thr$Meff) < 1e-12),
       paste(sprintf("%s %.3g", thr$model, thr$p_cnv_wide), collapse = "; "))
mx <- meff[, .(Meff_max = max(Meff)), by = model][thr, on = "model"]
report("3b. threshold Meff is the maximum over the traits", all(mx$Meff_max == mx$Meff))
exp_sig <- rbindlist(lapply(seq_len(nrow(meff)), function(i) {
  t <- thr[model == meff$model[i], p_cnv_wide]
  g <- gw[[paste(meff$model[i], meff$phenotype[i])]]
  ids <- g[!is.na(pval) & pval < t, probe_id]
  if (!length(ids)) return(NULL)
  data.table(model = meff$model[i], phenotype = meff$phenotype[i], probe_id = ids)
}))
got_sig <- unique(sig[, .(model = cnv_type, phenotype, probe_id)])
report("3c. significant-probe file = GWAS probes with P below the threshold",
       nrow(fsetdiff(exp_sig, got_sig)) == 0 && nrow(fsetdiff(got_sig, exp_sig)) == 0,
       sprintf("expected %d, file %d, missing %d, extra %d", nrow(exp_sig), nrow(got_sig),
               nrow(fsetdiff(exp_sig, got_sig)), nrow(fsetdiff(got_sig, exp_sig))))
cat("    Significant probes per trait and model:\n"); print(dcast(exp_sig[, .N, by = .(phenotype, model)], phenotype ~ model, value.var = "N", fill = 0))

# ---- 4. membership ----
mm <- memb[, .(n_cnvr = uniqueN(cnvr_id)), by = .(model = cnv_type, phenotype, probe_id)]
report("4a. every significant probe is in a CNVR", nrow(fsetdiff(exp_sig, mm[, .(model, phenotype, probe_id)])) == 0,
       sprintf("%d not assigned", nrow(fsetdiff(exp_sig, mm[, .(model, phenotype, probe_id)]))))
report("4b. no probe in two CNVRs of the same trait and model", all(mm$n_cnvr == 1), sprintf("%d in more than one", sum(mm$n_cnvr > 1)))
report("4c. no CNVR member is non-significant", nrow(fsetdiff(mm[, .(model, phenotype, probe_id)], exp_sig)) == 0)

# ---- 5. leads ----
mlead <- memb[, .(min_member_p = min(probe_pval)), by = .(cnv_type, phenotype, cnvr_id)]
cv <- merge(cnvr, mlead, by = c("cnv_type", "phenotype", "cnvr_id"), all.x = TRUE)
cv[, `:=`(gwas_p = NA_real_, gwas_beta = NA_real_, gwas_ncar = NA_integer_, gwas_n = NA_integer_)]
for (i in seq_len(nrow(cv))) {
  g <- gw[[paste(cv$cnv_type[i], cv$phenotype[i])]][probe_id == cv$lead_probe[i]]
  if (nrow(g)) set(cv, i, c("gwas_p", "gwas_beta", "gwas_ncar", "gwas_n"), list(g$pval, g$beta, g$n_carriers, g$n_analyzed))
}
report("5a. lead = smallest P among its CNVR members", all(abs(cv$pval - cv$min_member_p) <= 1e-12 * pmax(cv$pval, 1e-300), na.rm = FALSE),
       sprintf("%d of %d CNVRs", sum(abs(cv$pval - cv$min_member_p) <= 1e-12 * pmax(cv$pval, 1e-300), na.rm = TRUE), nrow(cv)))
report("5b. lead P and beta equal the GWAS", all(abs(cv$pval / cv$gwas_p - 1) < 1e-6 & abs(cv$lead_beta - cv$gwas_beta) < 1e-6, na.rm = FALSE))

# ---- 6. carriers and N at the lead ----
cv[, ncar := fifelse(cnv_type == "deletion", n_del, fifelse(cnv_type == "duplication", n_dup, n_del + n_dup))]
report("6a. carriers at the lead = GWAS carrier count", all(cv$ncar == cv$gwas_ncar), sprintf("%d of %d", sum(cv$ncar == cv$gwas_ncar, na.rm = TRUE), nrow(cv)))
report("6b. total N at the lead = GWAS n_analyzed", all(cv$total_n == cv$gwas_n), sprintf("%d of %d", sum(cv$total_n == cv$gwas_n, na.rm = TRUE), nrow(cv)))
report("6c. at least the minimum carriers at every lead", all(cv$ncar >= mc))
if (any(cv$ncar != cv$gwas_ncar | cv$total_n != cv$gwas_n, na.rm = TRUE))
  print(cv[ncar != gwas_ncar | total_n != gwas_n, .(phenotype, cnv_type, chr, lead_probe, ncar, gwas_ncar, total_n, gwas_n)])

# ---- 7. flags ----
cv[thr, on = .(cnv_type = model), `:=`(Meff = i.Meff, p_cnv = i.p_cnv_wide)]
report("7a. every CNVR lead passes the CNV-wide threshold", all(cv$pval <= cv$p_cnv))
report("7b. experiment-wide flag = P <= 0.05 / (Meff x n_effective_traits)", all(cv$pass_gw_experiment == (cv$pval <= 0.05 / (cv$Meff * ntr))),
       sprintf("n_effective_traits = %s; %d CNVRs pass", ntr, sum(cv$pval <= 0.05 / (cv$Meff * ntr))))

# ---- signal table ----
cat("\n== Signals (", rel, "): ", nrow(cv), " CNVRs (trait x model)\n", sep = "")
setorder(cv, chr, event_start, phenotype, cnv_type)
print(cv[, .(cytoband = locus_cytoband_final, chr, start = event_start, end = event_end, phenotype = sub("_rint$", "", phenotype),
             model = cnv_type, lead_probe, P = signif(pval, 3), beta = round(lead_beta, 3), n_del, n_dup, N = total_n,
             exp_wide = pass_gw_experiment)])
cat(sprintf("\n== %d of %d checks passed\n", sum(res_all), length(res_all)))
if (!all(res_all)) cat("FAILED:", paste(names(res_all)[!res_all], collapse = "; "), "\n")
