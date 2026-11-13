#!/usr/bin/env Rscript
# Stage 10: Table 1, CNV-wide significant locus-trait associations (all participants), from the stage 5 signal
# table and the stage 6 unrelated table (Table S9). Columns:
#   Locus | Trait | CNVR (Mb) | Model(s) | Carriers (del/dup) | beta (95% CI) | P | Exp.-wide | Unrelated
# Rules:
#   - one row per locus x trait; the model shown (CNVR, carriers, beta, P) is the mirror model when the mirror model is
#     significant for that locus-trait, otherwise the model with the smallest P;
#   - CNVR = event boundaries where available, otherwise the stage 5 CNVR (Mb, 2 decimals);
#   - carriers (del/dup) = deletion and duplication carriers at the lead probe of the model shown, in that model's
#     sample; for single-type models the type not tested comes from 00_table1_carriers_both_types.R;
#   - locus label = stage 5 cytoband, with "1q21.2" read as 1q21.1 (the sitting-height CNVR is the 1q21.1 distal CNV);
#     when one cytoband holds CNVRs that do not overlap, they are lettered a, b, ... by position;
#   - Exp.-wide and Unrelated: "Yes" if every significant model passes, "No" if none does, otherwise "Yes (<models>)".
# Aggregate output only. Output: results/10_paper/<version>/table1.csv (display text) and table1_numeric.csv
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
od <- file.path(pa$results_root, "10_paper", ver); dir.create(od, recursive = TRUE, showWarnings = FALSE)
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
s9f <- file.path(pa$results_root, "06_followup", ver, "unrelated_leads", "table_s9_unrelated_sensitivity.csv")
s9 <- if (file.exists(s9f)) fread(s9f) else { cat("NOTE: Table S9 not found yet; Unrelated column left blank\n"); NULL }
cbf <- file.path(od, "table1_carriers_both_types.csv")
if (!file.exists(cbf)) stop("Run scripts/10_paper/00_table1_carriers_both_types.R first (", cbf, " not found)")
cb <- fread(cbf)
x <- merge(x, cb[, .(phenotype, cnv_type, lead_probe, n_del_at_lead, n_dup_at_lead)], by = c("phenotype", "cnv_type", "lead_probe"), all.x = TRUE)
if (x[, anyNA(n_del_at_lead) | anyNA(n_dup_at_lead)]) stop("Both-type carrier counts missing for some lead probes")

trait_names <- c(bmi_calc = "BMI", bmi_at_25 = "BMI at age 25", standing_height_mm = "Standing height", sitting_height_mm = "Sitting height",
                 weight_kg = "Weight", waist_mm = "Waist circumference", hip_mm = "Hip circumference", waist_hip_ratio = "Waist-hip ratio",
                 fat_percent = "Fat percentage", sbp_mean = "SBP", dbp_mean = "DBP", heart_rate_mean = "Resting heart rate", random_glucose = "RPG")
model_names <- c(deletion = "Del", duplication = "Dup", mirror = "Mirror")
x[, `:=`(ph = sub("_rint$", "", phenotype), s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start),
         e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end))]
miss <- setdiff(unique(x$ph), names(trait_names)); if (length(miss)) stop("No display name for trait(s): ", paste(miss, collapse = ", "))
x[, band := sub("^1q21\\.2$", "1q21.1", locus_cytoband_final)]
# letter non-overlapping CNVR clusters within a cytoband
x[, cl := {
  o <- order(s); grp <- integer(.N); g <- 0L; hi <- -Inf
  for (i in o) { if (s[i] > hi) { g <- g + 1L; hi <- e[i] } else hi <- max(hi, e[i]); grp[i] <- g }
  grp }, by = band]
x[, label := if (uniqueN(cl) > 1) paste0(band, letters[cl]) else band, by = band]
lab_chk <- unique(x[, .(locus_cytoband_final, label, cluster_Mb = sprintf("%.2f-%.2f", s / 1e6, e / 1e6))])
cat("== Locus labels\n"); print(lab_chk[order(label)])

if (!is.null(s9)) s9[, phenotype := sub("_rint$", "", phenotype)]
yn <- function(models, pass) {
  if (all(pass)) "Yes" else if (!any(pass)) "No" else paste0("Yes (", paste(model_names[models[pass]], collapse = ", "), ")")
}
fmt_p <- function(p) { ex <- floor(log10(p)); m <- p / 10^ex; if (round(m, 2) >= 10) { m <- m / 10; ex <- ex + 1 }
  sprintf("%.2f × 10%s", m, gsub("-", "−", as.character(ex))) }
sgn <- function(v, plus = TRUE) { t <- sprintf("%.2f", abs(v)); ifelse(v < 0, paste0("−", t), if (plus) paste0("+", t) else t) }
sgn_ci <- function(v) { t <- sprintf("%.2f", abs(v)); ifelse(v < 0, paste0("−", t), t) }
tab <- x[, {
  shown <- if ("mirror" %in% cnv_type) which(cnv_type == "mirror")[1] else which.min(pval)
  lo <- lead_beta[shown] - 1.96 * lead_se[shown]; hi <- lead_beta[shown] + 1.96 * lead_se[shown]
  ord <- order(match(cnv_type, names(model_names)))
  cur_ph <- ph[1]; cur_key <- paste(cnv_type, lead_probe)   # match Table S9 by model and lead probe (16p11.2a and b share a cytoband)
  un <- if (is.null(s9)) "" else { u <- s9[phenotype == cur_ph & paste(model, lead_probe) %in% cur_key]
    if (nrow(u) != length(cur_key)) stop("Table S9 rows not found for ", .BY$Locus, " ", cur_ph)
    yn(u$model, u$pass_unrel) }
  .(chr = chr[1], pos = s[shown], Trait = trait_names[ph[1]], `CNVR (Mb)` = sprintf("%.2f–%.2f", s[shown] / 1e6, e[shown] / 1e6),
    `Model(s)` = paste(model_names[cnv_type[ord]], collapse = ", "), `Carriers (del/dup)` = sprintf("%d/%d", n_del_at_lead[shown], n_dup_at_lead[shown]),
    `beta (95% CI)` = sprintf("%s (%s, %s)", sgn(lead_beta[shown]), sgn_ci(lo), sgn_ci(hi)), P = fmt_p(pval[shown]),
    `Exp.-wide` = yn(cnv_type, pass_gw_experiment), Unrelated = un,
    model_shown = cnv_type[shown], lead_probe = lead_probe[shown], lead_pos = lead_pos[shown], beta = lead_beta[shown], se = lead_se[shown],
    p = pval[shown], n_del = n_del[shown], n_dup = n_dup[shown], n_del_at_lead = n_del_at_lead[shown], n_dup_at_lead = n_dup_at_lead[shown],
    N = total_n[shown])
}, by = .(Locus = label, ph)]
trait_order <- c("BMI", "BMI at age 25", "Fat percentage", "Hip circumference", "Waist circumference", "Waist-hip ratio", "Weight",
                 "Standing height", "Sitting height", "SBP", "DBP", "Resting heart rate", "RPG")   # as the manuscript table
tab[, tord := match(Trait, trait_order)]
setorder(tab, chr, Locus, tord); tab[, tord := NULL]
disp <- tab[, .(Locus, Trait, `CNVR (Mb)`, `Model(s)`, `Carriers (del/dup)`, `beta (95% CI)`, P, `Exp.-wide`, Unrelated)]
disp[duplicated(Locus), Locus := ""]
cat("\n== Table 1 (", nrow(tab), "locus-trait associations from", nrow(x), "model-specific rows;", uniqueN(tab$Locus), "loci )\n"); print(disp)
fwrite(disp, file.path(od, "table1.csv"), bom = TRUE)
fwrite(tab[, .(Locus, phenotype = ph, Trait, model_shown, lead_probe, chr, lead_pos, cnvr_start = pos, beta, se, p, n_del, n_dup, n_del_at_lead, n_dup_at_lead, N,
               exp_wide = `Exp.-wide`, unrelated = Unrelated)], file.path(od, "table1_numeric.csv"))
cat("\nSaved:", file.path(od, "table1.csv"), "and table1_numeric.csv\n")
