#!/usr/bin/env Rscript
# Stage 8, step 7: Table S2, baseline characteristics of the analysis population and N per trait.
# Analysis population: samples in the stage 2 matrices that pass CNV QC and have phenotype data.
# Checks that each trait N equals the N used in the association analysis (all participants) and prints PASS/FAIL.
# Aggregate output only. Output: results/08_descriptive/<version>/tables/table_s2_population.csv
source("R/config.R")
suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
out_dir <- file.path(pa$results_root, "08_descriptive", ver, "tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Load data ----
pheno <- readRDS(pa$phenotype_rds); setDT(pheno)
qc_master <- fread(pa$qc_gwas_file)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
key1 <- fread(pa$array1_key); key2 <- fread(pa$array2_key)
setnames(key1, c("ccvid", "xIID_cnv")); setnames(key2, c("ccvid", "xIID_cnv"))
key1[, gwas_array_type := "Axiom_CKB1"]; key2[, gwas_array_type := "Axiom_CKB2"]
key_all <- rbind(key1, key2)
mat_ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))

# ---- 2. Sample attrition ----
cat("=======================================================\n SAMPLE ATTRITION (", ver, ")\n=======================================================\n")
cat(sprintf("Samples in the %s matrices (callset after sample QC and <= 100 calls): %d\n", ver, length(mat_ids)))
qc_m <- qc_master[GID %in% mat_ids]
cat(sprintf("  Not in master QC file:         %d\n", sum(!mat_ids %in% qc_master$GID)))
cat(sprintf("  Fail LRR_SD > 0.35:            %d\n", sum(qc_m$LRR_SD > 0.35, na.rm = TRUE)))
cat(sprintf("  Fail |WF| > 0.05:              %d\n", sum(abs(qc_m$WF) > 0.05, na.rm = TRUE)))
cat(sprintf("  Fail BAF_drift > 0.01:         %d\n", sum(qc_m$BAF_drift > 0.01, na.rm = TRUE)))
cat(sprintf("  Fail raw NumCNV > 100:         %d\n", sum(qc_m$NumCNV > 100, na.rm = TRUE)))
pass_ids <- intersect(mat_ids, qc_master[qc_pass == TRUE]$GID)
cat(sprintf("Pass all CNV QC filters:         %d\n", length(pass_ids)))
# match each CNV ID to its participant through the key for its own array (as in the GWAS)
k <- key_all[xIID_cnv %in% pass_ids]
pheno_qc <- merge(pheno[gwas_array_type %in% c("Axiom_CKB1", "Axiom_CKB2")], k, by = c("ccvid", "gwas_array_type"))
cat(sprintf("  No phenotype data:             %d\n", length(pass_ids) - uniqueN(pheno_qc$xIID_cnv)))
if (anyDuplicated(pheno_qc$ccvid)) stop("Duplicate participants in the analysis population")
n_final <- nrow(pheno_qc)
cat(sprintf("Final analysis population:       %d\n", n_final))
cat("=======================================================\n\n")

# ---- 3. Reference objects ----
region_labels <- c("12" = "Qingdao (Urban)", "16" = "Harbin (Urban)", "26" = "Haikou (Urban)", "36" = "Suzhou (Urban)",
                   "46" = "Liuzhou (Urban)", "52" = "Sichuan (Rural)", "58" = "Gansu (Rural)", "68" = "Henan (Rural)",
                   "78" = "Zhejiang (Rural)", "88" = "Hunan (Rural)")
traits <- list("BMI (kg/m2)" = "bmi_calc", "BMI at age 25 (kg/m2)" = "bmi_at_25", "Standing height (mm)" = "standing_height_mm",
               "Sitting height (mm)" = "sitting_height_mm", "Weight (kg)" = "weight_kg", "Waist circumference (mm)" = "waist_mm",
               "Hip circumference (mm)" = "hip_mm", "Waist-hip ratio" = "waist_hip_ratio", "Body fat percentage (%)" = "fat_percent",
               "Systolic blood pressure (mmHg)" = "sbp_mean", "Diastolic blood pressure (mmHg)" = "dbp_mean",
               "Resting heart rate (bpm)" = "heart_rate_mean", "Random plasma glucose (mmol/L)" = "random_glucose")
fmt_mean_sd <- function(x) { x <- x[!is.na(x)]; sprintf("%.2f (%.2f)", mean(x), sd(x)) }

# ---- 4. Build table rows ----
rows <- list()
rows[["age"]] <- data.table(Characteristic = "Age (years)", Statistic = fmt_mean_sd(pheno_qc$age), N = sum(!is.na(pheno_qc$age)))
n_female <- sum(pheno_qc$sex == 1, na.rm = TRUE)
rows[["sex"]] <- data.table(Characteristic = "Sex (% female)", Statistic = sprintf("%.1f%%", 100 * n_female / n_final), N = n_female)
rows[["region_header"]] <- data.table(Characteristic = "Recruitment region, n (%)", Statistic = "", N = NA_integer_)
region_tab <- pheno_qc[, .N, by = region_code][order(region_code)]
region_tab[, label := region_labels[as.character(region_code)]]
for (i in seq_len(nrow(region_tab)))
  rows[[paste0("region_", i)]] <- data.table(Characteristic = paste0("  ", region_tab$label[i]),
                                             Statistic = sprintf("%.1f%%", 100 * region_tab$N[i] / n_final), N = region_tab$N[i])
rows[["array_header"]] <- data.table(Characteristic = "Genotyping array, n (%)", Statistic = "", N = NA_integer_)
array_tab <- pheno_qc[, .N, by = gwas_array_type][order(gwas_array_type)]
for (i in seq_len(nrow(array_tab)))
  rows[[paste0("array_", i)]] <- data.table(Characteristic = paste0("  ", array_tab$gwas_array_type[i]),
                                            Statistic = sprintf("%.1f%%", 100 * array_tab$N[i] / n_final), N = array_tab$N[i])
rows[["traits_header"]] <- data.table(Characteristic = "Quantitative traits, mean (SD)", Statistic = "", N = NA_integer_)
for (label in names(traits)) {
  col <- traits[[label]]; x <- pheno_qc[[col]]
  rows[[paste0("trait_", col)]] <- data.table(Characteristic = paste0("  ", label), Statistic = fmt_mean_sd(x), N = sum(!is.na(x)))
}
supp_table <- rbindlist(rows, fill = TRUE)

# ---- 5. Check trait N against the GWAS N (all participants; RINT column non-missing, as in the GWAS) ----
meff <- rbindlist(lapply(list.files(file.path(pa$results_root, "03_meff", ver, "tasks"), "^meff_summary_.*\\.csv$", full.names = TRUE), fread))
gN <- meff[Relatedness == "all" & CNV_Type == "mirror", .(col = sub("_rint$", "", Phenotype), N_gwas = N_Samples)]
chk <- data.table(col = unlist(traits))[, N_table := sapply(col, function(cc) sum(!is.na(pheno_qc[[paste0(cc, "_rint")]])))][gN, on = "col"]
chk[, N_raw := sapply(col, function(cc) sum(!is.na(pheno_qc[[cc]])))]
ok <- all(chk$N_table == chk$N_gwas) && all(chk$N_raw == chk$N_gwas)
cat(sprintf("[%s] trait N in Table S2 = GWAS N for all %d traits\n", if (ok) "PASS" else "FAIL", nrow(chk)))
if (!ok) print(chk)

# ---- 6. Print and write ----
cat("=======================================================\n TABLE S2: BASELINE CHARACTERISTICS\n")
cat(sprintf(" Analysis population: n = %d\n", n_final))
cat("=======================================================\n")
print(supp_table, nrows = Inf)
f <- file.path(out_dir, "table_s2_population.csv")
fwrite(supp_table, f, na = "")
cat("\nSaved:", f, "\n")
