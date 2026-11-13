#!/usr/bin/env Rscript
# Stage 7, step 2: CKB and UK Biobank carrier frequencies at the reported CNVRs (Table S7).
# For each CNVR, takes the UK Biobank array probe positions inside it (Auwerx et al. 2022, SD1, hg19,
# N = 331,522 unrelated white British participants) and counts CKB participants whose deletion or duplication
# call covers each position; summarised by the maximum (and median) across positions.
# CKB denominator: unrelated participants in the analysis sample (passing CNV QC, unrelated == 1).
# Aggregate output only. Check small-count disclosure rules before sharing cells with fewer than 5 carriers.
# Output: results/07_external/<version>/ukb/
source("R/config.R"); suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
if (is.null(pa$ukb_auwerx_dir)) stop("Set ukb_auwerx_dir in config/paths.yaml (folder with the Auwerx SD1 and SD3-SD5 files)")
od <- file.path(pa$results_root, "07_external", ver, "ukb"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
flank_bp <- 0L

# ---- CNVRs from stage 5: per locus, the union of the CNVRs of all its signals ----
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x[, `:=`(s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end))]
x[, trait := sub("_rint$", "", phenotype)]
# locus labels as in Table 1 / Figure 3 (loci sharing a cytoband split by position; sitting height "1q21.2" is the 1q21.1 distal CNV)
x[, locus := locus_cytoband_final]
x[locus == "16p11.2", locus := fifelse(s < 29.5e6, "16p11.2 BP2-BP3", "16p11.2 BP4-BP5")]
x[locus %in% c("1q21.1", "1q21.2") & lead_pos > 145e6, locus := "1q21.1 distal"]
x[locus == "1q21.1" & lead_pos < 145e6, locus := "1q21.1 proximal"]
cnvr <- x[, .(traits = paste(sort(unique(trait)), collapse = ", "), models = paste(sort(unique(cnv_type)), collapse = "+"),
              start = min(s), end = max(e)), by = .(locus, chr)]
cnvr[, `:=`(qstart = as.integer(start - flank_bp), qend = as.integer(end + flank_bp))]
cat("CNVRs:", nrow(cnvr), "\n"); print(cnvr[, .(locus, traits, models, chr, start, end)])

# ---- CKB denominator: unrelated participants in the analysis sample ----
ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))
qc  <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
ph  <- readRDS(pa$phenotype_rds); setDT(ph)
k1 <- fread(pa$array1_key); k2 <- fread(pa$array2_key); setnames(k1, c("ccvid", "xIID_cnv")); setnames(k2, c("ccvid", "xIID_cnv"))
k1[, gwas_array_type := "Axiom_CKB1"]; k2[, gwas_array_type := "Axiom_CKB2"]
unrel <- merge(ph[unrelated == 1, .(ccvid, gwas_array_type)], rbind(k1, k2), by = c("ccvid", "gwas_array_type"))$xIID_cnv
denom <- Reduce(intersect, list(ids, qc, unrel)); N_ckb <- length(denom)
cat("\nCKB denominator (unrelated, passing CNV QC, in the", ver, "analysis sample):", N_ckb, "\n")

cs <- as.data.table(readRDS(file.path(pa$results_root, "01_callset", ver, "final", "callset.rds")))
cnv <- cs[CKID %in% denom, .(chr = suppressWarnings(as.integer(sub("^chr", "", Chromosome))), start = as.integer(StartPos), end = as.integer(EndPos),
                             type = fifelse(as.numeric(CN) < 2, "del", "dup"), ckid = CKID)][!is.na(chr) & chr %in% 1:22]
cat("CKB calls used:", nrow(cnv), "| carriers:", uniqueN(cnv$ckid), "\n")

# ---- UKB probes ----
ukb <- fread(file.path(pa$ukb_auwerx_dir, "SD1_Auwerx_2021_CNV_Frequency_UKBB.txt"), select = c("ID", "CHR", "POS", "NumCNV", "NumDup", "NumDel", "NumNeutral"))
ukb[, CHR := suppressWarnings(as.integer(sub("^chr", "", CHR)))]; ukb <- ukb[!is.na(CHR)]; ukb[, N := NumCNV + NumNeutral]
cat("UKB probes:", nrow(ukb), "| N per probe:", paste(range(ukb$N), collapse = "-"), "(expect 331522)\n")
chk <- ukb[CHR == 16 & POS >= 29650000 & POS <= 30200000]
cat(sprintf("Build check, UKB 16p11.2 BP4-BP5: %d probes, median dup %.4f%% (expect about 0.03%%)\n", nrow(chk), 100 * median(chk$NumDup / chk$N)))

# ---- positions inside each CNVR, CKB carriers covering each ----
pos <- ukb[, .(chr = CHR, p0 = POS, p1 = POS, ID, ukb_del = NumDel, ukb_dup = NumDup, ukb_N = N)]; setkey(pos, chr, p0, p1)
reg <- cnvr[, .(locus, chr, qstart, qend)]; setkey(reg, chr, qstart, qend)
pos_in <- foverlaps(pos, reg, by.x = c("chr", "p0", "p1"), type = "within", nomatch = 0L)
cnv_k <- cnv[, .(chr, start, end, type, ckid)]; setkey(cnv_k, chr, start, end)
hit <- foverlaps(pos_in[, .(locus, ID, chr, p0, p1)], cnv_k, by.x = c("chr", "p0", "p1"), type = "within", nomatch = 0L)
ckb_at <- dcast(hit[, .(n = uniqueN(ckid)), by = .(locus, ID, type)], locus + ID ~ type, value.var = "n", fill = 0L)
for (v in c("del", "dup")) if (!v %in% names(ckb_at)) ckb_at[, (v) := 0L]
setnames(ckb_at, c("del", "dup"), c("ckb_del", "ckb_dup"))
probe_tab <- merge(pos_in, ckb_at, by = c("locus", "ID"), all.x = TRUE)
probe_tab[is.na(ckb_del), ckb_del := 0L]; probe_tab[is.na(ckb_dup), ckb_dup := 0L]
cnvr_any <- foverlaps(cnv_k, reg, type = "any", nomatch = 0L)[, .(ckb_any_del = uniqueN(ckid[type == "del"]), ckb_any_dup = uniqueN(ckid[type == "dup"])), by = locus]

res <- probe_tab[, .(n_ukb_probes = .N,
  ckb_del_max_n = max(ckb_del), ukb_del_max_n = max(ukb_del), ckb_dup_max_n = max(ckb_dup), ukb_dup_max_n = max(ukb_dup),
  ckb_del_max_pct = 100 * max(ckb_del) / N_ckb, ukb_del_max_pct = 100 * max(ukb_del / ukb_N),
  ckb_dup_max_pct = 100 * max(ckb_dup) / N_ckb, ukb_dup_max_pct = 100 * max(ukb_dup / ukb_N),
  ckb_del_med_pct = 100 * median(ckb_del) / N_ckb, ukb_del_med_pct = 100 * median(ukb_del / ukb_N),
  ckb_dup_med_pct = 100 * median(ckb_dup) / N_ckb, ukb_dup_med_pct = 100 * median(ukb_dup / ukb_N)), by = locus]
res <- merge(cnvr[, .(locus, traits, models, chr, start, end)], res, by = "locus", all.x = TRUE)
res <- merge(res, cnvr_any, by = "locus", all.x = TRUE)
res[, `:=`(ratio_del_max = ckb_del_max_pct / ukb_del_max_pct, ratio_dup_max = ckb_dup_max_pct / ukb_dup_max_pct, N_ckb = N_ckb)]
num <- names(res)[sapply(res, is.double)]; res[, (num) := lapply(.SD, signif, 3), .SDcols = num]
setorder(res, chr, start)
cat("\n== Table S7 (max across shared UKB positions; NA = no UKB probe inside the CNVR)\n")
print(res[, .(locus, traits, models, n_ukb_probes, ckb_del_max_n, ukb_del_max_n, ckb_del_max_pct, ukb_del_max_pct, ratio_del_max,
              ckb_dup_max_n, ukb_dup_max_n, ckb_dup_max_pct, ukb_dup_max_pct, ratio_dup_max)])
fwrite(res, file.path(od, "table_s7_cnvr_freq_ckb_vs_ukb.csv")); fwrite(probe_tab, file.path(od, "table_s7_per_probe.csv"))
cat("\nSaved:", file.path(od, "table_s7_cnvr_freq_ckb_vs_ukb.csv"), "\nNB: check CKB small-count disclosure rules before sharing cells with < 5 carriers.\n")
