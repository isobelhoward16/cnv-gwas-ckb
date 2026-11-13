#!/usr/bin/env Rscript
# Stage 1, step 3: call-level filters.
#   (a) Conf >= min_conf, NumSNP >= min_snps and min_length <= Length <= max_length;
#   (b) exclusion regions: a call is removed when >= exclusion_min_fraction of its length lies in
#       one exclusion region, after overlapping regions are merged (PennCNV scan_region.pl rule).
#
# Input:  02_sample_qc/calls_qc.rds; config/exclusion_regions_hg19.txt
# Output: 03_filtered/calls_filtered.rds
#         03_filtered/excluded_calls.tsv           calls removed by the exclusion regions
#         03_filtered/exclusion_region_counts.tsv  calls removed per exclusion region
source("R/config.R"); source("R/callset_utils.R")
cfg <- load_config(); fp <- cfg$pipeline$cnv_filter
STEP <- "03_filter_calls"
start_step(cfg, STEP)
out_dir <- callset_dir(cfg, "03_filtered")

x <- readRDS(callset_dir(cfg, "02_sample_qc", "calls_qc.rds"))
record(cfg, STEP, "calls in", nrow(x))

# ---- (a) metric filters ----
# a missing value counts as a fail
x[, `:=`(fail_conf   = !((Conf   >= fp$min_conf)   %in% TRUE), fail_snps   = !((NumSNP >= fp$min_snps)   %in% TRUE),
         fail_minlen = !((Length >= fp$min_length) %in% TRUE), fail_maxlen = !((Length <= fp$max_length) %in% TRUE))]
x[, pass_metrics := !(fail_conf | fail_snps | fail_minlen | fail_maxlen)]
record(cfg, STEP, sprintf("fail Conf < %s", fp$min_conf), sum(x$fail_conf))
record(cfg, STEP, sprintf("fail NumSNP < %s", fp$min_snps), sum(x$fail_snps))
record(cfg, STEP, sprintf("fail Length < %s", fp$min_length), sum(x$fail_minlen))
record(cfg, STEP, sprintf("fail Length > %s", fp$max_length), sum(x$fail_maxlen))
big <- x[fail_maxlen & !fail_conf & !fail_snps]
record(cfg, STEP, "calls > max_length that pass Conf and NumSNP (removed)", nrow(big))
record(cfg, STEP, "samples with such a call", uniqueN(big$CKID))
record(cfg, STEP, "calls passing metric filters", sum(x$pass_metrics))

# ---- (b) exclusion regions ----
reg_file <- file.path(cfg$root, cfg$pipeline$exclusion_regions_file)
reg <- read_regions(reg_file)
regc <- condense_regions(reg)
record(cfg, STEP, "exclusion file", cfg$pipeline$exclusion_regions_file)
record(cfg, STEP, "exclusion regions (lines)", nrow(reg))
record(cfg, STEP, "exclusion regions after merging overlaps", nrow(regc))

ov <- max_region_fraction(x$Chromosome, x$StartPos, x$EndPos, regc)
x[, `:=`(excl_frac = ov$frac, excl_hit = ov$hit)]
x[, in_excluded_region := excl_frac >= fp$exclusion_min_fraction]
x[, keep := pass_metrics & !in_excluded_region]
record(cfg, STEP, "calls passing metrics but removed by exclusion regions", sum(x$pass_metrics & x$in_excluded_region))
record(cfg, STEP, "calls kept", sum(x$keep))

excl <- x[pass_metrics & in_excluded_region]
fwrite(excl[, .(CKID, Region, CN_type, NumSNP, Conf, n_merged, excl_frac = round(excl_frac, 3), excl_hit)],
       file.path(out_dir, "excluded_calls.tsv"), sep = "\t")
reg_counts <- excl[, .(calls_removed = .N, del = sum(CN_type == "del"), dup = sum(CN_type == "dup")), by = .(region = excl_hit)][order(-calls_removed)]
fwrite(reg_counts, file.path(out_dir, "exclusion_region_counts.tsv"), sep = "\t")
cat("\n   Calls removed per exclusion region (top 15):\n"); print(head(reg_counts, 15))

drop_cols <- c("fail_conf", "fail_snps", "fail_minlen", "fail_maxlen", "pass_metrics", "in_excluded_region", "keep", "excl_hit")
saveRDS(x[keep == TRUE][, (drop_cols) := NULL][], file.path(out_dir, "calls_filtered.rds"))
