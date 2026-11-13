#!/usr/bin/env Rscript
# Stage 1, step 2: sample QC. Samples pass if LRR_SD <= max_lrr_sd, |BAF_drift| <= max_baf_drift
# and |WF| <= max_wf in the PennCNV QC summary; calls of other samples are dropped.
#
# Input:  01_merged/calls_merged.rds; paths$master_qc_summary
# Output: 02_sample_qc/calls_qc.rds
#         02_sample_qc/qc_passing_ckids.txt
#         02_sample_qc/qc_summary_by_region_phase.tsv
source("R/config.R")
cfg <- load_config(); q <- cfg$pipeline$sample_qc
STEP <- "02_sample_qc"
start_step(cfg, STEP)
out_dir <- callset_dir(cfg, "02_sample_qc")

qc <- fread(cfg$paths$master_qc_summary, sep = "\t")
check_cols(qc, c("CKID", "LRR_SD", "BAF_drift", "WF", "region_phase"), "master QC summary")
record(cfg, STEP, "samples in master QC summary", nrow(qc))
record(cfg, STEP, "duplicated sample IDs in master QC summary", sum(duplicated(qc$CKID)))

f_lrr <- qc$LRR_SD > q$max_lrr_sd
f_baf <- !f_lrr & abs(qc$BAF_drift) > q$max_baf_drift
f_wf  <- !f_lrr & !f_baf & abs(qc$WF) > q$max_wf
na_qc <- is.na(qc$LRR_SD) | is.na(qc$BAF_drift) | is.na(qc$WF)
record(cfg, STEP, sprintf("failed LRR_SD > %s", q$max_lrr_sd), sum(f_lrr, na.rm = TRUE))
record(cfg, STEP, sprintf("then failed |BAF_drift| > %s", q$max_baf_drift), sum(f_baf, na.rm = TRUE))
record(cfg, STEP, sprintf("then failed |WF| > %s", q$max_wf), sum(f_wf, na.rm = TRUE))
record(cfg, STEP, "missing a QC metric (fail)", sum(na_qc))
signal_pass <- qc[!na_qc & LRR_SD <= q$max_lrr_sd & abs(BAF_drift) <= q$max_baf_drift & abs(WF) <= q$max_wf]
record(cfg, STEP, "passed signal QC", nrow(signal_pass))

passing <- sort(unique(signal_pass$CKID))
record(cfg, STEP, "QC-passing samples", length(passing))

summ <- qc[, .(initial = .N), by = region_phase][
  signal_pass[, .(passing = .N), by = region_phase], on = "region_phase", passing := i.passing][
  is.na(passing), passing := 0L][, pct_retained := round(100 * passing / initial, 1)][order(region_phase)]
print(summ)

calls <- readRDS(callset_dir(cfg, "01_merged", "calls_merged.rds"))
record(cfg, STEP, "calls in", nrow(calls))
record(cfg, STEP, "samples with calls", uniqueN(calls$CKID))
record(cfg, STEP, "samples with calls but not in the QC summary", length(setdiff(unique(calls$CKID), qc$CKID)))
calls <- calls[CKID %in% passing]
record(cfg, STEP, "calls kept (QC-passing samples)", nrow(calls))
record(cfg, STEP, "QC-passing samples with at least one call", uniqueN(calls$CKID))
record(cfg, STEP, "QC-passing samples with no raw calls at all", length(setdiff(passing, calls$CKID)))

saveRDS(calls, file.path(out_dir, "calls_qc.rds"))
writeLines(passing, file.path(out_dir, "qc_passing_ckids.txt"))
fwrite(summ, file.path(out_dir, "qc_summary_by_region_phase.tsv"), sep = "\t")
