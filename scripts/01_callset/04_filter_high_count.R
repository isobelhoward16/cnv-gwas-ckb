#!/usr/bin/env Rscript
# Stage 1, step 4: remove samples with more than max_cnv calls after call filtering
# (calls counted per sample across all files; exact duplicate calls counted once).
#
# Input:  03_filtered/calls_filtered.rds
# Output: 04_sample_filtered/calls_sample_filtered.rds
#         04_sample_filtered/high_cnv_excluded_ckids.txt
source("R/config.R")
cfg <- load_config(); max_cnv <- cfg$pipeline$sample_qc$max_cnv
STEP <- "04_filter_high_count"
start_step(cfg, STEP)
out_dir <- callset_dir(cfg, "04_sample_filtered")

x <- readRDS(callset_dir(cfg, "03_filtered", "calls_filtered.rds"))
record(cfg, STEP, "calls in", nrow(x))
record(cfg, STEP, "samples with calls", uniqueN(x$CKID))
record(cfg, STEP, "samples with calls in more than one batch file", x[, uniqueN(batch), by = CKID][V1 > 1L, .N])

n_per <- unique(x[, .(CKID, Chromosome, StartPos, EndPos, CN)])[, .N, by = CKID]
high <- sort(n_per[N > max_cnv]$CKID)
record(cfg, STEP, sprintf("samples with > %s calls (removed)", max_cnv), length(high))
record(cfg, STEP, "calls removed with them", sum(x$CKID %in% high))
x <- x[!CKID %in% high]
record(cfg, STEP, "calls kept", nrow(x))
record(cfg, STEP, "samples with calls kept", uniqueN(x$CKID))

saveRDS(x, file.path(out_dir, "calls_sample_filtered.rds"))
writeLines(high, file.path(out_dir, "high_cnv_excluded_ckids.txt"))
