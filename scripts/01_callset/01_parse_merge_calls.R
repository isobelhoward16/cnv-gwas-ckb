#!/usr/bin/env Rscript
# Stage 1, step 1: read the raw PennCNV calls and merge fragmented calls.
#
# Within one sample, chromosome and CNV type, calls are taken in order of start position and the
# next call is merged into the current one when gap / span <= gap_fraction (config/pipeline.yaml).
# A merged call takes its coordinates from the outermost fragments; NumSNP and Conf are combined
# across fragments and State/CN are taken from the fragment with most probes.
#
# Input:  PennCNV call files (*_cnv_calls_batch_*.txt) under paths$raw_calls_dir
# Output: 01_merged/calls_merged.rds    one row per call after merging
#         01_merged/merge_details.tsv   fragments of each merged call
#         01_merged/batch_summary.tsv   calls per file before and after merging
source("R/config.R"); source("R/callset_utils.R")
cfg <- load_config(); mp <- cfg$pipeline$merge
STEP <- "01_parse_merge"
start_step(cfg, STEP)
out_dir <- callset_dir(cfg, "01_merged")
raw_dir <- normalizePath(cfg$paths$raw_calls_dir)

files <- list.files(raw_dir, pattern = "_cnv_calls_batch_.*\\.txt$", recursive = TRUE, full.names = TRUE)
rel <- substring(files, nchar(raw_dir) + 2)                     # skip non-call subfolders
files <- files[!grepl("logs|scripts|cnv_plots_all_regions|identical_cnv_lines", rel)]
dup <- basename(files)[duplicated(basename(files))]
if (length(dup)) stop("The same batch file name occurs in more than one folder: ", paste(head(dup, 10), collapse = ", "))
record(cfg, STEP, "batch files found", length(files))
record(cfg, STEP, "merge rule: gap_fraction", mp$gap_fraction)
record(cfg, STEP, "merge rule: NumSNP of merged call", mp$numsnp_rule)
record(cfg, STEP, "merge rule: Conf of merged call", mp$conf_rule)

res <- vector("list", length(files)); det <- vector("list", length(files)); bs <- vector("list", length(files))
n_empty <- 0L; n_bad_coord <- 0L; n_no_type <- 0L
for (k in seq_along(files)) {
  f <- files[k]; b <- basename(f)
  d <- parse_penncnv(f)
  if (is.null(d)) { n_empty <- n_empty + 1L; next }
  n_raw <- nrow(d)
  bad <- is.na(d$Chromosome) | is.na(d$StartPos) | is.na(d$EndPos)
  n_bad_coord <- n_bad_coord + sum(bad); d <- d[!bad]
  d[, CN_type := fifelse(CN %in% c(0L, 1L), "del", fifelse(!is.na(CN) & CN >= 3L, "dup", NA_character_))]
  n_no_type <- n_no_type + sum(is.na(d$CN_type)); d <- d[!is.na(CN_type)]
  if (isTRUE(mp$enabled)) {
    m <- merge_fragments(d, mp$gap_fraction, mp$numsnp_rule, mp$conf_rule)
  } else {
    d[, `:=`(n_merged = 1L, mixed_cn = FALSE)]
    m <- list(calls = d, details = NULL)
  }
  rp <- sub("_cnv_calls_batch.*$", "", b)
  m$calls[, `:=`(CKID = extract_ckid(File), region_phase = rp, batch = b)]
  res[[k]] <- m$calls
  if (!is.null(m$details) && nrow(m$details)) det[[k]] <- m$details[, batch := b]
  bs[[k]] <- data.table(batch = b, region_phase = rp, raw_calls = n_raw, calls_after_merge = nrow(m$calls),
                        merged_calls = sum(m$calls$n_merged > 1L))
  if (k %% 100 == 0) cat("   processed", k, "of", length(files), "files\n")
}
calls <- rbindlist(res); details <- rbindlist(det); bsum <- rbindlist(bs)

record(cfg, STEP, "empty batch files skipped", n_empty)
record(cfg, STEP, "raw calls read", sum(bsum$raw_calls))
record(cfg, STEP, "raw calls dropped: unreadable coordinates", n_bad_coord)
record(cfg, STEP, "raw calls dropped: CN not 0, 1 or >= 3 (no del/dup type)", n_no_type)
record(cfg, STEP, "calls after merging", nrow(calls))
record(cfg, STEP, "merged calls (n_merged > 1)", sum(calls$n_merged > 1L))
record(cfg, STEP, "fragments absorbed into merged calls", sum(calls$n_merged) - nrow(calls))
record(cfg, STEP, "merged calls whose fragments differ in CN", sum(calls$mixed_cn))
record(cfg, STEP, "calls with no sample ID extracted from File", sum(is.na(calls$CKID)))
record(cfg, STEP, "region_phase groups", uniqueN(calls$region_phase))
multi_file <- calls[, .(n = uniqueN(File)), by = CKID][n > 1L]
record(cfg, STEP, "sample IDs with calls under more than one File path", nrow(multi_file))

saveRDS(calls, file.path(out_dir, "calls_merged.rds"))
fwrite(details, file.path(out_dir, "merge_details.tsv"), sep = "\t")
fwrite(bsum, file.path(out_dir, "batch_summary.tsv"), sep = "\t")
cat("   saved:", file.path(out_dir, "calls_merged.rds"), "\n")
