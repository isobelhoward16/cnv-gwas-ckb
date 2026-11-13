#!/usr/bin/env Rscript
# Stage 1, step 6: freeze the call set.
#   - removes exact duplicate calls (same sample, chromosome, start, end and copy number);
#   - writes one call file per region_phase stratum, with a header row;
#   - writes the sample lists, including QC-passing samples with no CNV;
#   - writes a verification report and MD5 checksums.
#
# Input:  05_annotated/calls_annotated.rds; 02_sample_qc/qc_passing_ckids.txt;
#         04_sample_filtered/high_cnv_excluded_ckids.txt; logs/step_counts.tsv
# Output: final/calls/merged_gene_mapped_<region_phase>.txt
#         final/callset.rds
#         final/samples/{qc_passing,high_cnv_excluded,carriers,zero_cnv}_ckids.txt
#         final/verification_report.txt, final/MD5SUMS.txt
source("R/config.R")
cfg <- load_config()
STEP <- "06_freeze"
fin <- callset_dir(cfg, "final")
if (length(list.files(fin, recursive = TRUE)) && Sys.getenv("FORCE") != "1")
  stop("final/ already contains files. Callset ", cfg$pipeline$callset_version, " is frozen; set FORCE=1 only if you mean to rebuild it.")
start_step(cfg, STEP)
dir.create(file.path(fin, "calls"), showWarnings = FALSE); dir.create(file.path(fin, "samples"), showWarnings = FALSE)

x <- readRDS(callset_dir(cfg, "05_annotated", "calls_annotated.rds"))
record(cfg, STEP, "calls in", nrow(x))
key <- c("CKID", "Chromosome", "StartPos", "EndPos", "CN")
dup <- duplicated(x, by = key)
dups <- x[x[dup, ..key], on = key]
record(cfg, STEP, "exact duplicate calls removed", sum(dup))
record(cfg, STEP, "of which the copies are in different batch files", dups[, uniqueN(batch) > 1L, by = key][V1 == TRUE, .N])
x <- x[!dup]

out_cols <- c("Region", "NumSNP", "Length", "File", "StartSNP", "EndSNP", "Conf", "State", "CN",
             "Chromosome", "StartPos", "EndPos", "CN_type", "n_merged", "CKID", "Genes", "Gene_distance")
setorder(x, region_phase, CKID, Chromosome, StartPos)
for (rp in sort(unique(x$region_phase)))
  fwrite(x[region_phase == rp, ..out_cols], file.path(fin, "calls", paste0("merged_gene_mapped_", rp, ".txt")), sep = "\t")
saveRDS(x[, c(out_cols, "region_phase", "batch", "mixed_cn"), with = FALSE], file.path(fin, "callset.rds"))

passing <- readLines(callset_dir(cfg, "02_sample_qc", "qc_passing_ckids.txt"))
high    <- readLines(callset_dir(cfg, "04_sample_filtered", "high_cnv_excluded_ckids.txt"))
carriers <- sort(unique(x$CKID))
zero <- sort(setdiff(passing, c(carriers, high)))
if (length(setdiff(carriers, passing))) stop("Some carriers are not QC-passing samples")
if (length(intersect(carriers, high))) stop("Some high-count samples still have calls")
if (length(passing) != length(carriers) + length(zero) + length(intersect(high, passing))) stop("Sample counts do not add up")
writeLines(passing, file.path(fin, "samples", "qc_passing_ckids.txt"))
writeLines(high, file.path(fin, "samples", "high_cnv_excluded_ckids.txt"))
writeLines(carriers, file.path(fin, "samples", "carriers_ckids.txt"))
writeLines(zero, file.path(fin, "samples", "zero_cnv_ckids.txt"))

autos <- sub("^chr", "", x$Chromosome) %in% as.character(1:22)
record(cfg, STEP, "FINAL calls", nrow(x))
record(cfg, STEP, "FINAL calls on autosomes", sum(autos))
record(cfg, STEP, "FINAL deletions", sum(x$CN_type == "del"))
record(cfg, STEP, "FINAL duplications", sum(x$CN_type == "dup"))
for (cn in sort(unique(x$CN))) record(cfg, STEP, sprintf("FINAL calls with CN = %s", cn), sum(x$CN == cn))
record(cfg, STEP, "FINAL merged calls (n_merged > 1)", sum(x$n_merged > 1L))
record(cfg, STEP, "FINAL median length (bp)", median(x$Length))
record(cfg, STEP, "FINAL region_phase files", uniqueN(x$region_phase))
record(cfg, STEP, "SAMPLES QC-passing", length(passing))
record(cfg, STEP, "SAMPLES removed for high CNV count", length(high))
record(cfg, STEP, "SAMPLES carriers (>= 1 call)", length(carriers))
record(cfg, STEP, "SAMPLES carriers of an autosomal call", uniqueN(x$CKID[autos]))
record(cfg, STEP, "SAMPLES zero-CNV (QC-passing, no call)", length(zero))

# ---- verification report and checksums ----
cnt <- fread(counts_file(cfg), sep = "\t", colClasses = "character")
rep_file <- file.path(fin, "verification_report.txt")
con <- file(rep_file, "w")
writeLines(c(sprintf("Callset %s verification report, %s", cfg$pipeline$callset_version, format(Sys.time(), "%Y-%m-%d %H:%M")),
             "", "Settings (config/pipeline.yaml):", paste0("  ", strsplit(as.yaml(cfg$pipeline), "\n")[[1]]), ""), con)
for (s in unique(cnt$step_name)) {
  writeLines(paste0("[", s, "]"), con)
  y <- cnt[step_name == s]
  writeLines(sprintf("  %-62s %s", y$metric, y$value), con)
  writeLines("", con)
}
close(con)
fl <- list.files(fin, recursive = TRUE, full.names = TRUE); fl <- fl[!basename(fl) %in% "MD5SUMS.txt"]
md5 <- tools::md5sum(fl)
writeLines(paste(md5, sub(paste0("^", fin, "/"), "", names(md5))), file.path(fin, "MD5SUMS.txt"))
cat("\n   Report:", rep_file, "\n   Checksums:", file.path(fin, "MD5SUMS.txt"), "\n")
