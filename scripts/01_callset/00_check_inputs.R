#!/usr/bin/env Rscript
# Checks that every input needed by stage 1 exists and has the expected format.
# Quick; safe to run on a login node. Prints no sample IDs or call records.
source("R/config.R"); source("R/callset_utils.R")
cfg <- load_config(); pa <- cfg$paths; ok <- TRUE
say <- function(pass, what, detail = "") { cat(sprintf("  [%s] %s %s\n", if (pass) " OK " else "FAIL", what, detail)); if (!pass) ok <<- FALSE }

cat("R", R.version$major, ".", R.version$minor, " | data.table ", as.character(packageVersion("data.table")),
    " | yaml ", as.character(packageVersion("yaml")), "\n\nFiles and folders:\n", sep = "")
for (p in c("raw_calls_dir", "master_qc_summary", "refgene", "reflink", "penncnv_dir"))
  say(!is.null(pa[[p]]) && file.exists(pa[[p]]), p, pa[[p]])
say(file.exists(file.path(pa$penncnv_dir, "scan_region.pl")), "scan_region.pl")
say(!is.null(pa$results_root), "results_root", pa$results_root)

cat("\nRaw calls:\n")
raw_dir <- normalizePath(pa$raw_calls_dir)
files <- list.files(raw_dir, pattern = "_cnv_calls_batch_.*\\.txt$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("logs|scripts|cnv_plots_all_regions|identical_cnv_lines", substring(files, nchar(raw_dir) + 2))]
say(length(files) > 0, "batch files", length(files))
say(!anyDuplicated(basename(files)), "batch file names unique")
rp <- sub("_cnv_calls_batch.*$", "", basename(files))
say(uniqueN(rp) == 20, "region_phase groups", uniqueN(rp))
say(length(strsplit(readLines(files[1], n = 1), "[[:space:]]+")[[1]]) >= 8, "first batch file has the PennCNV call format")

cat("\nSample QC summary:\n")
qc <- fread(pa$master_qc_summary, sep = "\t", nrows = 5)
miss <- setdiff(c("CKID", "LRR_SD", "BAF_drift", "WF", "region_phase"), names(qc))
say(!length(miss), "QC columns present", if (length(miss)) paste("missing:", paste(miss, collapse = ", ")) else "")

cat("\nExclusion regions:\n")
reg <- tryCatch(read_regions(file.path(cfg$root, cfg$pipeline$exclusion_regions_file)), error = function(e) { say(FALSE, conditionMessage(e)); NULL })
if (!is.null(reg)) {
  say(nrow(reg) == 70, "regions parsed", nrow(reg))
  imm <- c("chr14:22090057-23021075", "chr14:106052774-107288051", "chr2:89156874-90274235", "chr22:22380474-23265085")
  say(all(imm %in% reg$region), "hg19 TRA, IGH, IGK, IGL rows present")
  say(nrow(condense_regions(reg)) <= nrow(reg), "regions after merging overlaps", nrow(condense_regions(reg)))
}

cat("\nGene annotation files (first 2 lines):\n")
cat("  refgene: ", paste(substr(readLines(pa$refgene, n = 2), 1, 110), collapse = "\n           "), "\n")
cat("  reflink: ", paste(substr(readLines(pa$reflink, n = 2), 1, 110), collapse = "\n           "), "\n")
cat("  reflink columns:", length(strsplit(readLines(pa$reflink, n = 1), "\t")[[1]]), "\n")

out <- file.path(pa$results_root, "01_callset", cfg$pipeline$callset_version, "final")
say(!(dir.exists(out) && length(list.files(out, recursive = TRUE))), "callset not yet frozen", out)
cat(if (ok) "\nAll checks passed: ready to submit.\n" else "\nSome checks FAILED: fix these before submitting.\n")
