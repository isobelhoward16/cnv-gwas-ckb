#!/usr/bin/env Rscript
# Stage 8, step 5: Table S4, carriers of 54 recurrent genomic-disorder CNVs in participants passing CNV QC,
# using the hg19 locus definitions of Crawford et al. (expected CNV type, >= 50% overlap between call and
# locus interval, and the locus-specific criteria). Also writes the CKB denominator used by step 6.
suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(yaml)
  library(ggplot2)
  library(scales)
})

# -----------------------------
# Config / paths
# -----------------------------
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)

out_root <- paths$out_root

# Annotated CNVs (from 04_exon_annotation.R)
annotated_dir <- paths$annotated_cnv_dir

# Paper1 denominator manifest (includes zero-CNV QC-pass samples)
manifest_path <- file.path(out_root, "data/derived/manifests/postQC_denominator_samples.tsv")

# Outputs (Paper1)
tables_derived <- file.path(out_root, "data/derived/tables")
tables_main    <- file.path(out_root, "results/tables/main")
tables_supp    <- file.path(out_root, "results/tables/supp")
logs_dir       <- file.path(out_root, "logs")

dir.create(tables_derived, recursive = TRUE, showWarnings = FALSE)
dir.create(tables_main,    recursive = TRUE, showWarnings = FALSE)
dir.create(tables_supp,    recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,       recursive = TRUE, showWarnings = FALSE)

run_ts <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
log_file <- file.path(logs_dir, paste0("03_pathogenic_identify_", run_ts, ".log"))
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)

cat("=== Paper1: UKB pathogenic CNVs (Crawford-style) in CKB ===\n")
cat("Started: ", format(Sys.time()), "\n")
cat("paths.yaml: ", paths_yaml, "\n")
cat("annotated_dir: ", annotated_dir, "\n")
cat("manifest: ", manifest_path, "\n")
cat("log: ", log_file, "\n\n")

# -----------------------------
# Denominators (Paper1)
# -----------------------------
denom_dt <- fread(manifest_path)
stopifnot("phase" %in% names(denom_dt))
# ---------------------------------------------------------
# Apply FINAL sample QC: exclude samples with >100 CNVs
# ---------------------------------------------------------

cat("--- Applying final NumCNV <= 100 QC ---\n")

stopifnot("NumCNV" %in% names(denom_dt))
stopifnot("ckid" %in% names(denom_dt))

cat("Manifest before NumCNV filter: ", nrow(denom_dt), "\n")

denom_dt <- denom_dt[
  !is.na(NumCNV) &
  NumCNV <= 100
]

cat("Manifest after NumCNV <= 100: ", nrow(denom_dt), "\n")

# CKB denominator, read by 06_figure_s6_ckb_vs_ukb.R
writeLines(as.character(nrow(denom_dt)), file.path(tables_derived, "ckb_denominator_full_qc.txt"))

final_qc_ckids <- unique(denom_dt$ckid)

cat("Final full-QC denominator confirmed: ",
    length(final_qc_ckids), "\n\n")

denom_overall <- nrow(denom_dt)
denom_by_phase <- denom_dt[, .N, by = phase]
setnames(denom_by_phase, "N", "denom")

cat("--- Denominators (Paper1 manifest) ---\n")
cat("Overall: ", denom_overall, "\n")
print(denom_by_phase)
cat("\n")

# -----------------------------
# Load annotated CNVs (EXPLICIT columns; no guessing)
# -----------------------------
files <- list.files(annotated_dir, pattern = "annotated.*\\.txt$", full.names = TRUE)
stopifnot(length(files) > 0)
cat("Annotated files found: ", length(files), "\n")

cnv_list <- lapply(files, function(f) {
  phase_tag <- gsub(".*annotated_merged_gene_mapped_(.*)\\.txt$", "\\1", basename(f))
  dt <- fread(f, showProgress = FALSE)

  # Explicit column usage (as confirmed):
  # chr,start,end,cnv_type,ckid,exonic_genes,genic_genes
  out <- dt[, .(
    chr         = as.character(chr),
    start       = as.integer(start),
    end         = as.integer(end),
    cnv_type    = as.character(cnv_type),
    ckid        = as.character(ckid),
    exonic_genes= as.character(exonic_genes),
    genic_genes = as.character(genic_genes),
    phase_tag   = phase_tag
  )]

  out[, phase := fifelse(grepl("Phase12", phase_tag), "Phase12",
                  fifelse(grepl("Phase3", phase_tag), "Phase3", NA_character_))]
  out
})

cnv <- rbindlist(cnv_list, use.names = TRUE, fill = TRUE)
cat("CNVs loaded (raw): ", nrow(cnv), "\n")

# Sanitize CNVs
cnv[, chr := sub("^chr", "", as.character(chr))]
cnv[, `:=`(start = as.integer(start), end = as.integer(end))]

# swap reversed intervals
swap_idx <- which(!is.na(cnv$start) & !is.na(cnv$end) & cnv$start > cnv$end)
if (length(swap_idx)) cnv[swap_idx, `:=`(start = end, end = start)]

# drop invalids
cnv <- cnv[!is.na(chr) & !is.na(start) & !is.na(end)]
cnv <- cnv[!is.na(phase)]

# normalize type
cnv[, cnv_type := tolower(cnv_type)]
cnv[cnv_type %chin% c("del","deletion","loss"), cnv_type := "deletion"]
cnv[cnv_type %chin% c("dup","duplication","gain"), cnv_type := "duplication"]
cnv <- cnv[cnv_type %chin% c("deletion","duplication")]

cat("CNVs after sanitize/type: ", nrow(cnv), "\n")
cat("Unique CKIDs in CNV table: ", uniqueN(cnv$ckid), "\n\n")

# Restrict CNV calls to the exact same final-QC population
cat("--- Restricting CNVs to final full-QC denominator ---\n")

cat("CNVs before full-QC restriction: ", nrow(cnv), "\n")
cat("Samples before full-QC restriction: ", uniqueN(cnv$ckid), "\n")

cnv <- cnv[ckid %in% final_qc_ckids]

cat("CNVs after full-QC restriction: ", nrow(cnv), "\n")
cat("Samples with >=1 CNV after restriction: ", uniqueN(cnv$ckid), "\n\n")

# -----------------------------
# Genomic-disorder CNV regions (Crawford et al., hg19)
# -----------------------------
cat("Setting up UKB CNV regions...\n")
ukb_regions <- data.table(
  cnv_locus = c(
    "TAR_del", "TAR_dup",
    "1q21.1del", "1q21.1dup",
    "NRXN1del",
    "2q11.2del", "2q11.2dup",
    "2q13del(NPHP1)", "2q13dup(NPHP1)",
    "2q13del", "2q13dup",
    "2q21.1del", "2q21.1dup",
    "3q29del", "3q29dup",
    "WBS_dup",
    "7q11.23dup_distal",
    "8p23.1dup",
    "10q11.21q11.23del", "10q11.21q11.23dup",
    "10q23dup",
    "13q12del(CRYL1)", "13q12dup(CRYL1)",
    "13q12.12del", "13q12.12dup",
    "15q11.2del", "15q11.2dup",
    "PWS_dup",
    "15q11q13del_BP3-BP4", "15q11q13dup_BP3-BP4",
    "15q11q13dup_BP3-BP5",
    "15q13.3del", "15q13.3dup",
    "15q13.3del(CHRNA7)", "15q13.3dup(CHRNA7)",
    "15q24dup",
    "16p13.11del", "16p13.11dup",
    "16p12.1del", "16p12.1dup",
    "16p11.2distal_del", "16p11.2distal_dup",
    "16p11.2del", "16p11.2dup",
    "17p12del(HNPP)", "17p12dup(CMT1A)",
    "Potocki-Lupski",
    "17q11.2del(NF1)",
    "17q12del", "17q12dup",
    "22q11.2del", "22q11.2dup",
    "22q11.2distal_del", "22q11.2distal_dup"
  ),
  chr = c(
    "1", "1",
    "1", "1",
    "2",
    "2", "2",
    "2", "2",
    "2", "2",
    "2", "2",
    "3", "3",
    "7",
    "7",
    "8",
    "10", "10",
    "10",
    "13", "13",
    "13", "13",
    "15", "15",
    "15",
    "15", "15",
    "15",
    "15", "15",
    "15", "15",
    "15",
    "16", "16",
    "16", "16",
    "16", "16",
    "16", "16",
    "17", "17",
    "17",
    "17",
    "17", "17",
    "22", "22",
    "22", "22"
  ),
  start = as.integer(c(
    145390000, 145390000,
    146530000, 146530000,
    50140000,
    96740000, 96740000,
    110860000, 110860000,
    111390000, 111390000,
    131480000, 131480000,
    195720000, 195720000,
    72740000,
    75140000,
    8100000,
    49390000, 49390000,
    82050000,
    20980000, 20980000,
    23560000, 23560000,
    22810000, 22810000,
    23680000,
    29160000, 29160000,
    29160000,
    31080000, 31080000,
    32020000, 32020000,
    72900000,
    15510000, 15510000,
    21950000, 21950000,
    28820000, 28820000,
    29650000, 29650000,
    14140000, 14140000,
    16810000,
    29120000,
    34810000, 34810000,
    19040000, 19040000,
    21920000, 21920000
  )),
  end = as.integer(c(
    145810000, 145810000,
    147390000, 147390000,
    51260000,
    97680000, 97680000,
    110980000, 110980000,
    112010000, 112010000,
    131930000, 131930000,
    197350000, 197350000,
    74140000,
    76060000,
    11870000,
    51060000, 51060000,
    88930000,
    21100000, 21100000,
    24880000, 24880000,
    23090000, 23090000,
    28390000,
    30380000, 30380000,
    32460000,
    32460000, 32460000,
    32450000, 32450000,
    78150000,
    16290000, 16290000,
    22430000, 22430000,
    29050000, 29050000,
    30200000, 30200000,
    15430000, 15430000,
    20210000,
    30270000,
    36220000, 36220000,
    21470000, 21470000,
    23650000, 23650000
  )),
  expected_type = c(
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "duplication",
    "duplication",
    "duplication",
    "deletion", "duplication",
    "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "duplication",
    "deletion", "duplication",
    "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication",
    "duplication",
    "deletion",
    "deletion", "duplication",
    "deletion", "duplication",
    "deletion", "duplication"
  )
)

cat("UKB regions defined: ", nrow(ukb_regions), " loci\n\n")

# -----------------------------
# Apply Overlap Filter
# -----------------------------
cat("Finding overlaps between CNVs and UKB regions...\n")
setkey(cnv, chr, start, end)
setkey(ukb_regions, chr, start, end)

hits <- foverlaps(cnv, ukb_regions, nomatch = 0L)
cat("Initial overlaps found: ", nrow(hits), "\n")

# Standardize column names for downstream processing
if ("i.start" %in% names(hits)) {
  setnames(hits, c("start", "end", "i.start", "i.end"),
           c("reg_start", "reg_end", "cnv_start", "cnv_end"))
} else {
  hits[, `:=`(cnv_start = start, cnv_end = end, reg_start = start, reg_end = end)]
}

# Filter by expected CNV type
cat("Applying CNV type filter...\n")
hits <- hits[cnv_type == expected_type]
cat("CNVs after type filtering: ", nrow(hits), "\n")

# -------- Apply 50% Overlap Rule --------
cat("Applying 50% overlap rule...\n")
hits[, `:=`(
  overlap_start = pmax(cnv_start, reg_start),
  overlap_end   = pmin(cnv_end, reg_end),
  region_length = reg_end - reg_start
)]
hits[, overlap_length := pmax(0, overlap_end - overlap_start)]
hits[, region_coverage := overlap_length / region_length]
hits_50pct <- hits[region_coverage >= 0.5]
cat("CNVs after 50% overlap filter: ", nrow(hits_50pct), "\n")

# -------- Apply the locus-specific criteria of Crawford et al. --------
cat("Applying UKB gene and size criteria...\n")

simple_overlap_cnvs <- c(
  "TAR_del", "TAR_dup", "1q21.1del", "1q21.1dup", "2q13del", "2q13dup",
  "2q21.1del", "2q21.1dup", "3q29del", "3q29dup", "WBS_dup", "7q11.23dup_distal",
  "10q11.21q11.23del", "10q11.21q11.23dup", "13q12.12del", "13q12.12dup",
  "15q11.2del", "15q11.2dup", "15q11q13del_BP3-BP4", "15q11q13dup_BP3-BP4",
  "15q11q13dup_BP3-BP5", "15q13.3del", "15q13.3dup", "16p13.11del", "16p13.11dup",
  "16p12.1del", "16p12.1dup", "16p11.2distal_del", "16p11.2distal_dup",
  "16p11.2del", "16p11.2dup", "Potocki-Lupski", "17q12del", "17q12dup",
  "22q11.2del", "22q11.2dup", "22q11.2distal_del", "22q11.2distal_dup", "15q24dup"
)

hits_50pct[, cnv_size := cnv_end - cnv_start]

hits_50pct[, `:=`(
  genic_list = lapply(genic_genes, function(x) {
    if (is.na(x) || x == "NOT_FOUND" || x == "") character(0)
    else trimws(strsplit(x, ",")[[1]])
  }),
  exonic_list = lapply(exonic_genes, function(x) {
    if (is.na(x) || x == "NOT_FOUND" || x == "") character(0)
    else trimws(strsplit(x, ",")[[1]])
  })
)]

hits_50pct[, passes_ukb := FALSE]

hits_50pct[cnv_locus %in% simple_overlap_cnvs, passes_ukb := TRUE]
hits_50pct[cnv_locus == "8p23.1dup" & cnv_size >= 1000000, passes_ukb := TRUE]

hits_50pct[cnv_locus == "PWS_dup", passes_ukb := {
  pws_region_length <- 28390000 - 23680000
  full_coverage <- (cnv_start <= 23680000 & cnv_end >= 28390000)
  covers_full_region <- overlap_length >= pws_region_length * 0.95
  full_coverage & covers_full_region
}]

hits_50pct[cnv_locus == "NRXN1del", passes_ukb := sapply(exonic_list, function(x) "NRXN1" %in% x)]
hits_50pct[cnv_locus == "13q12del(CRYL1)", passes_ukb := sapply(exonic_list, function(x) "CRYL1" %in% x)]

hits_50pct[cnv_locus == "13q12dup(CRYL1)", passes_ukb := sapply(genic_list, function(x) "CRYL1" %in% x)]
hits_50pct[cnv_locus == "2q11.2del", passes_ukb := sapply(genic_list, function(x) all(c("LMAN2L", "ARID5A") %in% x))]
hits_50pct[cnv_locus == "2q11.2dup", passes_ukb := sapply(genic_list, function(x) all(c("LMAN2L", "ARID5A") %in% x))]
hits_50pct[cnv_locus == "2q13del(NPHP1)", passes_ukb := sapply(genic_list, function(x) "NPHP1" %in% x)]
hits_50pct[cnv_locus == "2q13dup(NPHP1)", passes_ukb := sapply(genic_list, function(x) "NPHP1" %in% x)]
hits_50pct[cnv_locus == "15q13.3del(CHRNA7)", passes_ukb := sapply(genic_list, function(x) "CHRNA7" %in% x)]
hits_50pct[cnv_locus == "15q13.3dup(CHRNA7)", passes_ukb := sapply(genic_list, function(x) "CHRNA7" %in% x)]
hits_50pct[cnv_locus == "17p12del(HNPP)", passes_ukb := sapply(genic_list, function(x) "PMP22" %in% x)]
hits_50pct[cnv_locus == "17p12dup(CMT1A)", passes_ukb := sapply(genic_list, function(x) "PMP22" %in% x)]
hits_50pct[cnv_locus == "17q11.2del(NF1)", passes_ukb := sapply(genic_list, function(x) "NF1" %in% x)]

hits_50pct[cnv_locus == "10q23dup", passes_ukb := {
  size_ok <- cnv_size >= 1000000
  genes_ok <- sapply(genic_list, function(x) all(c("NRG3", "GRID1") %in% x))
  size_ok & genes_ok
}]

final_hits <- hits_50pct[passes_ukb == TRUE]
cat("CNVs after UKB criteria: ", nrow(final_hits), "\n\n")

# -----------------------------
# Paper1 summary tables with correct denominators
# -----------------------------
# Carrier counts (overall + by phase) using Paper1 denominators
ckb_counts_phase <- final_hits[, .(carriers = uniqueN(ckid)), by = .(phase, cnv_locus, chr, reg_start, reg_end, expected_type)]
ckb_counts_overall <- final_hits[, .(carriers = uniqueN(ckid)), by = .(cnv_locus, chr, reg_start, reg_end, expected_type)]
ckb_counts_overall[, phase := "Overall"]

ckb_counts <- rbindlist(list(
  ckb_counts_overall[, .(phase, cnv_locus, chr, reg_start, reg_end, expected_type, carriers)],
  ckb_counts_phase[,   .(phase, cnv_locus, chr, reg_start, reg_end, expected_type, carriers)]
), use.names = TRUE)

denom_map <- rbind(
  data.table(phase = "Overall", denom = denom_overall),
  denom_by_phase[, .(phase, denom)]
)

ckb_counts <- merge(ckb_counts, denom_map, by = "phase", all.x = TRUE)
ckb_counts[, freq := carriers / denom]
ckb_counts[, freq_percent := 100 * freq]

ckb_counts[, location_hg19 := sprintf("chr%s:%s–%s",
  chr,
  format(round(reg_start/1e6, 2), nsmall=2),
  format(round(reg_end/1e6, 2), nsmall=2)
)]

# ---- Ensure all 54 loci appear (fill zeros for loci with no passing CNVs) ----
all_loci <- unique(ukb_regions[, .(
  cnv_locus, expected_type, chr,
  reg_start = start,
  reg_end   = end
)])

all_loci[, location_hg19 := sprintf("chr%s:%s–%s",
                                    chr,
                                    format(round(reg_start/1e6, 2), nsmall=2),
                                    format(round(reg_end/1e6, 2), nsmall=2))]

# full grid of phases x loci
all_grid <- CJ(
  phase = c("Overall", "Phase12", "Phase3"),
  cnv_locus = all_loci$cnv_locus,
  unique = TRUE
)

all_grid <- merge(
  all_grid,
  all_loci[, .(cnv_locus, expected_type, chr, reg_start, reg_end, location_hg19)],
  by = "cnv_locus",
  all.x = TRUE
)

# join observed counts; fill missing carriers with 0
ckb_counts <- merge(
  all_grid,
  ckb_counts[, .(phase, cnv_locus, expected_type, chr, reg_start, reg_end, location_hg19, carriers)],
  by = c("phase","cnv_locus","expected_type","chr","reg_start","reg_end","location_hg19"),
  all.x = TRUE
)

ckb_counts[is.na(carriers), carriers := 0L]

# re-attach denominators + recompute freqs (safe)
ckb_counts <- merge(ckb_counts, denom_map, by = "phase", all.x = TRUE)
ckb_counts[, freq := carriers / denom]
ckb_counts[, freq_percent := 100 * freq]

# Wide table for paper (Overall + Phase12 + Phase3)
wide <- dcast(
  ckb_counts,
  cnv_locus + expected_type + location_hg19 ~ phase,
  value.var = c("carriers", "freq_percent"),
  fill = 0
)

# Order by genomic position (use Overall rows)
ord_ref <- unique(ckb_counts[phase == "Overall", .(cnv_locus, chr_num = as.integer(chr), reg_start)])
setorder(ord_ref, chr_num, reg_start)
wide <- merge(ord_ref[, .(cnv_locus)], wide, by = "cnv_locus", all.x = TRUE)

# -----------------------------
# Save outputs (Paper1)
# -----------------------------
out_hits        <- file.path(tables_derived, "ukb_pathogenic_cnvs_ckb.csv")
out_counts_long <- file.path(tables_derived, "ukb_pathogenic_counts_ckb_long.tsv")
out_table_main  <- file.path(tables_main,   "Table_pathogenic_loci_ckb.tsv")
out_table_supp  <- file.path(tables_supp,   "STable_pathogenic_loci_ckb_byphase.tsv")

fwrite(final_hits, out_hits)
fwrite(ckb_counts, out_counts_long, sep = "\t")
fwrite(wide, out_table_main, sep = "\t")
fwrite(ckb_counts, out_table_supp, sep = "\t")

cat("Saved:\n")
cat("  Raw passing CNVs: ", out_hits, "\n")
cat("  Long counts:      ", out_counts_long, "\n")
cat("  Main table:       ", out_table_main, "\n")
cat("  Supp table:       ", out_table_supp, "\n\n")

cat("Completed: ", format(Sys.time()), "\n")
sink()
close(zz)
