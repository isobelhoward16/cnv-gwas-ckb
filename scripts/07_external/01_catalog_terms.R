#!/usr/bin/env Rscript
# Stage 7, step 1: GWAS Catalog trait names used to annotate the CNVRs (Table S3). Checks the exact trait
# names in config/catalog_terms.tsv against the Catalog file and counts their associations (all, and autosomal
# with P < 5e-8, the filter used for the CNVR look-ups). Stops if any name is not in the Catalog.
# Also writes, for each phenotype, every other Catalog trait name that contains a related word, so that missed
# spelling variants can be found and added to the term list by hand.
# Usage: Rscript scripts/07_external/01_catalog_terms.R [catalog_file]
#        (default: gwas_catalog in config/paths.yaml)
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
args <- commandArgs(trailingOnly = TRUE)
f <- if (length(args)) args[1] else pa$gwas_catalog
if (is.null(f) || !file.exists(f)) stop("GWAS Catalog file not found: give it as an argument or set gwas_catalog in config/paths.yaml")
out_dir <- file.path(pa$results_root, "07_external", ver, "catalog")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- read the Catalog (UCSC gwasCatalog layout as distributed with ANNOVAR, with or without a header) ----
has_header <- grepl("trait", readLines(f, n = 1), ignore.case = TRUE)
gc <- fread(f, sep = "\t", quote = "", fill = TRUE, header = has_header, showProgress = FALSE)
if (!has_header) {
  ucsc <- c("bin","chrom","chromStart","chromEnd","name","pubMedID","author","pubDate","journal","title",
            "trait","initSample","replSample","region","genes","riskAllele","riskAlFreq","pValue","pValueDesc",
            "orOrBeta","ci95","platform","cnv")
  k <- min(ncol(gc), length(ucsc)); setnames(gc, seq_len(k), ucsc[seq_len(k)])
}
pick <- function(re) { x <- grep(re, names(gc), ignore.case = TRUE, value = TRUE); if (!length(x)) stop("No column matching ", re); x[1] }
tcol <- pick("^trait$|^DISEASE/TRAIT$|^DISEASE\\.TRAIT$|^disease_trait$")
pcol <- pick("^pValue$|^P-VALUE$|^P\\.VALUE$|^p_value$")
ccol <- pick("^chrom$|^CHR_ID$|^chr$")
gc[, `:=`(trait_raw = trimws(get(tcol)), p = suppressWarnings(as.numeric(get(pcol))),
          chr = sub("^chr", "", as.character(get(ccol))))]
gc[, key := tolower(trait_raw)]
cat("Catalog file:", f, "| modified", format(file.mtime(f)), "| rows:", nrow(gc), "| distinct trait names:", uniqueN(gc$trait_raw), "\n")
if ("pubDate" %in% names(gc)) cat("Latest publication date in file:", max(as.character(gc$pubDate), na.rm = TRUE), "\n")
cat("\n")

by_name <- gc[, .(spelling_in_catalog = paste(sort(unique(trait_raw)), collapse = " | "),
                  n_assoc = .N, n_assoc_gw_autosomal = sum(p < 5e-8 & chr %in% as.character(1:22), na.rm = TRUE)), by = key]

# ---- the term list ----
terms <- fread(file.path(cfg$root, "config", "catalog_terms.tsv"), sep = "\t", header = TRUE)
terms[, key := tolower(trimws(catalog_trait))]
st2 <- merge(terms, by_name, by = "key", all.x = TRUE, sort = FALSE)
missing <- st2[is.na(n_assoc)]
if (nrow(missing)) { print(missing[, .(phenotype, catalog_trait)]); stop("The names above are not in the Catalog file") }
expected <- sub("_rint$", "", cfg$pipeline$analysis$phenotypes)
if (length(setdiff(expected, terms$phenotype))) stop("Phenotypes with no Catalog terms: ", paste(setdiff(expected, terms$phenotype), collapse = ", "))

st2[, phenotype := factor(phenotype, levels = unique(terms$phenotype))]
setorder(st2, phenotype, -n_assoc)
cat("== Table S3 terms: all", nrow(st2), "names found\n")
print(st2[, .(phenotype, catalog_trait = spelling_in_catalog, n_assoc, n_assoc_gw_autosomal)])
summ <- st2[, .(n_names = .N, n_assoc = sum(n_assoc), n_assoc_gw_autosomal = sum(n_assoc_gw_autosomal),
                catalog_traits = paste(spelling_in_catalog, collapse = "; ")), by = phenotype]
cat("\n== Per phenotype\n"); print(summ[, .(phenotype, n_names, n_assoc, n_assoc_gw_autosomal)])

# ---- names not in the list that contain a related word (for review by hand) ----
words <- c(bmi_calc = "body mass|\\bbmi\\b", bmi_at_25 = "body mass|\\bbmi\\b", weight_kg = "weight",
           standing_height_mm = "height|stature|body size", sitting_height_mm = "sitting|height",
           waist_mm = "waist", hip_mm = "hip circ|hip ratio|waist-hip|waist-to-hip", waist_hip_ratio = "waist-hip|waist-to-hip|waist to hip",
           fat_percent = "body fat|fat mass|fat percentage|adiposity", sbp_mean = "systolic|blood pressure", dbp_mean = "diastolic|blood pressure",
           heart_rate_mean = "heart rate|pulse rate|rr interval", random_glucose = "glucose|blood sugar")
rev <- rbindlist(lapply(names(words), function(ph) {
  x <- by_name[grepl(words[[ph]], key, perl = TRUE) & !key %in% terms[phenotype == ph, key]][order(-n_assoc)]
  if (!nrow(x)) return(NULL)
  x[, .(phenotype = ph, catalog_trait = spelling_in_catalog, n_assoc, n_assoc_gw_autosomal)] }))
cat("\nOther Catalog names with a related word (not in the list):", nrow(rev), "rows; see the review file.\n")

fwrite(st2[, .(phenotype, catalog_trait = spelling_in_catalog, n_assoc, n_assoc_gw_autosomal)], file.path(out_dir, "table_s3_catalog_terms_long.csv"))
fwrite(summ, file.path(out_dir, "table_s3_catalog_terms_by_phenotype.csv"))
fwrite(rev, file.path(out_dir, "table_s3_review_other_names.csv"))
cat("Saved in", out_dir, ":\n  table_s3_catalog_terms_long.csv\n  table_s3_catalog_terms_by_phenotype.csv\n  table_s3_review_other_names.csv\n")
