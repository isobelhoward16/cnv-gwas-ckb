#!/usr/bin/env Rscript
# Stage 7, step 4: trait-matched GWAS Catalog SNP associations near each reported association.
# Autosomal Catalog associations with P < 5e-8 are matched to each trait by the exact Catalog trait names of
# Table S3 (config/catalog_terms.tsv; case-insensitive) and counted within the CNVR, the CNVR +/-250 kb, +/-500 kb
# and +/-1 Mb, and the union of the CNVRs of all models significant at the locus ("Auwerx" column). Counts are unique rsIDs, approximately independent loci (rsIDs
# more than 250 kb apart start a new locus) and publications (PMIDs). The CNVR is that of the model shown in
# Table 1 (the mirror CNVR when both a single-type and the mirror model are significant).
# Also writes every hit with its distance to the CNVR, PMID, first author and year.
# Output: results/07_external/<version>/catalog/
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
CLUMP_KB <- 250
FLANKS <- c(CNVR = 0, CNVR_250kb = 2.5e5, CNVR_500kb = 5e5, CNVR_1Mb = 1e6)
od <- file.path(pa$results_root, "07_external", ver, "catalog"); dir.create(od, recursive = TRUE, showWarnings = FALSE)

# ---- Catalog (read as in 01_catalog_terms.R) ----
f <- pa$gwas_catalog
has_header <- grepl("trait", readLines(f, n = 1), ignore.case = TRUE)
gc <- fread(f, sep = "\t", quote = "", fill = TRUE, header = has_header, showProgress = FALSE)
if (!has_header) {
  ucsc <- c("bin","chrom","chromStart","chromEnd","name","pubMedID","author","pubDate","journal","title","trait","initSample",
            "replSample","region","genes","riskAllele","riskAlFreq","pValue","pValueDesc","orOrBeta","ci95","platform","cnv")
  k <- min(ncol(gc), length(ucsc)); setnames(gc, seq_len(k), ucsc[seq_len(k)])
}
gc <- gc[, .(CHR_NUM = suppressWarnings(as.integer(sub("^chr", "", chrom))), BP = as.integer(floor((as.numeric(chromStart) + as.numeric(chromEnd)) / 2)),
             RSID = name, PMID = pubMedID, FIRST_AUTHOR = author, PUB_DATE = pubDate, TRAIT = trimws(trait), PVALUE = suppressWarnings(as.numeric(pValue)))]
cat("Catalog:", f, "| modified", format(file.mtime(f)), "| rows", nrow(gc), "\n")
gc <- gc[!is.na(PVALUE) & PVALUE < 5e-8 & CHR_NUM %in% 1:22]
cat("Autosomal entries with P < 5e-8:", nrow(gc), "\n")
terms <- fread(file.path(cfg$root, "config", "catalog_terms.tsv"), sep = "\t", header = TRUE)
terms[, key := tolower(trimws(catalog_trait))]
gc[, key := tolower(TRAIT)]

# ---- locus-phenotype pairs from stage 5, labelled as in Table 1 ----
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x[, `:=`(s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end),
         ph = sub("_rint$", "", phenotype))]
x[, label := locus_cytoband_final]
x[label == "16p11.2", label := fifelse(s < 29.5e6, "16p11.2a", "16p11.2b")]
x[label %in% c("1q21.1", "1q21.2"), label := fifelse(lead_pos < 145e6, "1q21.1a", "1q21.1b")]
pairs <- x[, {
  shown <- if ("mirror" %in% cnv_type) which(cnv_type == "mirror")[1] else which.min(pval)
  .(chr = chr[1], models = paste(sort(unique(cnv_type)), collapse = "+"), lead_pos = lead_pos[which.min(pval)],
    cnvr_start = s[shown], cnvr_end = e[shown], auwerx_start = min(s), auwerx_end = max(e))
}, by = .(label, ph)]
setorder(pairs, chr, cnvr_start, ph)
cat("Locus-phenotype pairs:", nrow(pairs), "(Table 1 has 14)\n\n")

clump <- function(bp) { if (!length(bp)) return(0L); bp <- sort(unique(bp)); 1L + sum(diff(bp) > CLUMP_KB * 1e3) }
summ <- list(); hits_out <- list()
for (i in seq_len(nrow(pairs))) {
  r <- pairs[i]
  tk <- terms[phenotype == r$ph, key]
  if (!length(tk)) stop("No Table S3 terms for ", r$ph)
  hit_chr <- gc[key %in% tk & CHR_NUM == r$chr]
  res <- data.table(label = r$label, phenotype = r$ph, models = r$models, chr = r$chr, lead_pos = r$lead_pos,
                    cnvr_start = r$cnvr_start, cnvr_end = r$cnvr_end, cnvr_kb = round((r$cnvr_end - r$cnvr_start) / 1e3))
  wins <- c(lapply(FLANKS, function(fl) c(r$cnvr_start - fl, r$cnvr_end + fl)), list(Auwerx = c(r$auwerx_start, r$auwerx_end)))
  for (w in names(wins)) {
    h <- hit_chr[BP >= wins[[w]][1] & BP <= wins[[w]][2]]
    res[, paste0(c("n_rsid_", "n_loci_", "n_pmid_"), w) := list(uniqueN(h$RSID), clump(h$BP), uniqueN(h$PMID))]
  }
  wmax <- wins[["CNVR_1Mb"]]
  res[, n_rsid_alltraits_CNVR_1Mb := uniqueN(gc[CHR_NUM == r$chr & BP >= wmax[1] & BP <= wmax[2], RSID])]
  summ[[i]] <- res
  h <- hit_chr[BP >= wmax[1] & BP <= wmax[2]]
  if (nrow(h)) {
    h[, dist_to_cnvr_kb := round(pmax(0, r$cnvr_start - BP, BP - r$cnvr_end) / 1e3)]
    h[, window := fifelse(dist_to_cnvr_kb == 0, "CNVR", fifelse(dist_to_cnvr_kb <= 250, "CNVR_250kb", fifelse(dist_to_cnvr_kb <= 500, "CNVR_500kb", "CNVR_1Mb")))]
    hits_out[[i]] <- h[, .(label = r$label, phenotype = r$ph, window, dist_to_cnvr_kb, RSID, chr = CHR_NUM, BP, TRAIT, PVALUE, PMID, FIRST_AUTHOR, PUB_DATE)][order(dist_to_cnvr_kb, PVALUE)]
  }
}
summ <- rbindlist(summ)
cat("== Trait-matched GWAS Catalog associations (unique rsIDs / approx. independent loci / PMIDs)\n")
print(summ[, .(label, phenotype, models, CNVR = sprintf("%.2f-%.2f Mb", cnvr_start / 1e6, cnvr_end / 1e6), cnvr_kb,
               CNVR_ = sprintf("%d/%d/%d", n_rsid_CNVR, n_loci_CNVR, n_pmid_CNVR),
               pm250kb = sprintf("%d/%d/%d", n_rsid_CNVR_250kb, n_loci_CNVR_250kb, n_pmid_CNVR_250kb),
               pm500kb = sprintf("%d/%d/%d", n_rsid_CNVR_500kb, n_loci_CNVR_500kb, n_pmid_CNVR_500kb),
               pm1Mb = sprintf("%d/%d/%d", n_rsid_CNVR_1Mb, n_loci_CNVR_1Mb, n_pmid_CNVR_1Mb),
               Auwerx = sprintf("%d/%d/%d", n_rsid_Auwerx, n_loci_Auwerx, n_pmid_Auwerx), alltraits_1Mb = n_rsid_alltraits_CNVR_1Mb)])
fwrite(summ, file.path(od, "catalog_counts_summary.csv"))
if (length(hits_out)) fwrite(rbindlist(hits_out), file.path(od, "catalog_counts_hits.csv"))
cat("\nSaved:", file.path(od, "catalog_counts_summary.csv"), "and catalog_counts_hits.csv\n")
