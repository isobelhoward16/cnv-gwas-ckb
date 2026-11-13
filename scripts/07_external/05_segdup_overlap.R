#!/usr/bin/env Rscript
# Stage 7, step 5: segmental-duplication (SD) context of each reported locus, i.e. whether the locus lies
# in a region flanked by SDs, as expected for recurrent CNVs mediated by non-allelic homologous recombination.
# Inputs: stage 5 signal table; callset restricted to participants passing CNV QC; UCSC genomicSuperDups (hg19).
# Per locus:
#   - CNVR (event boundaries where available, as in Table 1, otherwise the stage 5 CNVR; union over the
#     associations at the locus), its size, the fraction covered by SDs, and whether the lead probe lies in an SD;
#   - carrier calls covering the lead probe: number (del/dup), median start and end, and the fraction of calls whose
#     start and end both lie within 50 kb of the medians (a simple measure of recurrent breakpoints);
#   - flanking SD pairs: an SD whose two copies lie on the same chromosome, one within 250 kb of each breakpoint
#     (median call breakpoints), reported with identity, length and orientation (direct = same strand); pairs with
#     at least 95% identity and 10 kb length are also counted separately.
# Known recurrent loci (16p11.2 BP2-BP3 and BP4-BP5, 1q21.1 distal) act as positive controls.
# Output: results/07_external/<version>/segdups/segdup_by_locus.csv
source("R/config.R"); suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
options(width = 250)
FLANK <- 250e3; BP_TOL <- 50e3; MIN_ID <- 0.95; MIN_LEN <- 10e3
od <- file.path(pa$results_root, "07_external", ver, "segdups"); dir.create(od, recursive = TRUE, showWarnings = FALSE)
if (is.null(pa$segdups_hg19) || !file.exists(pa$segdups_hg19)) stop("Set paths$segdups_hg19 to UCSC hg19 genomicSuperDups.txt.gz")

# ---- SDs ----
sd <- fread(pa$segdups_hg19, header = FALSE, sep = "\t", select = c(2:4, 7:10, 27),
            col.names = c("chrom", "start", "end", "strand", "oChrom", "oStart", "oEnd", "fracMatch"))
sd <- sd[chrom %in% paste0("chr", 1:22)]
sd[, `:=`(chr = as.integer(sub("chr", "", chrom)), start = start + 1L, len = end - start + 1L)]
cat("SDs (autosomal, both copies listed):", nrow(sd), "| identity range", sprintf("%.3f-%.3f", min(sd$fracMatch), max(sd$fracMatch)), "\n")

# ---- loci from stage 5 (labels as in Table 1) ----
x <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x[, `:=`(s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end),
         ph = sub("_rint$", "", phenotype))]
# locus labels exactly as Table 1 (10_paper/01_table1.R): stage 5 "1q21.2" read as 1q21.1; non-overlapping CNVR
# clusters within one cytoband lettered a, b, ... by position
x[, band := sub("^1q21\\.2$", "1q21.1", locus_cytoband_final)]
x[, cl := { o <- order(s); grp <- integer(.N); g <- 0L; hi <- -Inf
  for (i in o) { if (s[i] > hi) { g <- g + 1L; hi <- e[i] } else hi <- max(hi, e[i]); grp[i] <- g }
  grp }, by = band]
x[, label := if (uniqueN(cl) > 1) paste0(band, letters[cl]) else band, by = band]
loci <- x[, .(chr = chr[1], cnvr_start = min(s), cnvr_end = max(e), lead_pos = lead_pos[which.min(pval)],
              phenotypes = paste(sort(unique(ph)), collapse = ", "), models = paste(sort(unique(cnv_type)), collapse = "+")), by = label]
setorder(loci, chr, cnvr_start)
cat("Loci:", nrow(loci), "\n")

# ---- calls in participants passing CNV QC ----
mat_ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))
qc_ids  <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
cs <- readRDS(file.path(pa$results_root, "01_callset", ver, "final", "callset.rds")); setDT(cs)
cs <- cs[CKID %in% intersect(mat_ids, qc_ids), .(chr = suppressWarnings(as.integer(sub("^chr", "", Chromosome))), start = StartPos, end = EndPos, type = CN_type)]
cat("Calls in participants passing CNV QC:", nrow(cs), "\n\n")

cover_frac <- function(ch, a, b) {           # fraction of [a, b] covered by the union of SDs
  y <- sd[chr == ch & end >= a & start <= b, .(s = pmax(start, a), e = pmin(end, b))]
  if (!nrow(y)) return(0)
  setorder(y, s); y[, g := cumsum(c(1L, s[-1] > cummax(e)[-.N]))]
  sum(y[, .(l = max(e) - min(s) + 1), by = g]$l) / (b - a + 1)
}
res <- rbindlist(lapply(seq_len(nrow(loci)), function(i) {
  L <- loci[i]
  cl <- cs[chr == L$chr & start <= L$lead_pos & end >= L$lead_pos]
  if ("mirror" %in% strsplit(L$models, "\\+")[[1]]) tp <- c("del", "dup") else tp <- sub("deletion", "del", sub("duplication", "dup", strsplit(L$models, "\\+")[[1]]))
  cl <- cl[type %in% tp]
  ms <- if (nrow(cl)) median(cl$start) else NA_real_; me <- if (nrow(cl)) median(cl$end) else NA_real_
  rec <- if (nrow(cl)) mean(abs(cl$start - ms) <= BP_TOL & abs(cl$end - me) <= BP_TOL) else NA_real_
  bs <- if (is.na(ms)) L$cnvr_start else ms; be <- if (is.na(me)) L$cnvr_end else me
  # SD pairs: one copy near each breakpoint, both copies on this chromosome
  pr <- sd[chr == L$chr & oChrom == paste0("chr", L$chr) & end >= bs - FLANK & start <= bs + FLANK &
           oEnd >= be - FLANK & oStart <= be + FLANK & oStart > end]
  best <- if (nrow(pr)) pr[order(-fracMatch, -len)][1] else NULL
  strong <- pr[fracMatch >= MIN_ID & len >= MIN_LEN]
  data.table(locus = L$label, chr = L$chr, phenotypes = L$phenotypes, models = L$models,
             CNVR_Mb = sprintf("%.2f-%.2f", L$cnvr_start / 1e6, L$cnvr_end / 1e6), CNVR_kb = round((L$cnvr_end - L$cnvr_start) / 1e3),
             SD_cover_CNVR_pct = round(100 * cover_frac(L$chr, L$cnvr_start, L$cnvr_end), 1),
             lead_in_SD = nrow(sd[chr == L$chr & start <= L$lead_pos & end >= L$lead_pos]) > 0,
             carrier_calls = nrow(cl), del = sum(cl$type == "del"), dup = sum(cl$type == "dup"),
             call_median_Mb = if (nrow(cl)) sprintf("%.2f-%.2f", ms / 1e6, me / 1e6) else "",
             recurrent_bp_pct = round(100 * rec, 0),
             SD_cover_250kb_start_pct = round(100 * cover_frac(L$chr, bs - FLANK, bs + FLANK), 1),
             SD_cover_250kb_end_pct = round(100 * cover_frac(L$chr, be - FLANK, be + FLANK), 1),
             flank_pairs = nrow(pr), flank_pairs_95pc_10kb = nrow(strong),
             flank_pairs_95pc_10kb_direct = nrow(strong[strand == "+"]),
             best_pair_identity = if (!is.null(best)) round(best$fracMatch, 4) else NA_real_,
             best_pair_len_kb = if (!is.null(best)) round(best$len / 1e3, 1) else NA_real_,
             best_pair_orientation = if (!is.null(best)) fifelse(best$strand == "+", "direct", "inverted") else "")
}))
cat("== Segmental-duplication context per locus\n")
print(res)
fwrite(res, file.path(od, "segdup_by_locus.csv"))
cat("\nSaved:", file.path(od, "segdup_by_locus.csv"), "\n")
cat("Notes: 'flank_pairs' = SD pairs with one copy within 250 kb of the median start and the partner within 250 kb of the\n",
    "median end of carrier calls; 'direct' = both copies on the same strand (UCSC strand '+'). Positive controls: 16p11.2a/b, 1q21.1.\n")
