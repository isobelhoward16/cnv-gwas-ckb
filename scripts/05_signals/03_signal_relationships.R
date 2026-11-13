#!/usr/bin/env Rscript
# Stage 5, step 3: relationships between signals. For every pair of signals on the same chromosome
# whose locus boundaries overlap, reports the probe spans, lead probes, P values and betas, carrier
# counts and the overlap of carriers (Jaccard index), and classifies the pair as:
#   (A) same signal under different models (same trait, overlapping locus);
#   (B) same locus, different trait (overlapping locus and Jaccard > 0.5);
#   (C) overlapping locus but likely distinct signals (Jaccard <= 0.5).
# Signals are then grouped into loci.
# Usage: Rscript scripts/05_signals/03_signal_relationships.R [all|unrelated]
suppressPackageStartupMessages({ library(data.table); library(Matrix) })

# ============================================================
# PATHS
# ============================================================
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
args <- commandArgs(trailingOnly = TRUE)
rel_arg <- if (length(args)) args[1] else "all"
relatedness_group <- if (rel_arg == "all") "related" else "unrelated"
cnvr_path <- file.path(pa$results_root, "05_signals", ver, rel_arg, "cnvr_outputs",
                       paste0("cnvr_by_model_supplemental_unpruned_", relatedness_group, ".csv"))
mat_base <- file.path(pa$results_root, "02_matrices", ver)
probe_map_path <- file.path(mat_base, "combined_probe_map_shared_probes.rds")

# ============================================================
# LOAD DATA
# ============================================================
cat("Loading CNVR table...\n")
cnvr <- fread(cnvr_path)
cnvr[, sig_id := .I]
cnvr[, label  := paste0(phenotype, " / ", cnv_type, " / ", locus_cytoband_final)]
cat("  ->", nrow(cnvr), "signals\n\n")

cat("Loading probe map...\n")
pm <- readRDS(probe_map_path)
setDT(pm)
pm[, chr := as.integer(chr)]
pm[, pos := as.integer(pos)]
cat("  ->", nrow(pm), "probes\n\n")

# Load deletion matrix for exact carrier overlap
# (deletion matrix covers both del and dup carriers via dosage encoding)
cat("Loading deletion matrix for carrier ID extraction...\n")
del_mat <- readRDS(file.path(mat_base, "combined_deletion_matrix_shared_probes_fixed.rds"))
dup_mat <- readRDS(file.path(mat_base, "combined_duplication_matrix_shared_probes_fixed.rds"))
cat("  -> matrices loaded\n\n")

# ============================================================
# HELPER: get carrier IDs at a locus from a matrix
# A carrier at a locus = sample with dosage != 0 at ANY probe
# in the Auwerx CNV region
# ============================================================
get_carriers_from_matrix <- function(mat, probe_map, chr_val, region_start, region_end) {
  # Find probes in this region
  probes_in_region <- probe_map[chr == chr_val &
                                 pos >= region_start &
                                 pos <= region_end, probe_id]
  probes_in_region <- intersect(probes_in_region, colnames(mat))
  if (!length(probes_in_region)) return(character(0))

  # Subset matrix to those probes
  sub <- mat[, probes_in_region, drop = FALSE]

  # Carriers = rows with any non-zero entry
  # For sparse matrices
  if (inherits(sub, "sparseMatrix")) {
    carrier_idx <- unique(sub@i + 1L)  # 0-indexed -> 1-indexed
  } else {
    carrier_idx <- which(rowSums(sub != 0) > 0)
  }
  rownames(mat)[carrier_idx]
}

# ============================================================
# HELPER: exact Jaccard from two character vectors
# ============================================================
jaccard_exact <- function(a, b) {
  if (!length(a) && !length(b)) return(NA_real_)
  inter <- length(intersect(a, b))
  union <- length(union(a, b))
  list(jaccard = round(inter / union, 3),
       intersection = inter,
       union = union,
       n_a = length(a),
       n_b = length(b))
}

# ============================================================
# HELPER: span overlap in kb
# ============================================================
overlap_kb <- function(s1, e1, s2, e2) {
  round(max(0, min(e1, e2) - max(s1, s2)) / 1e3, 1)
}

# ============================================================
# PRE-COMPUTE CARRIER SETS FOR EACH SIGNAL
# (do this once - expensive to redo per pair)
# ============================================================
cat("Extracting carrier sets for all 27 signals...\n")
carrier_sets <- vector("list", nrow(cnvr))

for (i in seq_len(nrow(cnvr))) {
  r <- cnvr[i]
  mat <- if (r$cnv_type == "duplication") dup_mat else del_mat

  carriers <- get_carriers_from_matrix(
    mat, pm,
    chr_val      = as.integer(r$chr),
    region_start = as.integer(r$auwerx_event_start),
    region_end   = as.integer(r$auwerx_event_end)
  )
  carrier_sets[[i]] <- carriers
  cat(sprintf("  [%d] %s -> %d carriers\n", i, r$label, length(carriers)))
}
cat("\n")

# ============================================================
# PAIRWISE COMPARISONS - only for pairs on same chromosome
# ============================================================
n <- nrow(cnvr)
results <- vector("list", n * (n-1) / 2)
k <- 0L

for (i in seq_len(n - 1)) {
  for (j in (i + 1):n) {
    a <- cnvr[i]; b <- cnvr[j]
    if (a$chr != b$chr) next  # different chromosomes - skip

    # Geometric checks
    aw_overlap    <- overlap_kb(a$auwerx_event_start, a$auwerx_event_end,
                                b$auwerx_event_start, b$auwerx_event_end)
    sig_overlap   <- overlap_kb(a$event_start, a$event_end,
                                b$event_start, b$event_end)
    lead_dist_kb  <- round(abs(a$lead_pos - b$lead_pos) / 1e3, 1)

    if (aw_overlap == 0) next  # no overlap at all - skip

    # Exact carrier overlap
    jac <- jaccard_exact(carrier_sets[[i]], carrier_sets[[j]])

    # Relationship classification
    same_pheno <- a$phenotype == b$phenotype
    same_model <- a$cnv_type  == b$cnv_type

    rel_type <- if (same_pheno && !same_model) {
      "A_same_signal_diff_model"
    } else if (!same_pheno && jac$jaccard > 0.5) {
      "B_same_locus_diff_phenotype"
    } else if (!same_pheno && jac$jaccard <= 0.5) {
      "C_overlapping_low_carrier_overlap"
    } else {
      "D_same_pheno_same_model_overlapping"
    }

    k <- k + 1L
    results[[k]] <- data.table(
      rel_type      = rel_type,
      sig_i = i, sig_j = j,
      label_i = a$label, label_j = b$label,
      # Spans
      sig_start_i   = a$event_start,      sig_end_i   = a$event_end,
      sig_start_j   = b$event_start,      sig_end_j   = b$event_end,
      sig_span_kb_i = round((a$event_end - a$event_start)/1e3, 1),
      sig_span_kb_j = round((b$event_end - b$event_start)/1e3, 1),
      sig_overlap_kb = sig_overlap,
      aw_start_i    = a$auwerx_event_start, aw_end_i  = a$auwerx_event_end,
      aw_start_j    = b$auwerx_event_start, aw_end_j  = b$auwerx_event_end,
      aw_span_kb_i  = round((a$auwerx_event_end - a$auwerx_event_start)/1e3, 1),
      aw_span_kb_j  = round((b$auwerx_event_end - b$auwerx_event_start)/1e3, 1),
      aw_overlap_kb = aw_overlap,
      # Lead probes
      lead_pos_i    = a$lead_pos, lead_pos_j = b$lead_pos,
      lead_dist_kb  = lead_dist_kb,
      lead_same_side = (a$lead_pos >= b$auwerx_event_start &
                        a$lead_pos <= b$auwerx_event_end),
      # Stats
      pval_i   = a$pval,        pval_j   = b$pval,
      beta_i   = round(a$lead_beta, 3), beta_j = round(b$lead_beta, 3),
      same_direction = sign(a$lead_beta) == sign(b$lead_beta),
      # Carriers
      n_del_i  = a$n_del,  n_dup_i  = a$n_dup,
      n_del_j  = b$n_del,  n_dup_j  = b$n_dup,
      n_carrier_i = a$n_carrier, n_carrier_j = b$n_carrier,
      carrier_jaccard    = jac$jaccard,
      carrier_intersect  = jac$intersection,
      carrier_union      = jac$union
    )
  }
}

dt <- rbindlist(results[seq_len(k)], fill = TRUE)
setorder(dt, rel_type, sig_i, sig_j)

# ============================================================
# PRINT FULL DIAGNOSTICS
# ============================================================

print_pair <- function(r) {
  cat(sprintf(
"  Signal i: [%d] %s
           sig span: %s-%s (%g kb) | Auwerx: %s-%s (%g kb) | lead: %s | p=%s | beta=%g
           n_del=%d  n_dup=%d  n_carrier=%d
  Signal j: [%d] %s
           sig span: %s-%s (%g kb) | Auwerx: %s-%s (%g kb) | lead: %s | p=%s | beta=%g
           n_del=%d  n_dup=%d  n_carrier=%d
  Overlap:  Auwerx=%.1fkb | SigProbe=%.1fkb | Lead dist=%.1fkb | Same direction=%s
  Carriers: Jaccard=%.3f | Intersection=%d | Union=%d (n_i=%d, n_j=%d)\n\n",
    r$sig_i, r$label_i,
    format(r$sig_start_i, big.mark=",", scientific=FALSE),
    format(r$sig_end_i,   big.mark=",", scientific=FALSE),
    r$sig_span_kb_i,
    format(r$aw_start_i,  big.mark=",", scientific=FALSE),
    format(r$aw_end_i,    big.mark=",", scientific=FALSE),
    r$aw_span_kb_i,
    format(r$lead_pos_i,  big.mark=",", scientific=FALSE),
    formatC(r$pval_i, format="e", digits=2),
    r$beta_i,
    r$n_del_i, r$n_dup_i, r$n_carrier_i,
    r$sig_j, r$label_j,
    format(r$sig_start_j, big.mark=",", scientific=FALSE),
    format(r$sig_end_j,   big.mark=",", scientific=FALSE),
    r$sig_span_kb_j,
    format(r$aw_start_j,  big.mark=",", scientific=FALSE),
    format(r$aw_end_j,    big.mark=",", scientific=FALSE),
    r$aw_span_kb_j,
    format(r$lead_pos_j,  big.mark=",", scientific=FALSE),
    formatC(r$pval_j, format="e", digits=2),
    r$beta_j,
    r$n_del_j, r$n_dup_j, r$n_carrier_j,
    r$aw_overlap_kb, r$sig_overlap_kb, r$lead_dist_kb,
    r$same_direction,
    r$carrier_jaccard, r$carrier_intersect, r$carrier_union,
    r$n_carrier_i, r$n_carrier_j
  ))
}

sep <- paste(rep("-", 100), collapse="")
SEP <- paste(rep("=", 100), collapse="")

cat(SEP, "\n")
cat("  (A) SAME SIGNAL - different models (same phenotype, overlapping locus)\n")
cat(SEP, "\n\n")
grp_a <- dt[rel_type == "A_same_signal_diff_model"]
if (nrow(grp_a)) for (i in seq_len(nrow(grp_a))) { cat(sep,"\n"); print_pair(grp_a[i]) }
cat(sprintf("  -> %d pairs\n\n", nrow(grp_a)))

cat(SEP, "\n")
cat("  (B) SAME LOCUS - different phenotypes (Auwerx overlap >0, Jaccard >0.5)\n")
cat(SEP, "\n\n")
grp_b <- dt[rel_type == "B_same_locus_diff_phenotype"]
if (nrow(grp_b)) for (i in seq_len(nrow(grp_b))) { cat(sep,"\n"); print_pair(grp_b[i]) }
cat(sprintf("  -> %d pairs\n\n", nrow(grp_b)))

cat(SEP, "\n")
cat("  (C) OVERLAPPING AUWERX but LOW carrier overlap (Jaccard <=0.5) - likely distinct\n")
cat(SEP, "\n\n")
grp_c <- dt[rel_type == "C_overlapping_low_carrier_overlap"]
if (nrow(grp_c)) for (i in seq_len(nrow(grp_c))) { cat(sep,"\n"); print_pair(grp_c[i]) }
cat(sprintf("  -> %d pairs\n\n", nrow(grp_c)))

cat(SEP, "\n")
cat("  (D) SAME PHENOTYPE + SAME MODEL + OVERLAPPING - check for duplicates\n")
cat(SEP, "\n\n")
grp_d <- dt[rel_type == "D_same_pheno_same_model_overlapping"]
if (nrow(grp_d)) for (i in seq_len(nrow(grp_d))) { cat(sep,"\n"); print_pair(grp_d[i]) }
cat(sprintf("  -> %d pairs\n\n", nrow(grp_d)))

# ============================================================
# UNION-FIND LOCUS GROUPING using exact Jaccard
# ============================================================
cat(SEP, "\n")
cat("  LOCUS GROUP ASSIGNMENTS (collapsing A + B)\n")
cat(SEP, "\n\n")

parent <- seq_len(n)
find  <- function(x) { while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }; x }
unite <- function(x, y) { px <- find(x); py <- find(y); if (px != py) parent[px] <<- py }

for (i in seq_len(nrow(dt))) {
  r <- dt[i]
  if (r$rel_type %in% c("A_same_signal_diff_model", "B_same_locus_diff_phenotype"))
    unite(r$sig_i, r$sig_j)
}

groups <- data.table(sig_id = seq_len(n), label = cnvr$label,
                     group  = sapply(seq_len(n), find))
setorder(groups, group, sig_id)

for (g in unique(groups$group)) {
  m <- groups[group == g]
  if (nrow(m) > 1) {
    cat(sprintf("  *** Group %d (%d signals - same locus):\n", g, nrow(m)))
    for (r in seq_len(nrow(m))) cat(sprintf("      [%d] %s\n", m$sig_id[r], m$label[r]))
  } else {
    cat(sprintf("  Group %d: [%d] %s\n", g, m$sig_id, m$label))
  }
  cat("\n")
}

cat(sprintf("\n  Total signals: %d  |  Unique loci after collapsing: %d\n\n",
            n, uniqueN(sapply(seq_len(n), find))))
