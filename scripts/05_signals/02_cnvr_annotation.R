#!/usr/bin/env Rscript
# Stage 5, step 2: CNV regions (CNVRs) and their annotation.
# For each trait and model, the most significant remaining probe is taken as a lead, and significant
# probes within +/-3 Mb whose CNV dosage has r2 >= 0.5 with the lead form one CNVR; this is repeated
# until every significant probe is assigned. Locus boundaries follow Auwerx et al. (2022): the most
# distal probes within +/-3 Mb with r2 >= 0.5 with the lead, using the regression sample of the trait.
# Adds carrier counts at the lead probe (from the matrix, as in the regression), overlapping genes
# (GENCODE v19), cytobands and the experiment-wide flag, and merges CNVRs across models.
# Output: results/05_signals/<version>/<all|unrelated>/cnvr_outputs/
#         cnvr_by_model_supplemental_unpruned_*.csv, cnvr_merged_loci_unpruned_*.csv,
#         cnvr_probe_membership_unpruned_*.csv
# Usage: Rscript scripts/05_signals/02_cnvr_annotation.R <all|unrelated>
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(GenomicRanges)
  library(IRanges)
  library(rtracklayer)
})

cat("=== UNPRUNED: CNVR + SUPPLEMENTAL + CYTOBAND ===\n")
cat("Working directory:", getwd(), "\n\n")

source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; an <- cfg$pipeline$analysis; ver <- cfg$pipeline$callset_version
args <- commandArgs(trailingOnly = TRUE)
rel_arg <- if (length(args)) args[1] else "all"
stopifnot(rel_arg %in% c("all", "unrelated"))
relatedness_group <- if (rel_arg == "all") "related" else "unrelated"   # file-name label for each sample set

sig_dir <- file.path(pa$results_root, "05_signals", ver, rel_arg)
sig_probes_file <- file.path(sig_dir, paste0("significant_probes_unpruned_", relatedness_group, ".csv"))

mat_base_dir   <- file.path(pa$results_root, "02_matrices", ver)
probe_map_file <- file.path(mat_base_dir, "combined_probe_map_shared_probes.rds")

gwas_dir <- file.path(pa$results_root, "04_gwas", ver, rel_arg)

cnv_dir   <- file.path(pa$results_root, "01_callset", ver, "final", "calls")
pheno_rds <- pa$phenotype_rds
key1_path <- pa$array1_key
key2_path <- pa$array2_key
gtf_path  <- pa$gtf_gencode_v19
cyto_file <- pa$cytoband_hg19

tag_kb    <- 3000
tag_r2    <- 0.5
min_pairs <- 10

use_rint           <- FALSE
del_types          <- c("del","DEL","deletion")
dup_types          <- c("dup","DUP","duplication")
mirror_matches_dup <- TRUE
n_traits_tested    <- an$n_effective_traits   # effective number of traits (experiment-wide threshold)

# =============================================================================
# meff_map: maximum Meff across traits per model (stage 3 thresholds table)
# =============================================================================

thr_tab <- fread(file.path(pa$results_root, "03_meff", ver, "meff_thresholds.csv"))
meff_map <- setNames(thr_tab[relatedness == rel_arg, Meff], thr_tab[relatedness == rel_arg, model])
meff_map <- meff_map[c("deletion", "duplication", "mirror")]

cat("meff_map (maximum Meff across traits):\n")
print(meff_map)
cat("\n")

# =============================================================================
# SAFETY CHECKS
# =============================================================================

for (p in c(sig_probes_file, probe_map_file, pheno_rds, key1_path, key2_path, gtf_path, cyto_file)) {
  if (!file.exists(p)) stop("Missing file: ", p)
}
if (!dir.exists(cnv_dir))      stop("Missing cnv_dir: ", cnv_dir)
if (!dir.exists(mat_base_dir)) stop("Missing mat_base_dir: ", mat_base_dir)
if (!dir.exists(gwas_dir))     stop("Missing gwas_dir: ", gwas_dir)

# =============================================================================
# HELPERS
# =============================================================================

find_matrix_file <- function(cnv_type, mat_base_dir) {
  f <- file.path(mat_base_dir, paste0("combined_", cnv_type, "_matrix_shared_probes_fixed.rds"))
  if (file.exists(f)) return(f)
  patt  <- paste0("^combined_", cnv_type, "_matrix_shared_probes.*\\.rds$")
  files <- list.files(mat_base_dir, pattern = patt, full.names = TRUE)
  files <- files[!grepl("pruned", files)]
  if (length(files) == 0) stop("No unpruned matrix found for cnv_type=", cnv_type)
  files[which.max(file.info(files)$mtime)]
}

safe_cor <- function(x, y, min_pairs = 10) {
  ok <- !is.na(x) & !is.na(y)
  if (sum(ok) < min_pairs) return(NA_real_)
  suppressWarnings(cor(x[ok], y[ok]))
}

get_sign <- function(x) {
  if (is.na(x)) return(NA_integer_)
  if (x > 0) return(1L)
  if (x < 0) return(-1L)
  0L
}

intervals_overlap <- function(a_start, a_end, b_start, b_end) {
  !(a_end < b_start | b_end < a_start)
}

directionally_concordant <- function(model_a, beta_a, model_b, beta_b,
                                     mirror_matches_dup = TRUE) {
  sa <- get_sign(beta_a); sb <- get_sign(beta_b)
  if (is.na(sa) || is.na(sb) || sa == 0 || sb == 0) return(FALSE)
  if (model_a == model_b) return(sa == sb)
  if (!("mirror" %in% c(model_a, model_b))) return(FALSE)
  if (model_a == "mirror") { s_m <- sa; other <- model_b; s_o <- sb
  } else { s_m <- sb; other <- model_a; s_o <- sa }
  if (mirror_matches_dup) {
    if (other == "duplication") return(s_m == s_o)
    if (other == "deletion")    return(s_m == -s_o)
  } else {
    if (other == "deletion")    return(s_m == s_o)
    if (other == "duplication") return(s_m == -s_o)
  }
  FALSE
}

extract_sample_ids <- function(obj, cnv_type) {
  rn <- tryCatch(rownames(obj), error=function(e) NULL)
  if (!is.null(rn) && length(rn) > 0) return(rn)
  dn <- tryCatch(dimnames(obj), error=function(e) NULL)
  if (!is.null(dn) && length(dn) >= 1 && !is.null(dn[[1]]) && length(dn[[1]]) > 0) return(dn[[1]])
  if (isS4(obj)) {
    dn2 <- tryCatch(obj@Dimnames[[1]], error=function(e) NULL)
    if (!is.null(dn2) && length(dn2) > 0) return(dn2)
  }
  stop("Cannot extract sample IDs for cnv_type=", cnv_type)
}

get_pheno_col <- function(ph_all, pheno_base, use_rint = FALSE) {
  if (use_rint) {
    col <- paste0(pheno_base, "_rint")
    if (!col %in% names(ph_all)) stop("Missing phenotype column: ", col)
    return(col)
  }
  if (pheno_base %in% names(ph_all)) return(pheno_base)
  col <- paste0(pheno_base, "_rint")
  if (col %in% names(ph_all)) return(col)
  stop("Missing phenotype columns: ", pheno_base, " and ", col)
}

safe_stats <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(list(mean=NA_real_, median=NA_real_, sd=NA_real_, n=0L))
  list(mean=mean(x), median=as.numeric(stats::median(x)), sd=stats::sd(x), n=length(x))
}

collapse_by_row <- function(hit_df, n_rows) {
  if (nrow(hit_df) == 0L) return(rep(NA_character_, n_rows))
  hit_df <- hit_df[!is.na(hit_df$gene) & hit_df$gene != "", , drop=FALSE]
  if (nrow(hit_df) == 0L) return(rep(NA_character_, n_rows))
  spl <- split(hit_df$gene, hit_df$query_idx)
  out <- rep(NA_character_, n_rows)
  idx <- as.integer(names(spl))
  out[idx] <- vapply(spl, function(x) paste(sort(unique(x)), collapse=";"), character(1))
  out
}

collapse_unique <- function(x) {
  x <- unique(x[!is.na(x) & x != ""])
  if (!length(x)) return(NA_character_)
  paste(sort(x), collapse=";")
}

genes_overlapping_interval_scalar <- function(chr_int, start, end, genes_gr, gene_names) {
  if (is.na(chr_int) || is.na(start) || is.na(end)) return(character(0))
  if (!is.finite(start) || !is.finite(end)) return(character(0))
  if (end < start) { tmp <- start; start <- end; end <- tmp }
  gr <- GRanges(seqnames=paste0("chr", chr_int), ranges=IRanges(start=start, end=end))
  ov <- findOverlaps(gr, genes_gr, ignore.strand=TRUE)
  s  <- subjectHits(ov)
  if (!length(s)) return(character(0))
  gene_names[s]
}

# =============================================================================
# GWAS FILE INDEX + n_analyzed LOOKUP
# =============================================================================

build_gwas_file_index <- function(gwas_dir) {
  files <- list.files(gwas_dir,
                      pattern = "^gwas_linear_unpruned_combined_.*\\.csv$",
                      full.names = TRUE)
  if (length(files) == 0) stop("No GWAS CSV files found in: ", gwas_dir)
  idx <- data.table(path = files)
  idx[, bn := basename(path)]
  idx[, cnv_type  := fifelse(grepl("_combined_deletion_",    bn), "deletion",
                     fifelse(grepl("_combined_duplication_", bn), "duplication",
                     fifelse(grepl("_combined_mirror_",      bn), "mirror", NA_character_)))]
  idx[, phenotype := sub(paste0("^gwas_linear_unpruned_combined_(?:deletion|duplication|mirror)_(.+?)_mincar[0-9]+_.*$"), "\\1", bn)]
  idx[, mtime := file.info(path)$mtime]
  idx <- idx[idx[, .I[which.max(mtime)], by = .(cnv_type, phenotype)]$V1]
  idx[, .(cnv_type, phenotype, path)]
}

gwas_index   <- build_gwas_file_index(gwas_dir)
gwas_n_cache <- list()

get_n_analyzed <- function(cnv_type_val, phenotype_val, lead_probe_val) {
  key <- paste0(cnv_type_val, "___", phenotype_val)
  if (is.null(gwas_n_cache[[key]])) {
    row <- gwas_index[cnv_type == cnv_type_val & phenotype == phenotype_val]
    if (nrow(row) == 0) {
      warning("No GWAS CSV found for cnv_type=", cnv_type_val, " phenotype=", phenotype_val)
      return(NA_integer_)
    }
    gwas_n_cache[[key]] <<- fread(row$path[1], select = c("probe_id", "n_analyzed"))
  }
  dt <- gwas_n_cache[[key]]
  n  <- dt[probe_id == lead_probe_val, n_analyzed]
  if (length(n) == 0 || is.na(n[1])) {
    warning("Lead probe ", lead_probe_val, " not found in GWAS CSV for ",
            cnv_type_val, "/", phenotype_val)
    return(NA_integer_)
  }
  as.integer(n[1])
}

cat("GWAS file index built:", nrow(gwas_index), "files indexed\n\n")

# =============================================================================
# LOAD INPUTS
# =============================================================================

sig <- fread(sig_probes_file)
probe_map <- readRDS(probe_map_file); setDT(probe_map)
probe_map[, chr := gsub("^chr", "", as.character(chr))]
probe_map[, chr_num := suppressWarnings(as.integer(chr))]

if (!all(c("chr","pos") %in% names(sig))) {
  sig <- merge(sig, probe_map[, .(probe_id, chr, pos)], by="probe_id", all.x=TRUE)
}
sig <- sig[!is.na(chr) & !is.na(pos)]
sig[, chr := gsub("^chr","", as.character(chr))]
sig[, chr_num := suppressWarnings(as.integer(chr))]
if (!("beta" %in% names(sig))) sig[, beta := NA_real_]
if (!("se"   %in% names(sig))) sig[, se   := NA_real_]

sig <- sig[cnv_type %in% c("deletion","duplication","mirror")]
cat("Loaded significant probes:", nrow(sig), "\n")
cat("  Phenotypes:", paste(unique(sig$phenotype), collapse=", "), "\n")
cat("  Models:", paste(unique(sig$cnv_type), collapse=", "), "\n\n")

# =============================================================================
# LOAD PHENOTYPES + KEYS
# =============================================================================

ph <- readRDS(pheno_rds); setDT(ph)
if (relatedness_group == "unrelated") {
  if ("unrelated" %in% names(ph)) ph <- ph[unrelated == 1]
}

k1 <- fread(key1_path); k2 <- fread(key2_path)
setnames(k1, c("ccvid","xIID_cnv"))
setnames(k2, c("ccvid","xIID_cnv"))

ph1 <- ph[gwas_array_type == "Axiom_CKB1"]; ph2 <- ph[gwas_array_type == "Axiom_CKB2"]
ph1 <- merge(ph1, k1, by="ccvid", all=FALSE)
ph2 <- merge(ph2, k2, by="ccvid", all=FALSE)

ph_all <- rbindlist(list(ph1, ph2), use.names=TRUE, fill=TRUE)
setkey(ph_all, xIID_cnv)

# Master-QC filter, identical to scripts/04_gwas/01_gwas_task.R and scripts/03_meff/01_meff_task.R
qc_master <- fread(pa$qc_gwas_file)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_samples <- qc_master[qc_pass == TRUE]$GID
cat("Master-QC passing samples:", length(qc_pass_samples), "\n")
cat("Loaded phenotype rows with xIID_cnv:", nrow(ph_all), "\n\n")

# =============================================================================
# CACHE MATRIX SAMPLE SET + N_TESTED_PROBES PER MODEL
# =============================================================================

models_present <- sort(unique(sig$cnv_type))
mat_row_cache  <- list()
tested_probe_counts <- data.table(cnv_type=models_present, n_tested_probes=NA_integer_)

for (ct in models_present) {
  f <- find_matrix_file(ct, mat_base_dir)
  cat("Caching matrix info:", ct, "->", basename(f), "\n")
  obj <- readRDS(f)
  mat_row_cache[[ct]] <- extract_sample_ids(obj, ct)
  tested_probe_counts[cnv_type == ct, n_tested_probes := ncol(obj)]
  rm(obj); gc()
}
cat("\n")

get_gwas_samples <- function(pheno_base, cnv_type) {
  col      <- get_pheno_col(ph_all, pheno_base, use_rint=use_rint)
  ph_sub   <- ph_all[!is.na(get(col)), xIID_cnv]
  mat_rows <- mat_row_cache[[cnv_type]]
  if (is.null(mat_rows)) stop("No cached matrix rownames for cnv_type: ", cnv_type)
  intersect(intersect(ph_sub, mat_rows), qc_pass_samples)
}

# =============================================================================
# LOAD GENE MODELS
# =============================================================================

cat("Loading gene annotation:\n  ", gtf_path, "\n")
gtf      <- import(gtf_path)
genes_gr <- gtf[gtf$type == "gene"]
seqlevelsStyle(genes_gr) <- "UCSC"
mcols_names     <- names(mcols(genes_gr))
gene_name_field <- if ("gene_name" %in% mcols_names) "gene_name" else
                   if ("gene_id"   %in% mcols_names) "gene_id"   else
                   stop("No gene_name or gene_id field in gene records.")
gene_names <- as.character(mcols(genes_gr)[[gene_name_field]])

# =============================================================================
# LOAD CYTOBANDS
# =============================================================================

cyto <- fread(cyto_file, header=FALSE)
setnames(cyto, c("chr_raw","start","end","band","stain"))
cyto[, chr := gsub("^chr","", chr_raw)]
cyto <- cyto[chr %in% as.character(1:22)]
cyto[, chr_num := as.integer(chr)]
setkey(cyto, chr_num, start, end)

# =============================================================================
# LOAD CNV SEGMENTS (used only for cnv_region gene annotation, not carrier counting)
# =============================================================================

cnv_files <- list.files(cnv_dir, pattern="\\.txt$", full.names=TRUE)
if (!length(cnv_files)) stop("No CNV files found in: ", cnv_dir)
cat("Loading CNV segments:", length(cnv_files), "files\n")
cnv <- rbindlist(lapply(cnv_files, fread), fill=TRUE)
if (ncol(cnv) < 17) stop("CNV table has unexpected column count: ", ncol(cnv))
setnames(cnv, c("region","length","n_probes","source_file","start_probe","end_probe",
                 "conf","state","sample_qc","chr","start","end","type","col14","ckid","gene","col17"))
cnv[, chr   := as.integer(gsub("^chr","", chr))]
cnv[, start := as.numeric(start)]
cnv[, end   := as.numeric(end)]
cnv[, ckid  := as.character(ckid)]
cnv[, type  := as.character(type)]
cat("CNV rows:", nrow(cnv), "\n\n")

# =============================================================================
# CNVR BUILDING
# =============================================================================

make_cnvr_for_pheno_model <- function(sig_dt, mat, ph, ct, tag_kb, tag_r2, min_pairs) {
  x <- copy(sig_dt)
  setorder(x, pval)
  x <- x[probe_id %in% colnames(mat)]
  if (nrow(x) == 0) return(list(cnvr_dt=data.table(), membership_dt=data.table()))

  cnvrs <- list(); memb <- list(); cnvr_id_counter <- 1L

  for (chr_i in unique(x$chr_num)) {
    xx <- x[chr_num == chr_i]
    if (nrow(xx) == 0) next
    remaining <- rep(TRUE, nrow(xx))

    while (any(remaining)) {
      idx    <- which(remaining)
      lead_i <- idx[which.min(xx$pval[idx])]
      lead   <- xx[lead_i]
      lead_probe <- lead$probe_id
      lead_pos   <- lead$pos

      win_start <- lead_pos - tag_kb*1000
      win_end   <- lead_pos + tag_kb*1000
      in_win    <- idx[xx$pos[idx] >= win_start & xx$pos[idx] <= win_end]
      pw        <- xx[in_win]
      lead_vec  <- as.numeric(mat[, lead_probe])

      tagged <- data.table(
        probe_id=pw$probe_id, pos=pw$pos, pval=pw$pval, beta=pw$beta,
        se=if ("se" %in% names(pw)) pw$se else NA_real_
      )
      tagged[, r2_to_lead := NA_real_]
      for (k in seq_len(nrow(tagged))) {
        p <- tagged$probe_id[k]
        v <- as.numeric(mat[, p])
        r <- safe_cor(lead_vec, v, min_pairs=min_pairs)
        tagged$r2_to_lead[k] <- ifelse(is.na(r), NA_real_, r^2)
      }

      keep <- tagged[probe_id == lead_probe | (!is.na(r2_to_lead) & r2_to_lead >= tag_r2)]
      if (nrow(keep) == 0) keep <- tagged[probe_id == lead_probe]

      cnvr_id <- sprintf("CNVR_UNPRUNED_%s_%s_chr%d_%05d", ph, ct, chr_i, cnvr_id_counter)
      cnvr_id_counter <- cnvr_id_counter + 1L

      cnvrs[[length(cnvrs)+1]] <- data.table(
        phenotype=ph, cnv_type=ct, cnvr_id=cnvr_id, chr=chr_i,
        event_start=min(keep$pos), event_end=max(keep$pos),
        n_sig_probes_in_event=nrow(keep),
        lead_probe=lead_probe, lead_pos=lead_pos,
        pval=lead$pval, lead_beta=lead$beta,
        lead_se=if ("se" %in% names(lead)) lead$se else NA_real_,
        mean_beta=mean(keep$beta, na.rm=TRUE),
        mean_se=mean(keep$se, na.rm=TRUE)
      )

      memb[[length(memb)+1]] <- data.table(
        phenotype=ph, cnv_type=ct, cnvr_id=cnvr_id,
        probe_id=keep$probe_id, probe_pos=keep$pos, probe_pval=keep$pval,
        probe_beta=keep$beta, probe_se=keep$se, r2_to_lead=keep$r2_to_lead
      )

      remaining[in_win] <- remaining[in_win] & !(xx$probe_id[in_win] %in% keep$probe_id)
      remaining[lead_i] <- FALSE
    }
  }

  list(cnvr_dt=rbindlist(cnvrs, fill=TRUE), membership_dt=rbindlist(memb, fill=TRUE))
}

cnvr_list <- list(); memb_list <- list()

for (ct in models_present) {
  mat_file <- find_matrix_file(ct, mat_base_dir)
  cat("Loading matrix for tagging:", ct, "->", basename(mat_file), "\n")
  mat <- readRDS(mat_file)

  for (phn in unique(sig$phenotype)) {
    ssub <- sig[phenotype == phn & cnv_type == ct]
    if (nrow(ssub) == 0) next
    reg_rows <- intersect(mat_row_cache[[ct]], get_gwas_samples(phn, ct))   # regression sample
    cat("  r2 tagging sample:", phn, ct, length(reg_rows), "\n")
    res <- make_cnvr_for_pheno_model(ssub, mat[reg_rows, , drop = FALSE], phn, ct, tag_kb, tag_r2, min_pairs)
    if (nrow(res$cnvr_dt) > 0) {
      cnvr_list[[length(cnvr_list)+1]] <- res$cnvr_dt
      memb_list[[length(memb_list)+1]] <- res$membership_dt
    }
  }

  rm(mat); gc()
  cat("\n")
}

cnvr_by_model <- rbindlist(cnvr_list, fill=TRUE)
cnvr_members  <- rbindlist(memb_list, fill=TRUE)
if (nrow(cnvr_by_model) == 0) stop("No CNVRs created.")

cnvr_by_model <- merge(cnvr_by_model, tested_probe_counts, by="cnv_type", all.x=TRUE)
cnvr_by_model[, pct_sig_of_tested := round(100 * n_sig_probes_in_event / n_tested_probes, 4)]
setorder(cnvr_by_model, phenotype, cnv_type, chr, event_start)

cat("=== CNVR BY MODEL: BASIC SUMMARY ===\n")
print(cnvr_by_model[, .N, by=.(phenotype, cnv_type)][order(phenotype, cnv_type)])
cat("\n")

# =============================================================================
# AUWERX LD LOCUS BOUNDARIES
# =============================================================================

cat("Computing Auwerx LD-based locus boundaries + genes...\n")
cnvr_by_model[, `:=`(auwerx_event_start=NA_real_, auwerx_event_end=NA_real_, auwerx_genes=NA_character_)]

for (ct in models_present) {
  mat_file   <- find_matrix_file(ct, mat_base_dir)
  cat("  Loading matrix for Auwerx boundaries:", ct, "->", basename(mat_file), "\n")
  mat        <- readRDS(mat_file)
  mat_probes <- colnames(mat)
  pm_sub     <- probe_map[probe_id %in% mat_probes, .(probe_id, chr_num, pos)]
  rows_ct    <- which(cnvr_by_model$cnv_type == ct)

  for (i in rows_ct) {
    lp     <- cnvr_by_model$lead_probe[i]
    lp_pos <- cnvr_by_model$lead_pos[i]
    chr_i  <- as.integer(cnvr_by_model$chr[i])
    if (is.na(lp) || is.na(lp_pos) || is.na(chr_i)) next
    if (!lp %in% mat_probes) next

    win_s      <- lp_pos - tag_kb*1000
    win_e      <- lp_pos + tag_kb*1000
    candidates <- pm_sub[chr_num == chr_i & pos >= win_s & pos <= win_e]
    if (nrow(candidates) == 0) next

    lead_vec <- as.numeric(mat[, lp])
    r2_vals  <- vapply(candidates$probe_id, function(pid) {
      v  <- as.numeric(mat[, pid])
      ok <- !is.na(lead_vec) & !is.na(v)
      if (sum(ok) < min_pairs) return(NA_real_)
      r <- suppressWarnings(cor(lead_vec[ok], v[ok]))
      if (is.na(r)) NA_real_ else r^2
    }, numeric(1))

    candidates[, r2 := r2_vals]
    in_ld <- candidates[!is.na(r2) & r2 >= tag_r2]

    if (nrow(in_ld) > 0) {
      a_start <- min(in_ld$pos); a_end <- max(in_ld$pos)
      cnvr_by_model[i, `:=`(auwerx_event_start=a_start, auwerx_event_end=a_end)]
      g <- genes_overlapping_interval_scalar(chr_i, a_start, a_end, genes_gr, gene_names)
      if (length(g) > 0) cnvr_by_model[i, auwerx_genes := paste(sort(unique(g)), collapse=";")]
    }
  }

  rm(mat); gc()
}
cat("Auwerx locus boundaries + genes added.\n\n")

# =============================================================================
# GENE ANNOTATIONS
# =============================================================================

cat("Annotating genes (lead pos + event span)...\n")

lead_ok   <- which(!is.na(cnvr_by_model$chr) & !is.na(cnvr_by_model$lead_pos))
lead_gr   <- GRanges(seqnames=paste0("chr", cnvr_by_model$chr[lead_ok]),
                     ranges=IRanges(start=cnvr_by_model$lead_pos[lead_ok],
                                    end=cnvr_by_model$lead_pos[lead_ok]))
lead_hits <- findOverlaps(lead_gr, genes_gr, ignore.strand=TRUE)
lead_df   <- data.frame(query_idx=lead_ok[queryHits(lead_hits)],
                         gene=gene_names[subjectHits(lead_hits)], stringsAsFactors=FALSE)
cnvr_by_model[, lead_genes   := collapse_by_row(lead_df, .N)]
cnvr_by_model[, lead_in_gene := !is.na(lead_genes)]

ev_ok      <- which(!is.na(cnvr_by_model$chr) & !is.na(cnvr_by_model$event_start) & !is.na(cnvr_by_model$event_end))
ev_start   <- pmin(cnvr_by_model$event_start[ev_ok], cnvr_by_model$event_end[ev_ok])
ev_end     <- pmax(cnvr_by_model$event_start[ev_ok], cnvr_by_model$event_end[ev_ok])
event_gr   <- GRanges(seqnames=paste0("chr", cnvr_by_model$chr[ev_ok]),
                      ranges=IRanges(start=ev_start, end=ev_end))
event_hits <- findOverlaps(event_gr, genes_gr, ignore.strand=TRUE)
event_df   <- data.frame(query_idx=ev_ok[queryHits(event_hits)],
                          gene=gene_names[subjectHits(event_hits)], stringsAsFactors=FALSE)
cnvr_by_model[, event_span_genes := collapse_by_row(event_df, .N)]

cat("  lead_in_gene TRUE:", cnvr_by_model[, sum(lead_in_gene, na.rm=TRUE)], "\n\n")

# =============================================================================
# CNV MAPPING + PHENOTYPE STATS
# Carrier counts use the lead-probe dosage in the matrix, as in the regression.
# Call (segment) overlap is used only for the genes in the CNV region.
# =============================================================================

cat("Counting carriers via lead probe dosage + computing phenotype stats...\n")

add_list <- vector("list", nrow(cnvr_by_model))

for (ct in models_present) {
  mat_file <- find_matrix_file(ct, mat_base_dir)
  cat("Loading matrix for carrier counting:", ct, "->", basename(mat_file), "\n")
  mat      <- readRDS(mat_file)
  mat_rows <- extract_sample_ids(mat, ct)

  rows_ct <- which(cnvr_by_model$cnv_type == ct)

  for (i in rows_ct) {
    r          <- cnvr_by_model[i]
    pheno      <- r[["phenotype"]]
    chr_val    <- as.integer(r[["chr"]])
    lead_probe <- r[["lead_probe"]]
    aw_start   <- as.numeric(r[["auwerx_event_start"]])
    aw_end     <- as.numeric(r[["auwerx_event_end"]])

    gwas_samples <- get_gwas_samples(pheno, ct)

    total_n <- get_n_analyzed(ct, pheno, lead_probe)
    if (is.na(total_n)) {
      warning("Falling back to matrix-based N for ", pheno, "/", ct, "/", lead_probe)
      total_n <- length(gwas_samples)
    }

    # ------------------------------------------------------------------
    # CARRIER COUNTING: lead probe dosage (matches regression exactly)
    # ------------------------------------------------------------------
    common <- intersect(gwas_samples, mat_rows)

    if (length(common) == 0 || !lead_probe %in% colnames(mat)) {
      warning("No common samples or lead probe missing for ", pheno, "/", ct)
      del_carriers <- character(0)
      dup_carriers <- character(0)
    } else {
      lead_vec        <- as.numeric(mat[common, lead_probe])
      names(lead_vec) <- common
      if (ct == "deletion") {
        # deletion matrix encodes deletions as +1
        del_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
        dup_carriers <- character(0)
      } else if (ct == "duplication") {
        # duplication matrix encodes duplications as +1
        del_carriers <- character(0)
        dup_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
      } else if (ct == "mirror") {
        # mirror matrix encodes deletions as -1, duplications as +1
        del_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec < 0])
        dup_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
      }
    }

    any_carriers       <- union(del_carriers, dup_carriers)
    noncarriers_strict <- setdiff(common, any_carriers)

    ph_col <- get_pheno_col(ph_all, pheno, use_rint=use_rint)

    pull_vals <- function(ids) {
      ids <- intersect(ids, gwas_samples)
      if (!length(ids)) return(numeric(0))
      ph_all[J(ids), get(ph_col)]
    }

    st_all    <- safe_stats(pull_vals(gwas_samples))
    st_del    <- safe_stats(pull_vals(del_carriers))
    st_dup    <- safe_stats(pull_vals(dup_carriers))
    st_noncar <- safe_stats(pull_vals(noncarriers_strict))

    n_focal <- switch(ct,
      deletion    = length(del_carriers),
      duplication = length(dup_carriers),
      mirror      = length(any_carriers),
      length(any_carriers)
    )

    # ------------------------------------------------------------------
    # CNV REGION BOUNDARIES: segment overlap for gene annotation only
    # ------------------------------------------------------------------
    ov_sub <- cnv[chr == chr_val & start <= aw_end & end >= aw_start & ckid %in% gwas_samples]

    cnv_region_start <- if (nrow(ov_sub) > 0) suppressWarnings(min(ov_sub$start, na.rm=TRUE)) else NA_real_
    cnv_region_end   <- if (nrow(ov_sub) > 0) suppressWarnings(max(ov_sub$end,   na.rm=TRUE)) else NA_real_
    if (!is.finite(cnv_region_start)) cnv_region_start <- NA_real_
    if (!is.finite(cnv_region_end))   cnv_region_end   <- NA_real_

    cnv_region_genes <- if (!is.na(cnv_region_start) && !is.na(cnv_region_end)) {
      collapse_unique(genes_overlapping_interval_scalar(
        chr_val, cnv_region_start, cnv_region_end, genes_gr, gene_names))
    } else NA_character_

    add_list[[i]] <- data.table(
      pheno_col_used   = ph_col,
      total_n          = total_n,
      n_del            = length(del_carriers),
      n_dup            = length(dup_carriers),
      n_carrier        = n_focal,
      n_noncarrier     = total_n - n_focal,
      carrier_freq_pct = if (total_n > 0) round(100 * n_focal / total_n, 3) else NA_real_,
      mean_all         = st_all$mean,    median_all        = st_all$median,    sd_all        = st_all$sd,
      mean_del         = st_del$mean,    median_del        = st_del$median,    sd_del        = st_del$sd,
      mean_dup         = st_dup$mean,    median_dup        = st_dup$median,    sd_dup        = st_dup$sd,
      mean_noncarrier  = st_noncar$mean, median_noncarrier = st_noncar$median, sd_noncarrier = st_noncar$sd,
      cnv_region_start = cnv_region_start,
      cnv_region_end   = cnv_region_end,
      cnv_region_genes = cnv_region_genes
    )

    cat(sprintf("  [%d/%d] %s | %s | n_del=%d n_dup=%d n_carrier=%d\n",
                i, nrow(cnvr_by_model), pheno, ct,
                length(del_carriers), length(dup_carriers), n_focal))
  }

  rm(mat); gc()
  cat("\n")
}

extra      <- rbindlist(add_list, fill=TRUE)
cnvr_annot <- cbind(cnvr_by_model, extra)

cnvr_annot[cnv_type == "deletion",    c("mean_dup","median_dup","sd_dup") := .(NA_real_, NA_real_, NA_real_)]
cnvr_annot[cnv_type == "duplication", c("mean_del","median_del","sd_del") := .(NA_real_, NA_real_, NA_real_)]

# =============================================================================
# CYTOBANDS
# =============================================================================

cnvr_annot[, chr_num := as.integer(chr)]
cnvr_annot[, row_id  := .I]
cnvr_annot[, lead_pos_int := as.integer(round(lead_pos))]

lead_dt <- cnvr_annot[!is.na(lead_pos_int) & !is.na(chr_num),
                      .(row_id, chr_num, start=lead_pos_int, end=lead_pos_int+1L)]
setkey(lead_dt, chr_num, start, end)
lead_hit <- foverlaps(lead_dt, cyto, by.x=c("chr_num","start","end"),
                      by.y=c("chr_num","start","end"), type="within", nomatch=NA)
cnvr_annot[, lead_cytoband := NA_character_]
cnvr_annot[lead_hit$row_id, lead_cytoband := paste0(lead_hit$chr_num, lead_hit$band)]

make_cytoband_range <- function(dt, start_col, end_col, out_col) {
  range_dt <- dt[!is.na(get(start_col)) & !is.na(get(end_col)) & !is.na(chr_num),
                 .(row_id, chr_num, start=as.integer(get(start_col)), end=as.integer(get(end_col)))]
  if (nrow(range_dt) == 0L) { dt[, (out_col) := NA_character_]; return(dt) }
  setkey(range_dt, chr_num, start, end)
  range_hits <- foverlaps(range_dt, cyto, by.x=c("chr_num","start","end"),
                          by.y=c("chr_num","start","end"), type="any", nomatch=NA)
  setorder(range_hits, row_id, start, end)
  range_lab <- range_hits[!is.na(band), .(band_first=band[1], band_last=band[.N]), by=row_id]
  dt[, (out_col) := NA_character_]
  dt[range_lab, on="row_id",
     (out_col) := fifelse(i.band_first == i.band_last,
                          paste0(chr_num, i.band_first),
                          paste0(chr_num, i.band_first, "–", i.band_last))]
  dt
}

cnvr_annot <- make_cytoband_range(cnvr_annot, "event_start",        "event_end",        "sig_probe_cytoband_range")
cnvr_annot <- make_cytoband_range(cnvr_annot, "auwerx_event_start", "auwerx_event_end", "auwerx_cytoband_range")
cnvr_annot <- make_cytoband_range(cnvr_annot, "cnv_region_start",   "cnv_region_end",   "cnv_region_cytoband_range")
cnvr_annot[, locus_cytoband_final := lead_cytoband]

# =============================================================================
# EXPERIMENT-WIDE THRESHOLD
# =============================================================================

cat("\n=== SENSITIVITY: EXPERIMENT-WIDE MULTIPLE TESTING ===\n")
cnvr_annot[, gw_thr_model      := 0.05 / meff_map[cnv_type]]
cnvr_annot[, gw_thr_experiment := 0.05 / (meff_map[cnv_type] * n_traits_tested)]
cnvr_annot[, pass_gw_model      := !is.na(pval) & pval <= gw_thr_model]
cnvr_annot[, pass_gw_experiment := !is.na(pval) & pval <= gw_thr_experiment]

sens_sum <- cnvr_annot[, .(
  n_signals=.N,
  n_pass_model=sum(pass_gw_model, na.rm=TRUE),
  n_pass_experiment=sum(pass_gw_experiment, na.rm=TRUE),
  thr_model=unique(gw_thr_model),
  thr_experiment=unique(gw_thr_experiment)
), by=cnv_type][order(cnv_type)]
print(sens_sum)

cat("\nSignals passing model-specific GW threshold:",
    sum(cnvr_annot$pass_gw_model, na.rm=TRUE), "of", nrow(cnvr_annot), "\n")
cat("Signals passing experiment-wide threshold:",
    sum(cnvr_annot$pass_gw_experiment, na.rm=TRUE), "of", nrow(cnvr_annot), "\n\n")

# =============================================================================
# PRINT FINAL TABLE
# =============================================================================

cat("\n=== FINAL: CNVRs BY MODEL ===\n")
print_cols <- intersect(c(
  "phenotype","cnv_type","chr","event_start","event_end",
  "auwerx_event_start","auwerx_event_end","auwerx_genes",
  "pval","lead_beta","lead_se","lead_pos","mean_beta","mean_se",
  "n_sig_probes_in_event","n_tested_probes","pct_sig_of_tested",
  "total_n","cnv_region_start","cnv_region_end",
  "n_del","n_dup","n_carrier","n_noncarrier","carrier_freq_pct",
  "lead_in_gene","lead_genes","event_span_genes","cnv_region_genes",
  "lead_cytoband","sig_probe_cytoband_range","auwerx_cytoband_range",
  "cnv_region_cytoband_range","locus_cytoband_final","cnvr_id","lead_probe"
), names(cnvr_annot))

setorder(cnvr_annot, phenotype, cnv_type, chr, event_start)
print(cnvr_annot[, ..print_cols])
cat("\nRows:", nrow(cnvr_annot), "\n")

# =============================================================================
# MERGED CNVRs ACROSS MODELS
# =============================================================================

cat("\n=== MERGED CNVRs ACROSS MODELS (mirror-only) ===\n")
tmp <- copy(cnvr_annot); tmp[, idx := .I]
edges <- list()
for (phn in unique(tmp$phenotype)) {
  dtp <- tmp[phenotype == phn]
  if (nrow(dtp) < 2) next
  for (i in 1:(nrow(dtp)-1)) {
    for (j in (i+1):nrow(dtp)) {
      a <- dtp[i]; b <- dtp[j]
      if (a$chr != b$chr) next
      if (!intervals_overlap(a$event_start, a$event_end, b$event_start, b$event_end)) next
      if (!directionally_concordant(a$cnv_type, a$lead_beta, b$cnv_type, b$lead_beta, mirror_matches_dup)) next
      edges[[length(edges)+1]] <- data.table(phenotype=phn, from=a$idx, to=b$idx)
    }
  }
}
edges_dt <- if (length(edges)) rbindlist(edges) else data.table()
cat("Merge edges found:", nrow(edges_dt), "\n")

tmp[, merged_comp := idx]
if (nrow(edges_dt) > 0) {
  for (phn in unique(edges_dt$phenotype)) {
    dtp    <- tmp[phenotype == phn]
    parent <- setNames(dtp$idx, dtp$idx)
    findp  <- function(x) { while (parent[[as.character(x)]] != x) x <- parent[[as.character(x)]]; x }
    unionp <- function(a,b) { ra <- findp(a); rb <- findp(b); if (ra != rb) parent[[as.character(rb)]] <<- ra }
    ee <- edges_dt[phenotype == phn]
    for (k in seq_len(nrow(ee))) unionp(ee$from[k], ee$to[k])
    comp <- sapply(dtp$idx, findp)
    tmp[phenotype == phn, merged_comp := comp]
  }
}

tmp[, merged_cnvr_id := paste0("MERGED_UNPRUNED_", phenotype, "_", merged_comp)]

merged <- tmp[, .(
  phenotype=phenotype[1], chr=chr[1],
  merged_start=min(event_start), merged_end=max(event_end),
  models_involved=paste(sort(unique(cnv_type)), collapse=","),
  n_source_cnvrs=.N, min_lead_p=min(pval, na.rm=TRUE),
  lead_probe_best=lead_probe[which.min(pval)],
  lead_beta_best=lead_beta[which.min(pval)],
  cytoband_final=locus_cytoband_final[which.min(pval)]
), by=merged_cnvr_id]

best_model_dt <- tmp[, .(best_model=cnv_type[which.min(pval)]), by=merged_cnvr_id]
merged <- merge(merged, best_model_dt, by="merged_cnvr_id", all.x=TRUE)
merged[, gw_thr_experiment  := 0.05 / (meff_map[best_model] * n_traits_tested)]
merged[, pass_gw_experiment := !is.na(min_lead_p) & min_lead_p <= gw_thr_experiment]

cat("\n=== SENSITIVITY ON MERGED LOCI ===\n")
cat("Merged loci passing experiment-wide threshold:",
    sum(merged$pass_gw_experiment, na.rm=TRUE), "of", nrow(merged), "\n")
setorder(merged, phenotype, chr, merged_start)
print(merged)
cat("\nDone.\n")

# =============================================================================
# SAVE OUTPUTS
# =============================================================================

save_dir <- file.path(sig_dir, "cnvr_outputs")
dir.create(save_dir, showWarnings=FALSE, recursive=TRUE)

premerge_file <- file.path(save_dir, paste0("cnvr_by_model_supplemental_unpruned_", relatedness_group, ".csv"))
fwrite(cnvr_annot, premerge_file)
cat("Saved pre-merge CNVR table:\n  ", premerge_file, "\n", sep="")

postmerge_file <- file.path(save_dir, paste0("cnvr_merged_loci_unpruned_", relatedness_group, ".csv"))
fwrite(merged, postmerge_file)
cat("Saved merged CNVR table:\n  ", postmerge_file, "\n", sep="")

if (exists("cnvr_members") && nrow(cnvr_members) > 0) {
  membership_file <- file.path(save_dir, paste0("cnvr_probe_membership_unpruned_", relatedness_group, ".csv"))
  fwrite(cnvr_members, membership_file)
  cat("Saved CNVR probe membership:\n  ", membership_file, "\n", sep="")
}

cat("All outputs saved in:\n  ", save_dir, "\n", sep="")

