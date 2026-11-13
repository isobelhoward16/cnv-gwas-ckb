#!/usr/bin/env Rscript
# Stage 9: boxplots of the trait by copy-number group at the lead probe of every stage 5 signal (Figures 4B, 4C,
# 5B, 5D and the rest in results/09_figures/<version>/boxplots/). Participants as in the association analysis
# (passing CNV QC, with phenotype data); carrier groups with up to 100 participants are also drawn as points.
# Usage: Rscript scripts/09_figures/04_boxplots.R [related | unrelated]  (default related = all participants)
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggbeeswarm)
  library(stringr)
  library(gridExtra)
  library(Matrix)
})

source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version

# ========= ARGUMENTS =========
args <- commandArgs(trailingOnly = TRUE)
relatedness <- ifelse(length(args) > 0, args[1], "related")

cat("=== CNV BOXPLOT GENERATOR (UNPRUNED, WITH SAMPLE POINTS) ===\n")
cat("Relatedness:", toupper(relatedness), "\n\n")

# ========= CONFIG =========
use_rint <- FALSE

# ========= PATHS =========
cnvr_table <- file.path(pa$results_root, "05_signals", ver, if (relatedness == "related") "all" else "unrelated", "cnvr_outputs", paste0("cnvr_by_model_supplemental_unpruned_", relatedness, ".csv"))

mat_base_dir   <- file.path(pa$results_root, "02_matrices", ver)
pheno_rds      <- pa$phenotype_rds
key1_path      <- pa$array1_key
key2_path      <- pa$array2_key
qc_master_path <- pa$qc_gwas_file

outdir <- file.path(pa$results_root, "09_figures", ver, "boxplots", paste0("unpruned_", relatedness, "_withpoints"))
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ========= HELPERS =========

get_pheno_col <- function(ph_all, pheno_base, use_rint = FALSE) {
  pheno_base <- gsub("_rint$", "", pheno_base)
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

find_matrix_file <- function(cnv_type, mat_base_dir) {
  f <- file.path(mat_base_dir, paste0("combined_", cnv_type, "_matrix_shared_probes_fixed.rds"))
  if (file.exists(f)) return(f)
  f2 <- file.path(mat_base_dir, paste0("combined_", cnv_type, "_matrix_shared_probes.rds"))
  if (file.exists(f2)) return(f2)
  patt  <- paste0("^combined_", cnv_type, "_matrix_.*shared_probes.*\\.rds$")
  files <- list.files(mat_base_dir, pattern = patt, full.names = TRUE)
  files <- files[!grepl("pruned", files)]
  if (length(files) == 0) stop("No unpruned matrix found for cnv_type=", cnv_type)
  files[which.max(file.info(files)$mtime)]
}

extract_sample_ids <- function(obj) {
  rn <- tryCatch(rownames(obj), error = function(e) NULL)
  if (!is.null(rn) && length(rn) > 0) return(rn)
  dn <- tryCatch(dimnames(obj), error = function(e) NULL)
  if (!is.null(dn) && length(dn) >= 1 && !is.null(dn[[1]]) && length(dn[[1]]) > 0) return(dn[[1]])
  if (isS4(obj)) {
    dn2 <- tryCatch(obj@Dimnames[[1]], error = function(e) NULL)
    if (!is.null(dn2) && length(dn2) > 0) return(dn2)
  }
  stop("Cannot extract sample IDs from matrix")
}

# ========= CLEAN PHENOTYPE LABELS =========
get_phenotype_label <- function(phenotype_code) {
  if (is.null(phenotype_code) || phenotype_code == "") return("Unknown phenotype")
  base <- gsub("_rint$", "", phenotype_code)
  label_map <- c(
    standing_height_mm = "Standing height",
    sitting_height_mm  = "Sitting height",
    weight_kg          = "Weight",
    bmi_calc           = "BMI",
    bmi_at_25          = "BMI at age 25",
    fat_percent        = "Body fat",
    waist_mm           = "Waist circumference",
    hip_mm             = "Hip circumference",
    waist_hip_ratio    = "Waist-hip ratio",
    dbp_mean           = "DBP",
    sbp_mean           = "SBP",
    heart_rate_mean    = "Resting heart rate",
    random_glucose     = "RPG",
    fasting_glucose    = "Fasting glucose"
  )
  if (!is.null(label_map[[base]])) return(label_map[[base]])
  stringr::str_to_title(gsub("_", " ", base))
}

get_phenotype_units <- function(phenotype_code) {
  units_map <- list(
    "bmi_calc" = "kg/m\u00b2", "bmi" = "kg/m\u00b2", "bmi_at_25" = "kg/m\u00b2",
    "dbp_mean" = "mmHg", "sbp_mean" = "mmHg", "fat_percent" = "%",
    "heart_rate_mean" = "bpm", "hip_mm" = "cm", "waist_mm" = "cm",
    "standing_height_mm" = "cm", "sitting_height_mm" = "cm", "weight_kg" = "kg",
    "random_glucose" = "mmol/L", "fasting_glucose" = "mmol/L"
  )
  search_code <- tolower(gsub("_rint$", "", phenotype_code))
  if (!is.null(units_map[[search_code]])) return(units_map[[search_code]])
  for (k in names(units_map)) {
    if (grepl(k, search_code, fixed = TRUE) || grepl(search_code, k, fixed = TRUE))
      return(units_map[[k]])
  }
  ""
}

# ========= CLEAN LOCUS LABEL =========
get_locus_label <- function(row) {
  cyto     <- as.character(row$locus_cytoband_final)
  if (identical(cyto, "1q21.2")) cyto <- "1q21.1"   # sitting height: same 1q21.1 distal CNV as standing height
  aw_start <- suppressWarnings(as.numeric(row$auwerx_event_start))
  model    <- as.character(row$cnv_type)

  if (!is.na(cyto) && grepl("^16p11\\.2", cyto)) {
    bp_label  <- if (!is.na(aw_start) && aw_start < 29300000) "BP2\u2013BP3" else "BP4\u2013BP5"
    locus_str <- paste0("16p11.2 ", bp_label)
  } else {
    locus_str <- cyto
  }

  model_str <- switch(model,
    deletion    = "Deletion",
    duplication = "Duplication",
    mirror      = "Mirror",
    model
  )

  list(locus = locus_str, model = model_str)
}

format_int_commas <- function(x) {
  format(as.integer(round(x)), big.mark = ",", scientific = FALSE, trim = TRUE)
}

get_box_whiskers <- function(x) {
  x <- x[!is.na(x)]
  if (!length(x)) return(c(NA_real_, NA_real_))
  q1  <- as.numeric(quantile(x, 0.25, na.rm = TRUE, type = 7))
  q3  <- as.numeric(quantile(x, 0.75, na.rm = TRUE, type = 7))
  iqr <- q3 - q1
  c(min(x[x >= q1 - 1.5 * iqr], na.rm = TRUE),
    max(x[x <= q3 + 1.5 * iqr], na.rm = TRUE))
}

fmt_mm <- function(mn, med) {
  if (is.na(mn) || is.na(med)) return(NA_character_)
  paste0(round(mn, 2), " (", round(med, 2), ")")
}

# ========= LOAD DATA =========

load_phenotype_data <- function() {
  cat("Loading phenotype data...\n")
  ph <- readRDS(pheno_rds); setDT(ph)
  if (relatedness == "unrelated" && "unrelated" %in% names(ph))
    ph <- ph[unrelated == 1]
  k1 <- fread(key1_path); k2 <- fread(key2_path)
  setnames(k1, c("ccvid", "xIID_cnv")); setnames(k2, c("ccvid", "xIID_cnv"))
  ph1 <- merge(ph[gwas_array_type == "Axiom_CKB1"], k1, by = "ccvid", all = FALSE)
  ph2 <- merge(ph[gwas_array_type == "Axiom_CKB2"], k2, by = "ccvid", all = FALSE)
  ph_all <- rbindlist(list(ph1, ph2), use.names = TRUE, fill = TRUE)

  # ---- QC filter, to match the association analysis sample ----
  cat("Loading QC master file...\n")
  qc <- fread(qc_master_path)
  qc_pass <- qc[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
  n_before <- nrow(ph_all)
  ph_all <- ph_all[xIID_cnv %in% qc_pass]
  cat(sprintf("  QC filter: %d -> %d samples (removed %d)\n",
              n_before, nrow(ph_all), n_before - nrow(ph_all)))

  setkey(ph_all, xIID_cnv)
  cat("\u2713 Phenotype data loaded:", nrow(ph_all), "samples\n\n")
  ph_all
}

load_matrices <- function(models_needed) {
  cat("Loading full matrices...\n")
  mat_cache <- list()
  for (ct in models_needed) {
    f <- find_matrix_file(ct, mat_base_dir)
    cat("  ", ct, "->", basename(f), "\n")
    mat_cache[[ct]] <- readRDS(f)
  }
  cat("\n")
  mat_cache
}

get_gwas_samples <- function(ph_all, pheno_base, cnv_type, mat_cache) {
  col      <- get_pheno_col(ph_all, pheno_base, use_rint = use_rint)
  ph_sub   <- ph_all[!is.na(get(col)), xIID_cnv]
  mat_rows <- extract_sample_ids(mat_cache[[cnv_type]])
  intersect(ph_sub, mat_rows)
}

# ========= FIND CARRIERS VIA LEAD PROBE DOSAGE =========

find_cnv_groups <- function(row, ph_all, mat_cache) {
  cnv_model  <- as.character(row$cnv_type)
  phenotype  <- as.character(row$phenotype)
  lead_probe <- as.character(row$lead_probe)

  mat <- mat_cache[[cnv_model]]
  if (is.null(mat)) return(NULL)
  if (!lead_probe %in% colnames(mat)) {
    cat("    -> Lead probe not in matrix\n")
    return(NULL)
  }

  gwas_samples <- tryCatch(
    get_gwas_samples(ph_all, phenotype, cnv_model, mat_cache),
    error = function(e) NULL
  )
  if (is.null(gwas_samples) || length(gwas_samples) == 0) return(NULL)

  mat_rows <- extract_sample_ids(mat)
  common   <- intersect(gwas_samples, mat_rows)
  if (length(common) == 0) return(NULL)

  lead_vec        <- as.numeric(mat[common, lead_probe])
  names(lead_vec) <- common

  if (cnv_model == "deletion") {
    del_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
    neutral      <- names(lead_vec[!is.na(lead_vec) & lead_vec == 0])
    group_dt <- data.table(
      ckid     = c(del_carriers, neutral),
      cn_group = c(rep("Del", length(del_carriers)), rep("Copy neutral", length(neutral)))
    )
    group_levels <- c("Copy neutral", "Del")
    n_focal      <- length(del_carriers)

  } else if (cnv_model == "duplication") {
    dup_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
    neutral      <- names(lead_vec[!is.na(lead_vec) & lead_vec == 0])
    group_dt <- data.table(
      ckid     = c(neutral, dup_carriers),
      cn_group = c(rep("Copy neutral", length(neutral)), rep("Dup", length(dup_carriers)))
    )
    group_levels <- c("Copy neutral", "Dup")
    n_focal      <- length(dup_carriers)

  } else if (cnv_model == "mirror") {
    del_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec < 0])
    dup_carriers <- names(lead_vec[!is.na(lead_vec) & lead_vec > 0])
    neutral      <- names(lead_vec[!is.na(lead_vec) & lead_vec == 0])
    group_dt <- data.table(
      ckid     = c(del_carriers, neutral, dup_carriers),
      cn_group = c(rep("Del", length(del_carriers)),
                   rep("Copy neutral", length(neutral)),
                   rep("Dup", length(dup_carriers)))
    )
    group_levels <- c("Del", "Copy neutral", "Dup")
    n_focal      <- length(del_carriers) + length(dup_carriers)

  } else return(NULL)

  group_dt <- unique(group_dt, by = "ckid")
  group_dt[, cn_group := factor(cn_group, levels = group_levels)]

  list(
    group_dt  = group_dt,
    pheno_col = get_pheno_col(ph_all, phenotype, use_rint = use_rint),
    total_n   = length(common),
    n_focal   = n_focal,
    del_n     = if (cnv_model == "deletion") length(del_carriers) else
                if (cnv_model == "mirror")   length(del_carriers) else 0L,
    dup_n     = if (cnv_model == "duplication") length(dup_carriers) else
                if (cnv_model == "mirror")      length(dup_carriers) else 0L,
    neutral_n = length(neutral)
  )
}

# ========= CREATE BOXPLOT =========

mm_vars <- c("standing_height_mm", "sitting_height_mm", "waist_mm", "hip_mm")

create_boxplot <- function(ph_all, row, mat_cache) {
  groups <- find_cnv_groups(row, ph_all, mat_cache)
  if (is.null(groups)) return(NULL)

  pheno_col   <- groups$pheno_col
  pheno_label <- get_phenotype_label(as.character(row$phenotype))
  units       <- get_phenotype_units(as.character(row$phenotype))
  locus_info  <- get_locus_label(row)

  # carrier n from the plot groups
  n_carriers_actual <- groups$del_n + groups$dup_n

  plot_dt <- merge(
    groups$group_dt,
    ph_all[, .(xIID_cnv, pheno_val = get(pheno_col))],
    by.x = "ckid", by.y = "xIID_cnv", all.x = FALSE
  )
  setnames(plot_dt, "pheno_val", "phenotype_value")

  # Convert mm -> cm for relevant phenotypes
  if (gsub("_rint$", "", as.character(row$phenotype)) %in% mm_vars) {
    plot_dt[, phenotype_value := phenotype_value / 10]
  }

  plot_dt <- plot_dt[!is.na(phenotype_value)]
  if (nrow(plot_dt) == 0) return(NULL)

  plot_dt[, cn_group := droplevels(cn_group)]

  full_name_map <- c(
    "Del"          = "Deletion",
    "Dup"          = "Duplication",
    "Copy neutral" = "Copy neutral"
  )
  count_dt <- plot_dt[, .(n = .N), by = cn_group]
  count_dt[, x_label := paste0(
    full_name_map[as.character(cn_group)],
    "\n(n=", format(n, big.mark = ",", trim = TRUE), ")"
  )]
  x_lab_map <- setNames(count_dt$x_label, as.character(count_dt$cn_group))

  # individual points only for carrier groups with <= point_max_n participants (per group);
  # larger groups are summarised by the box alone
  point_max_n <- 100
  small_carrier_groups <- count_dt[cn_group != "Copy neutral" & n <= point_max_n, cn_group]
  point_dt <- plot_dt[cn_group %in% small_carrier_groups]

  whisk_mat <- do.call(rbind, lapply(split(plot_dt$phenotype_value, plot_dt$cn_group), get_box_whiskers))
  ylow  <- min(whisk_mat[, 1], na.rm = TRUE)
  yhigh <- max(whisk_mat[, 2], na.rm = TRUE)
  # y range: whiskers of every group, extended to the points actually drawn
  if (nrow(point_dt) > 0) {
    ylow  <- min(ylow,  min(point_dt$phenotype_value, na.rm = TRUE))
    yhigh <- max(yhigh, max(point_dt$phenotype_value, na.rm = TRUE))
  }
  if (!is.finite(ylow) || !is.finite(yhigh)) {
    ylow  <- min(plot_dt$phenotype_value, na.rm = TRUE)
    yhigh <- max(plot_dt$phenotype_value, na.rm = TRUE)
  }
  yrange <- yhigh - ylow
  if (!is.finite(yrange) || yrange <= 0) yrange <- max(abs(yhigh), 1) * 0.1
  pad <- yrange * 0.08

  # title "{locus} - {phenotype}", subtitle "{model} - lead p - n carriers"
  title_txt    <- paste0(locus_info$locus, "  \u2022  ", pheno_label)
  subtitle_txt <- paste0(
    locus_info$model, "  \u00b7  ",
    "lead p\u2009=\u2009", format(as.numeric(row$pval), scientific = TRUE, digits = 2),
    "  \u00b7  n\u2009=\u2009", format(n_carriers_actual, big.mark = ","), " carriers"
  )

  fill_vals <- c("Del" = "#D95F0E", "Copy neutral" = "#BDBDBD", "Dup" = "#2C7FB8")
  fill_vals <- fill_vals[levels(plot_dt$cn_group)]
  point_vals <- c("Del" = "#8C3B06", "Copy neutral" = "#4D4D4D", "Dup" = "#134E75")
  point_vals <- point_vals[levels(plot_dt$cn_group)]

  # point size and transparency from the largest carrier group
  carrier_n <- count_dt[cn_group != "Copy neutral", n]
  n_max <- if (length(carrier_n) > 0) max(carrier_n) else max(count_dt$n)
  pt_size  <- if (n_max > 150) 1.0 else if (n_max > 40) 1.25 else 1.7
  pt_alpha <- if (n_max > 150) 0.28 else if (n_max > 40) 0.35 else 0.55

  ggplot(plot_dt, aes(x = cn_group, y = phenotype_value)) +
    geom_boxplot(aes(fill = cn_group), outlier.shape = NA, width = 0.6,
                 alpha = 0.5, linewidth = 0.7) +
    ggbeeswarm::geom_quasirandom(data = point_dt, aes(color = cn_group),
                 width = 0.22, size = pt_size, alpha = pt_alpha,
                 stroke = 0, method = "pseudorandom") +
    stat_summary(fun = mean, geom = "point", shape = 21, size = 3.0,
                 fill = "black", color = "white", stroke = 0.9) +   # mean
    scale_fill_manual(values = fill_vals, drop = TRUE) +
    scale_color_manual(values = point_vals, drop = TRUE) +
    scale_x_discrete(labels = x_lab_map, drop = TRUE) +
    scale_y_continuous(expand = c(0, 0)) +
    coord_cartesian(ylim = c(ylow - pad, yhigh + pad)) +
    labs(
      title    = title_txt,
      subtitle = subtitle_txt,
      x        = "Copy number state",
      y        = ifelse(units == "", pheno_label, paste0(pheno_label, " (", units, ")"))
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold", size = 17,
                                     lineheight = 1.05, margin = margin(b = 4)),
      plot.subtitle   = element_text(hjust = 0.5, size = 14, color = "grey30"),
      axis.title.x    = element_text(face = "bold", size = 17, margin = margin(t = 10)),
      axis.title.y    = element_text(face = "bold", size = 17, margin = margin(r = 10)),
      axis.text.x     = element_text(size = 15, margin = margin(t = 6)),
      axis.text.y     = element_text(size = 15),
      axis.line       = element_line(linewidth = 0.7),
      legend.position = "none",
      plot.margin     = margin(14, 20, 10, 14)
    )
}

# ========= BUILD SUMMARY ROWS =========

build_summary_row_single_type <- function(locus, phenotype, model,
                                          cnvr_start, cnvr_end, lead_pos,
                                          grp_stats) {
  carrier_grp   <- if (model == "deletion") "Del" else "Dup"
  carrier_stats <- grp_stats[[carrier_grp]]
  neutral_stats <- grp_stats[["Copy neutral"]]

  data.frame(
    Locus                        = locus,
    Phenotype                    = phenotype,
    "CNVR start"                 = format_int_commas(cnvr_start),
    "CNVR end"                   = format_int_commas(cnvr_end),
    "Lead probe pos"             = format_int_commas(lead_pos),
    "N carriers"                 = if (!is.null(carrier_stats)) carrier_stats$N  else 0L,
    "N non-carriers"             = if (!is.null(neutral_stats)) neutral_stats$N  else NA_integer_,
    "Carriers mean (median)"     = if (!is.null(carrier_stats) && carrier_stats$N > 0)
                                     fmt_mm(carrier_stats$Mean, carrier_stats$Median)
                                   else NA_character_,
    "Non-carriers mean (median)" = if (!is.null(neutral_stats))
                                     fmt_mm(neutral_stats$Mean, neutral_stats$Median)
                                   else NA_character_,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

build_summary_row_mirror <- function(locus, phenotype,
                                     cnvr_start, cnvr_end, lead_pos,
                                     grp_stats) {
  del_stats     <- grp_stats[["Del"]]
  dup_stats     <- grp_stats[["Dup"]]
  neutral_stats <- grp_stats[["Copy neutral"]]

  data.frame(
    Locus                       = locus,
    Phenotype                   = phenotype,
    "CNVR start"                = format_int_commas(cnvr_start),
    "CNVR end"                  = format_int_commas(cnvr_end),
    "Lead probe pos"            = format_int_commas(lead_pos),
    "N del"                     = if (!is.null(del_stats))     del_stats$N     else 0L,
    "N dup"                     = if (!is.null(dup_stats))     dup_stats$N     else 0L,
    "N non-carrier"             = if (!is.null(neutral_stats)) neutral_stats$N else NA_integer_,
    "Del mean (median)"         = if (!is.null(del_stats) && del_stats$N > 0)
                                    fmt_mm(del_stats$Mean, del_stats$Median)
                                  else "0 (0)",
    "Dup mean (median)"         = if (!is.null(dup_stats) && dup_stats$N > 0)
                                    fmt_mm(dup_stats$Mean, dup_stats$Median)
                                  else "0 (0)",
    "Non-carrier mean (median)" = if (!is.null(neutral_stats))
                                    fmt_mm(neutral_stats$Mean, neutral_stats$Median)
                                  else NA_character_,
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

# ========= MAIN =========

main <- function() {
  cat("CNVR table:", cnvr_table, "\n")
  cat("Output dir:", outdir, "\n\n")

  if (!file.exists(cnvr_table)) stop("CNVR table not found: ", cnvr_table)
  cnvr <- fread(cnvr_table)
  if (nrow(cnvr) == 0) stop("CNVR table is empty")

  req_cols <- c("phenotype", "cnv_type", "chr", "pval", "event_start", "event_end",
                "lead_probe", "lead_pos", "cnvr_id", "locus_cytoband_final",
                "auwerx_event_start", "auwerx_event_end", "n_carrier")
  missing_cols <- setdiff(req_cols, names(cnvr))
  if (length(missing_cols)) stop("Missing columns: ", paste(missing_cols, collapse = ", "))

  ph_all    <- load_phenotype_data()
  cnvr      <- cnvr[cnv_type %in% c("deletion", "duplication", "mirror") &
                    !is.na(chr) & chr %in% 1:22]
  mat_cache <- load_matrices(sort(unique(cnvr$cnv_type)))

  cat("Rows to plot:", nrow(cnvr), "\n\n")

  all_plots <- list()
  del_rows  <- list()
  dup_rows  <- list()
  mir_rows  <- list()

  for (i in seq_len(nrow(cnvr))) {
    r <- cnvr[i]
    cat(sprintf("[%d/%d] %s | %s | chr%s | %s\n",
                i, nrow(cnvr), r$phenotype, r$cnv_type, r$chr, r$cnvr_id))

    locus_info <- get_locus_label(r)
    phenotype  <- gsub("_rint$", "", as.character(r$phenotype))
    model      <- as.character(r$cnv_type)
    cnvr_start <- as.integer(r$auwerx_event_start)
    cnvr_end   <- as.integer(r$auwerx_event_end)
    lead_pos   <- as.integer(r$lead_pos)

    p <- tryCatch(
      create_boxplot(ph_all, r, mat_cache),
      error = function(e) { cat("    -> ERROR:", conditionMessage(e), "\n"); NULL }
    )
    if (!is.null(p)) {
      fname <- sprintf("boxplot_%s_%s_chr%s_%s.png",
                       r$phenotype, r$cnv_type, r$chr,
                       gsub("[^A-Za-z0-9_-]", "_", r$cnvr_id))
      ggsave(file.path(outdir, fname), p, width = 8, height = 6, dpi = 300, bg = "white")
      cat("    -> Saved", fname, "\n")
      all_plots[[fname]] <- p
    }

    grps <- tryCatch(find_cnv_groups(r, ph_all, mat_cache), error = function(e) NULL)
    if (!is.null(grps)) {
      plot_dt2 <- merge(
        grps$group_dt,
        ph_all[, .(xIID_cnv, pheno_val = get(grps$pheno_col))],
        by.x = "ckid", by.y = "xIID_cnv", all.x = FALSE
      )
      setnames(plot_dt2, "pheno_val", "phenotype_value")
      plot_dt2 <- plot_dt2[!is.na(phenotype_value)]

      if (nrow(plot_dt2) > 0) {
        grp_stats <- list()
        for (grp in levels(plot_dt2$cn_group)) {
          vals <- plot_dt2[cn_group == grp, phenotype_value]
          if (length(vals) > 0) {
            grp_stats[[grp]] <- list(
              N      = length(vals),
              Mean   = mean(vals,   na.rm = TRUE),
              Median = median(vals, na.rm = TRUE),
              SD     = sd(vals,     na.rm = TRUE)
            )
          }
        }

        if (model == "deletion") {
          del_rows[[length(del_rows) + 1]] <- build_summary_row_single_type(
            locus_info$locus, phenotype, model, cnvr_start, cnvr_end, lead_pos, grp_stats)
        } else if (model == "duplication") {
          dup_rows[[length(dup_rows) + 1]] <- build_summary_row_single_type(
            locus_info$locus, phenotype, model, cnvr_start, cnvr_end, lead_pos, grp_stats)
        } else if (model == "mirror") {
          mir_rows[[length(mir_rows) + 1]] <- build_summary_row_mirror(
            locus_info$locus, phenotype, cnvr_start, cnvr_end, lead_pos, grp_stats)
        }
      }
    }
  }

  if (length(all_plots) == 0) {
    cat("\n No plots created\n")
    return(invisible(NULL))
  }

  # ========= WRITE BOXPLOT PDF grouped by phenotype (up to 4 per page) =========
  pdf_file <- file.path(outdir, paste0("all_boxplots_unpruned_", relatedness, "_withpoints.pdf"))

  fname_list <- names(all_plots)
  plot_pheno_map <- data.table(
    fname     = fname_list,
    phenotype = vapply(fname_list, function(nm) {
      matched <- cnvr[sprintf("boxplot_%s_%s_chr%s_%s.png",
                               phenotype, cnv_type, chr,
                               gsub("[^A-Za-z0-9_-]", "_", cnvr_id)) == nm]
      if (nrow(matched) > 0) gsub("_rint$", "", matched$phenotype[1]) else "unknown"
    }, character(1))
  )
  setorder(plot_pheno_map, phenotype)

  pdf(pdf_file, width = 16, height = 12)
  for (ph in unique(plot_pheno_map$phenotype)) {
    ph_label  <- get_phenotype_label(ph)
    ph_fnames <- plot_pheno_map[phenotype == ph, fname]
    pages <- split(ph_fnames, ceiling(seq_along(ph_fnames) / 4))
    for (pg_plots in pages) {
      plot_list <- lapply(pg_plots, function(nm) all_plots[[nm]])
      while (length(plot_list) < 4) plot_list <- c(plot_list, list(ggplot() + theme_void()))
      gridExtra::grid.arrange(
        grobs = plot_list, ncol = 2, nrow = 2,
        top   = grid::textGrob(
          ph_label,
          gp = grid::gpar(fontsize = 18, fontface = "bold")
        )
      )
    }
  }
  dev.off()
  cat("\n\u2713 Boxplot PDF written to:", pdf_file, "\n")

  # ========= WRITE SUMMARY TABLE PDF =========
  tbl_theme <- gridExtra::ttheme_default(
    base_size = 11,
    padding   = unit(c(6, 4), "mm"),
    core      = list(fg_params = list(cex = 0.85)),
    colhead   = list(fg_params = list(cex = 0.85, fontface = "bold"))
  )

  make_tbl_page <- function(df, title_str) {
    tbl <- gridExtra::tableGrob(df, rows = NULL, theme = tbl_theme)
    gridExtra::grid.arrange(
      tbl,
      top = grid::textGrob(title_str, gp = grid::gpar(fontsize = 14, fontface = "bold"))
    )
  }

  tbl_file <- file.path(outdir, paste0("summary_table_unpruned_", relatedness, "_withpoints.pdf"))
  pdf(tbl_file, width = 24, height = 16)
  if (length(del_rows) > 0) make_tbl_page(do.call(rbind, del_rows), "Deletion model signals")
  if (length(dup_rows) > 0) make_tbl_page(do.call(rbind, dup_rows), "Duplication model signals")
  if (length(mir_rows) > 0) make_tbl_page(do.call(rbind, mir_rows), "Mirror model signals")
  dev.off()
  cat("\u2713 Summary table PDF written to:", tbl_file, "\n")
  cat("\u2713 Done\n")
}

main()
