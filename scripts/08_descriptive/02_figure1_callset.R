#!/usr/bin/env Rscript
# Stage 8, step 2: Figure 1 (CNV burden, length and copy-number state across the callset) and the
# summary tables behind it, overall and by array phase.

suppressPackageStartupMessages({
  library(data.table)
  library(gtable)
  library(ggplot2)
  library(stringr)
  library(yaml)
  library(scales)
  library(cowplot)
  library(grid)
})

# -----------------------------
# Helpers
# -----------------------------
dir_create <- function(x) {
  if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
}

timestamp_string <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")

normalize_cnv_type <- function(x) {
  x <- tolower(trimws(x))
  ifelse(
    x %in% c("del", "deletion", "loss"), "del",
    ifelse(x %in% c("dup", "duplication", "gain"), "dup", NA_character_)
  )
}

parse_chr_int <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^chr", "")
  suppressWarnings(as.integer(x))
}

save_plot_png <- function(p, out_prefix, width = 8, height = 6, dpi = 600) {
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = width,
    height = height,
    dpi = dpi,
    units = "in",
    bg = "white"
  )
}

format_bp_kb_mb <- function(x) {
  if (is.na(x)) return(NA_character_)
  if (x >= 1e6) {
    paste0(round(x / 1e6, 2), " Mb")
  } else {
    paste0(round(x / 1e3, 1), " kb")
  }
}

format_mb <- function(x) {
  if (is.na(x)) return(NA_character_)
  paste0(round(x, 2), " Mb")
}

format_kb_from_mb <- function(x) {
  if (is.na(x)) return(NA_character_)
  paste0(round(x * 1e3), " kb")
}

summ_stats <- function(dt, value_col, group_cols = NULL, label = NULL) {
  v <- value_col
  if (is.null(group_cols)) {
    dt[, .(
      metric = label,
      N      = .N,
      mean   = as.double(mean(get(v),                    na.rm = TRUE)),
      median = as.double(median(get(v),                  na.rm = TRUE)),
      p25    = as.double(quantile(get(v), 0.25,          na.rm = TRUE)),
      p75    = as.double(quantile(get(v), 0.75,          na.rm = TRUE)),
      p95    = as.double(quantile(get(v), 0.95,          na.rm = TRUE))
    )]
  } else {
    dt[, .(
      metric = label,
      N      = .N,
      mean   = as.double(mean(get(v),                    na.rm = TRUE)),
      median = as.double(median(get(v),                  na.rm = TRUE)),
      p25    = as.double(quantile(get(v), 0.25,          na.rm = TRUE)),
      p75    = as.double(quantile(get(v), 0.75,          na.rm = TRUE)),
      p95    = as.double(quantile(get(v), 0.95,          na.rm = TRUE))
    ), by = group_cols]
  }
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)

  w <- suppressWarnings(wilcox.test(x, y, exact = FALSE)$statistic)

  n1 <- as.numeric(length(x))
  n2 <- as.numeric(length(y))
  U <- as.numeric(w) - n1 * (n1 + 1) / 2

  (2 * U) / (n1 * n2) - 1
}

ks_p <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  if (length(x) == 0 || length(y) == 0) return(NA_real_)
  suppressWarnings(ks.test(x, y)$p.value)
}

# -----------------------------
# Config
# -----------------------------
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)

out_root             <- paths$out_root
gene_cleaned_merged_dir <- paths$gene_cleaned_merged_dir
qc_master_file          <- paths$qc_master_file
manifest_path <- file.path(out_root, "data/derived/manifests/postQC_denominator_samples.tsv")

tables_dir   <- file.path(out_root, "data/derived/tables")
fig_main_dir <- file.path(out_root, "results/figures/main")
fig_supp_dir <- file.path(out_root, "results/figures/supp")
logs_dir     <- file.path(out_root, "logs")

dir_create(tables_dir)
dir_create(fig_main_dir)
dir_create(fig_supp_dir)
dir_create(logs_dir)

run_ts   <- timestamp_string()
log_file <- file.path(logs_dir, paste0("landscape_burden_basics_", run_ts, ".log"))
zz       <- file(log_file, open = "wt")
sink(zz, split = TRUE)

on.exit({
  sink()
  close(zz)
}, add = TRUE)

cat("=== Paper1: CNV burden basics (aligned composite figure) ===\n")
cat("Started: ", format(Sys.time()), "\n")
cat("gene_cleaned_merged_dir: ", gene_cleaned_merged_dir, "\n")
cat("manifest: ", manifest_path, "\n")
cat("qc_master_file: ", qc_master_file, "\n")
cat("log: ", log_file, "\n\n")

# -----------------------------
# Denominator manifest
# -----------------------------
denom_dt <- fread(manifest_path)
stopifnot(all(c("ckid", "phase") %in% names(denom_dt)))
denom_dt[, phase := factor(phase, levels = c("Phase12", "Phase3"))]

cat("Denominator N from manifest (pre-QC filter): ", nrow(denom_dt), "\n")
print(denom_dt[, .N, by = phase])
cat("\n")

# ---- Apply QC filter from master QC file ----
cat("Loading QC master file:", qc_master_file, "\n")
qc_master <- fread(qc_master_file)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_ids <- qc_master[qc_pass == TRUE]$GID
cat("QC-passing samples in master file:", length(qc_pass_ids), "\n")

before_n <- nrow(denom_dt)
denom_dt <- denom_dt[ckid %in% qc_pass_ids]
cat("Removed", before_n - nrow(denom_dt), "QC-failing samples from denominator\n")
cat("Denominator N after QC filter:", nrow(denom_dt), "\n")
print(denom_dt[, .N, by = phase])
cat("\n")

# -----------------------------
# CNV ingest (final callset; headerless)
# -----------------------------
CHR_COL   <- 10L
START_COL <- 11L
END_COL   <- 12L
TYPE_COL  <- 13L
CKID_COL  <- 15L
CN_COL    <- 9L

cnv_files <- list.files(gene_cleaned_merged_dir, pattern = "\\.txt$", full.names = TRUE)
stopifnot(length(cnv_files) > 0)
cat("CNV files found: ", length(cnv_files), "\n")

peek <- fread(cnv_files[1], sep = "\t", header = FALSE, nrows = 3, fill = TRUE, quote = "")
cat("Peek ncol(first file): ", ncol(peek), "\n")
stopifnot(ncol(peek) >= CKID_COL)

read_cnv_file <- function(f) {
  dt <- fread(
    f,
    sep = "\t",
    header = FALSE,
    fill = TRUE,
    quote = "",
    showProgress = FALSE
  )

  phase <- if (str_detect(basename(f), "Phase12")) {
    "Phase12"
  } else if (str_detect(basename(f), "Phase3")) {
    "Phase3"
  } else {
    NA_character_
  }

  data.table(
    file = basename(f),
    phase = phase,
    ckid = as.character(dt[[CKID_COL]]),
    chr = parse_chr_int(dt[[CHR_COL]]),
    start = suppressWarnings(as.integer(dt[[START_COL]])),
    end = suppressWarnings(as.integer(dt[[END_COL]])),
    cnv_type = normalize_cnv_type(as.character(dt[[TYPE_COL]])),
    copy_number = suppressWarnings(as.integer(dt[[CN_COL]]))
  )
}

cnv_dt <- rbindlist(lapply(cnv_files, read_cnv_file), use.names = TRUE, fill = TRUE)

cat("CNV records loaded: ", nrow(cnv_dt), "\n")
cat("Unique CKIDs in CNV records: ", uniqueN(cnv_dt$ckid), "\n")

before_cnv_n <- nrow(cnv_dt)
cnv_dt <- cnv_dt[ckid %in% qc_pass_ids]
cat("CNV records after restricting to QC-passing samples:", nrow(cnv_dt),
    "(removed", before_cnv_n - nrow(cnv_dt), "records from QC-failing samples)\n")

cnv_dt <- cnv_dt[
  !is.na(phase) &
    !is.na(ckid) & ckid != "" &
    !is.na(chr) & chr >= 1 & chr <= 22 &
    !is.na(start) & !is.na(end) & start <= end &
    !is.na(cnv_type) &
    !is.na(copy_number)
]

cnv_dt[, cnv_len_bp := as.numeric(end - start)]
cnv_dt <- cnv_dt[cnv_len_bp >= 0]
cnv_dt[, len_bp_plot := pmax(cnv_len_bp, 1)]
cnv_dt[, len_mb_plot := pmax(cnv_len_bp / 1e6, 1e-6)]

cat("After filtering: ", nrow(cnv_dt), " CNVs\n")
print(cnv_dt[, .N, by = .(phase, cnv_type)][order(phase, cnv_type)])
cat("\n")

cat("Copy-number state counts overall:\n")
print(cnv_dt[, .N, by = copy_number][order(copy_number)])
cat("\n")

cat("Copy-number by CNV type:\n")
print(cnv_dt[, .N, by = .(cnv_type, copy_number)][order(cnv_type, copy_number)])
cat("\n")

# -----------------------------
# Per-sample burden (include zero-CNV)
# -----------------------------
burden <- cnv_dt[, .(
  n_cnv = .N,
  bp_burden = sum(cnv_len_bp, na.rm = TRUE)
), by = .(ckid, phase)]

burden_all <- merge(
  denom_dt[, .(ckid, phase)],
  burden,
  by = c("ckid", "phase"),
  all.x = TRUE
)

burden_all[is.na(n_cnv), n_cnv := 0L]
burden_all[is.na(bp_burden), bp_burden := 0]
burden_all[, mb_burden := bp_burden / 1e6]
burden_all[, is_carrier := n_cnv > 0]

burden_carriers      <- burden_all[is_carrier == TRUE]
burden_carriers_plot <- burden_carriers[mb_burden > 0]

cat("Carrier counts (from burden table):\n")
print(burden_all[, .(N = .N, carriers = sum(is_carrier)), by = phase])
cat("Overall carriers: ", sum(burden_all$is_carrier), " / ", nrow(burden_all), "\n\n")

# -----------------------------
# Save tables
# -----------------------------
out_burden <- file.path(tables_dir, paste0("per_sample_cnv_burden_", run_ts, ".tsv"))
fwrite(burden_all, out_burden, sep = "\t")
cat("Saved per-sample burden table (incl zeros): ", out_burden, "\n")

stats_all <- rbindlist(list(
  summ_stats(burden_all, "n_cnv", NULL, "n_cnv_per_person_all_samples"),
  summ_stats(burden_all, "n_cnv", "phase", "n_cnv_per_person_all_samples"),
  summ_stats(burden_all, "mb_burden", NULL, "mb_burden_per_person_all_samples"),
  summ_stats(burden_all, "mb_burden", "phase", "mb_burden_per_person_all_samples"),
  summ_stats(burden_carriers, "mb_burden", NULL, "mb_burden_per_person_carriers_only"),
  summ_stats(burden_carriers, "mb_burden", "phase", "mb_burden_per_person_carriers_only"),
  summ_stats(cnv_dt, "cnv_len_bp", NULL, "cnv_length_bp_per_call"),
  summ_stats(cnv_dt, "cnv_len_bp", c("phase", "cnv_type"), "cnv_length_bp_per_call")
), use.names = TRUE, fill = TRUE)

out_stats <- file.path(tables_dir, paste0("burden_basics_summary_stats_", run_ts, ".tsv"))
fwrite(stats_all, out_stats, sep = "\t")
cat("Saved summary stats: ", out_stats, "\n")

phase_comp <- rbindlist(list(
  data.table(
    metric = "n_cnv_per_person_all_samples",
    ks_p = ks_p(
      burden_all[phase == "Phase12"]$n_cnv,
      burden_all[phase == "Phase3"]$n_cnv
    ),
    cliffs_delta = cliffs_delta(
      burden_all[phase == "Phase12"]$n_cnv,
      burden_all[phase == "Phase3"]$n_cnv
    )
  ),
  data.table(
    metric = "mb_burden_per_person_carriers_only",
    ks_p = ks_p(
      burden_carriers[phase == "Phase12"]$mb_burden,
      burden_carriers[phase == "Phase3"]$mb_burden
    ),
    cliffs_delta = cliffs_delta(
      burden_carriers[phase == "Phase12"]$mb_burden,
      burden_carriers[phase == "Phase3"]$mb_burden
    )
  ),
  data.table(
    metric = "cnv_length_log10bp_per_call",
    ks_p = ks_p(
      log10(cnv_dt[phase == "Phase12"]$len_bp_plot),
      log10(cnv_dt[phase == "Phase3"]$len_bp_plot)
    ),
    cliffs_delta = cliffs_delta(
      log10(cnv_dt[phase == "Phase12"]$len_bp_plot),
      log10(cnv_dt[phase == "Phase3"]$len_bp_plot)
    )
  )
), use.names = TRUE, fill = TRUE)

out_phasecomp <- file.path(tables_dir, paste0("burden_phase_comparison_stats_", run_ts, ".tsv"))
fwrite(phase_comp, out_phasecomp, sep = "\t")
cat("Saved phase comparison stats: ", out_phasecomp, "\n\n")

copy_number_stats <- cnv_dt[, .N, by = .(copy_number)][order(copy_number)]
out_copy_number <- file.path(tables_dir, paste0("copy_number_state_counts_", run_ts, ".tsv"))
fwrite(copy_number_stats, out_copy_number, sep = "\t")
cat("Saved copy-number state counts: ", out_copy_number, "\n\n")

# -----------------------------
# Plot styling
# -----------------------------
phase_cols <- c("Phase12" = "#2E86AB", "Phase3" = "#E67E22")
type_cols  <- c("del" = "#2E86AB", "dup" = "#E67E22")

ANNOT_SIZE       <- 5.2
COUNT_LABEL_SIZE <- 3.8
PHASE_ANNOT_SIZE <- 5.0

AXIS_TITLE_SIZE <- 18
AXIS_TEXT_SIZE  <- 15

theme_pub <- theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.6),
    plot.title = element_blank(),
    legend.position = "none",
    axis.title = element_text(size = AXIS_TITLE_SIZE),
    axis.text = element_text(size = AXIS_TEXT_SIZE),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.margin = margin(10, 10, 8, 8)
  )

# -----------------------------
# Plot 1: CNVs per person (all samples; truncated at 30)
# -----------------------------
x_cap        <- 30
n_gt_cap     <- burden_all[n_cnv > x_cap, .N]
median_n_cnv <- median(burden_all$n_cnv, na.rm = TRUE)
mean_n_cnv   <- mean(burden_all$n_cnv, na.rm = TRUE)

hist_p1 <- ggplot_build(
  ggplot(burden_all[n_cnv <= x_cap], aes(x = n_cnv)) +
    geom_histogram(bins = x_cap, boundary = 0)
)$data[[1]]
ymax_p1 <- max(hist_p1$count, na.rm = TRUE)

p1 <- ggplot(burden_all[n_cnv <= x_cap], aes(x = n_cnv)) +
  geom_histogram(
    bins = x_cap,
    boundary = 0,
    fill = "#2E86AB",
    color = "#2E86AB",
    alpha = 0.75
  ) +
  theme_pub +
  labs(x = "CNVs per person", y = "Individuals") +
  annotate("text", x = Inf, y = Inf,
           label = paste0(">30: ", format(n_gt_cap, big.mark = ",")),
           hjust = 1.02, vjust = 1.4, size = ANNOT_SIZE, color = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median: ", round(median_n_cnv, 2)),
           hjust = 1.02, vjust = 2.8, size = ANNOT_SIZE, color = "#2E86AB") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mean: ", round(mean_n_cnv, 2)),
           hjust = 1.02, vjust = 4.2, size = ANNOT_SIZE, color = "#7B3294")

save_plot_png(p1, file.path(fig_main_dir, "Fig_burden_basics_main_A_CNVs_per_individual"),
              width = 7.2, height = 4.5)

# -----------------------------
# Plot 2: Total CNV burden among carriers (Mb; log10 x)
# -----------------------------
median_mb_burden <- median(burden_carriers_plot$mb_burden, na.rm = TRUE)
mean_mb_burden   <- mean(burden_carriers_plot$mb_burden, na.rm = TRUE)

hist_p2 <- ggplot_build(
  ggplot(burden_carriers_plot, aes(x = mb_burden)) +
    geom_histogram(bins = 60)
)$data[[1]]
ymax_p2 <- max(hist_p2$count, na.rm = TRUE)

p2 <- ggplot(burden_carriers_plot, aes(x = mb_burden)) +
  geom_histogram(bins = 60, fill = "#2E86AB", color = "#2E86AB", alpha = 0.75) +
  geom_vline(xintercept = median_mb_burden, color = "black",  linewidth = 1, linetype = "dotted") +
  geom_vline(xintercept = mean_mb_burden,   color = "red3",   linewidth = 1, linetype = "dotted") +
  scale_x_log10(labels = function(x) label_number()(x * 1e3)) +
  theme_pub +
  labs(x = "Total CNV burden (kb; log10)", y = "Individuals") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median: ", format_kb_from_mb(median_mb_burden)),
           hjust = 1.02, vjust = 1.4, size = ANNOT_SIZE, color = "black") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Mean: ", format_kb_from_mb(mean_mb_burden)),
           hjust = 1.02, vjust = 2.8, size = ANNOT_SIZE, color = "red3")

save_plot_png(p2, file.path(fig_main_dir, "Fig_burden_basics_main_B_total_burden_carriers"),
              width = 10.5, height = 4.5)

# -----------------------------
# Plot 3: CNV length distribution per call (del vs dup; Mb on x-axis)
# -----------------------------
del_median_len_mb     <- cnv_dt[cnv_type == "del", median(len_mb_plot, na.rm = TRUE)]
dup_median_len_mb     <- cnv_dt[cnv_type == "dup", median(len_mb_plot, na.rm = TRUE)]
overall_median_len_mb <- median(cnv_dt$len_mb_plot, na.rm = TRUE)

hist_p3 <- ggplot_build(
  ggplot(cnv_dt, aes(x = len_mb_plot, fill = cnv_type)) +
    geom_histogram(aes(y = after_stat(density)), bins = 80, position = "identity")
)$data[[1]]
ymax_p3 <- max(hist_p3$density, na.rm = TRUE)

p3 <- ggplot(cnv_dt, aes(x = len_mb_plot, color = cnv_type, fill = cnv_type)) +
  geom_histogram(aes(y = after_stat(density)), bins = 80, position = "identity",
                 alpha = 0.35, show.legend = FALSE) +
  geom_vline(xintercept = del_median_len_mb,     color = type_cols["del"], linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = dup_median_len_mb,     color = type_cols["dup"], linewidth = 1, linetype = "dashed") +
  geom_vline(xintercept = overall_median_len_mb, color = "#7B3294",        linewidth = 1, linetype = "solid") +
  scale_x_log10(labels = function(x) label_number()(x * 1e3)) +
  scale_color_manual(values = type_cols, guide = "none") +
  scale_fill_manual(values = type_cols, guide = "none") +
  theme_pub +
  labs(x = "CNV length (kb; log10)", y = "Density") +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median deletion length: ", format_kb_from_mb(del_median_len_mb)),
           hjust = 1.02, vjust = 1.4, size = ANNOT_SIZE, color = type_cols["del"]) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Median duplication length: ", format_kb_from_mb(dup_median_len_mb)),
           hjust = 1.02, vjust = 2.8, size = ANNOT_SIZE, color = type_cols["dup"]) +
  annotate("text", x = Inf, y = Inf,
           label = paste0("Overall median length: ", format_kb_from_mb(overall_median_len_mb)),
           hjust = 1.02, vjust = 4.2, size = ANNOT_SIZE, color = "#7B3294")

save_plot_png(p3, file.path(fig_main_dir, "Fig_burden_basics_main_C_cnv_length_distribution"),
              width = 10.5, height = 4.5)

# -----------------------------
# Plot 4: Copy-number state counts (main figure D)
# -----------------------------
copy_num_counts <- cnv_dt[, .N, by = .(copy_number)][order(copy_number)]

copy_num_counts <- merge(
  data.table(copy_number = sort(unique(cnv_dt$copy_number))),
  copy_num_counts,
  by = "copy_number",
  all.x = TRUE
)

copy_num_counts[is.na(N), N := 0L]
copy_num_counts[, copy_number := as.integer(copy_number)]

copy_num_counts[, cnv_class := fifelse(
  copy_number %in% c(0L, 1L), "del",
  fifelse(copy_number %in% c(3L, 4L), "dup", "other")
)]

copy_num_counts[, copy_number_f := factor(copy_number, levels = sort(unique(copy_number)))]

p4 <- ggplot(copy_num_counts, aes(x = copy_number_f, y = N, fill = cnv_class, color = cnv_class)) +
  geom_col(alpha = 0.85, show.legend = FALSE, width = 0.8) +
  geom_text(aes(label = paste0("n=", comma(N))),
            vjust = -0.4, color = "black", size = COUNT_LABEL_SIZE, show.legend = FALSE) +
  scale_fill_manual(values = c(del = "#2E86AB", dup = "#E67E22", other = "grey70"), guide = "none") +
  scale_color_manual(values = c(del = "#2E86AB", dup = "#E67E22", other = "grey70"), guide = "none") +
  expand_limits(y = max(copy_num_counts$N) * 1.12) +
  theme_pub +
  labs(x = "Copy-number state", y = "Number of CNVs")

save_plot_png(p4, file.path(fig_main_dir, "Fig_burden_basics_main_D_copy_number_state_counts"),
              width = 7.2, height = 4.5)

# -----------------------------
# Composite main figure (strict axis alignment)
# Layout:
#   top    = D | C  (side by side)
#   middle = A      (full width)
#   bottom = B      (full width)
# -----------------------------
gA <- ggplotGrob(p1)
gB <- ggplotGrob(p2)
gC <- ggplotGrob(p3)
gD <- ggplotGrob(p4)

# Align left margins of B and C (the two full-width rows)
left_cols <- 1:7
shared_left_widths <- unit.pmax(gB$widths[left_cols], gC$widths[left_cols])
gB$widths[left_cols] <- shared_left_widths
gC$widths[left_cols] <- shared_left_widths

panel_col_B <- gB$layout[gB$layout$name == "panel", "l"][1]
panel_col_C <- gC$layout[gC$layout$name == "panel", "l"][1]

shared_panel_width <- unit.pmax(gB$widths[panel_col_B], gC$widths[panel_col_C])
gB$widths[panel_col_B] <- shared_panel_width
gC$widths[panel_col_C] <- shared_panel_width

right_cols_B <- (panel_col_B + 1):length(gB$widths)
right_cols_C <- (panel_col_C + 1):length(gC$widths)

n_right <- min(length(right_cols_B), length(right_cols_C))

shared_right_widths <- unit.pmax(
  gB$widths[right_cols_B[seq_len(n_right)]],
  gC$widths[right_cols_C[seq_len(n_right)]]
)

gB$widths[right_cols_B[seq_len(n_right)]] <- shared_right_widths
gC$widths[right_cols_C[seq_len(n_right)]] <- shared_right_widths

# New layout:
#   Row 1: D (copy-number states, left) | A (CNVs per person, right) -> labelled A | B
#   Row 2: C (CNV length, full width)                                 -> labelled C
#   Row 3: B (total burden carriers, full width)                      -> labelled D

top_left  <- ggdraw(gD) + draw_plot_label("A", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 16)
top_right <- ggdraw(gA) + draw_plot_label("B", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 16)
mid_row   <- ggdraw(gC) + draw_plot_label("C", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 16)
bot_row   <- ggdraw(gB) + draw_plot_label("D", x = 0.01, y = 0.98, hjust = 0, vjust = 1, size = 16)

top_row <- plot_grid(top_left, top_right, nrow = 1, rel_widths = c(1, 1))

final_main <- plot_grid(
  top_row, NULL, mid_row, NULL, bot_row,
  ncol = 1,
  rel_heights = c(0.95, 0.06, 1, 0.06, 1)
)

save_plot_png(final_main, file.path(fig_main_dir, "Fig_burden_basics_main_composite"),
              width = 10.5, height = 12.2)

# -----------------------------
# Plot 5: Phase overlay CNVs per person
# -----------------------------
n_gt_cap_p12 <- burden_all[phase == "Phase12" & n_cnv > x_cap, .N]
n_gt_cap_p3  <- burden_all[phase == "Phase3"  & n_cnv > x_cap, .N]

p5 <- ggplot(burden_all[n_cnv <= x_cap], aes(x = n_cnv, color = phase, fill = phase)) +
  geom_histogram(bins = x_cap, boundary = 0, position = "identity",
                 alpha = 0.25, show.legend = FALSE) +
  scale_color_manual(values = phase_cols, guide = "none") +
  scale_fill_manual(values  = phase_cols, guide = "none") +
  theme_pub +
  labs(x = "CNVs per person", y = "Individuals") +
  annotate("text", x = x_cap * 0.72, y = ymax_p1 * 0.92,
           label = paste0(">30 Phase12: ", format(n_gt_cap_p12, big.mark = ","),
                          "\n>30 Phase3: ",  format(n_gt_cap_p3, big.mark = ",")),
           hjust = 0, size = PHASE_ANNOT_SIZE)

save_plot_png(p5, file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_A_CNVs_per_individual"),
              width = 7.5, height = 4.5)

# -----------------------------
# Plot 6: Phase overlay total burden among carriers
# -----------------------------
p6 <- ggplot(burden_carriers_plot, aes(x = mb_burden, color = phase, fill = phase)) +
  geom_histogram(bins = 60, position = "identity", alpha = 0.25, show.legend = FALSE) +
  scale_x_log10(labels = label_number()) +
  scale_color_manual(values = phase_cols, guide = "none") +
  scale_fill_manual(values  = phase_cols, guide = "none") +
  theme_pub +
  labs(x = "Total CNV burden (Mb; log10)", y = "Individuals")

save_plot_png(p6, file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_B_total_burden_carriers"),
              width = 10.5, height = 4.5)

# -----------------------------
# Plot 7: Phase overlay CNV length distribution
# -----------------------------
p7 <- ggplot(cnv_dt, aes(x = len_mb_plot, color = phase)) +
  geom_density(linewidth = 0.9, adjust = 1.2, show.legend = FALSE) +
  scale_x_log10(labels = label_number()) +
  scale_color_manual(values = phase_cols, guide = "none") +
  theme_pub +
  labs(x = "CNV length (Mb; log10)", y = "Density")

save_plot_png(p7, file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_C_cnv_length_distribution"),
              width = 10.5, height = 4.5)

cat("Saved PNG plots:\n")
cat("  ", file.path(fig_main_dir, "Fig_burden_basics_main_A_CNVs_per_individual.png"), "\n")
cat("  ", file.path(fig_main_dir, "Fig_burden_basics_main_B_total_burden_carriers.png"), "\n")
cat("  ", file.path(fig_main_dir, "Fig_burden_basics_main_C_cnv_length_distribution.png"), "\n")
cat("  ", file.path(fig_main_dir, "Fig_burden_basics_main_D_copy_number_state_counts.png"), "\n")
cat("  ", file.path(fig_main_dir, "Fig_burden_basics_main_composite.png"), "\n")
cat("  ", file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_A_CNVs_per_individual.png"), "\n")
cat("  ", file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_B_total_burden_carriers.png"), "\n")
cat("  ", file.path(fig_supp_dir, "SFig_burden_basics_phase_overlay_C_cnv_length_distribution.png"), "\n")

cat("\nCompleted: ", format(Sys.time()), "\n")

cat("\n=== CNV length medians shown in main Figure panel C ===\n")
cat("Median deletion length (Mb): ", format(del_median_len_mb, digits = 10), "\n", sep = "")
cat("Median duplication length (Mb): ", format(dup_median_len_mb, digits = 10), "\n", sep = "")
cat("Overall median length (Mb): ", format(overall_median_len_mb, digits = 10), "\n", sep = "")
