#!/usr/bin/env Rscript
# Stage 8, step 6: Figure S6 and the CKB vs UK Biobank comparison columns of Table S4 (carrier frequencies of
# the recurrent genomic-disorder CNVs; UK Biobank counts from Crawford et al.).

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggrepel)
  library(scales)
  library(yaml)
})

# -----------------------------
# Config / paths
# -----------------------------
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)

out_root <- paths$out_root

tables_main <- file.path(out_root, "results/tables/main")
fig_main    <- file.path(out_root, "results/figures/main")
fig_supp    <- file.path(out_root, "results/figures/supp")
logs_dir    <- file.path(out_root, "logs")

dir.create(tables_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_main,    recursive = TRUE, showWarnings = FALSE)
dir.create(fig_supp,    recursive = TRUE, showWarnings = FALSE)
dir.create(logs_dir,    recursive = TRUE, showWarnings = FALSE)

run_ts  <- format(Sys.time(), "%Y-%m-%d_%H%M%S")
log_file <- file.path(logs_dir, paste0("03_pathogenic_compare_", run_ts, ".log"))
zz <- file(log_file, open = "wt")
sink(zz, split = TRUE)

cat("=== Paper1: CKB vs UKB pathogenic CNV frequencies (Crawford loci) ===\n")
cat("Started: ", format(Sys.time()), "\n")
cat("paths.yaml: ", paths_yaml, "\n")
cat("log: ", log_file, "\n\n")

# CKB and UKB denominators
CKB_N <- as.integer(readLines(file.path(out_root, "data/derived/tables/ckb_denominator_full_qc.txt")))   # written by step 5
cat("CKB_N (from Table S4 step): ", CKB_N, "\n")
UKB_N <- 421274   # UK Biobank participants in Crawford et al.

# -----------------------------
# Inputs: CKB table (should be 54 loci)
# -----------------------------
ckb_table_path <- file.path(out_root, "results/tables/main/Table_pathogenic_loci_ckb.tsv")
stopifnot(file.exists(ckb_table_path))
ckb <- fread(ckb_table_path)

stopifnot(all(c("cnv_locus", "carriers_Overall", "freq_percent_Overall") %in% names(ckb)))

ckb_cmp <- ckb[, .(
  CNV_locus = cnv_locus,
  ckb_carriers = as.integer(carriers_Overall),
  ckb_freq_percent = as.numeric(freq_percent_Overall)
)]

# -----------------------------
# UKB (Crawford Supp Table 4 values)
# -----------------------------
ukb_complete <- data.table(
  CNV_locus = c(
    "TAR_del","TAR_dup","1q21.1del","1q21.1dup","NRXN1del","2q11.2del","2q11.2dup",
    "2q13del","2q13dup","2q13del(NPHP1)","2q13dup(NPHP1)","2q21.1del","2q21.1dup",
    "3q29del","3q29dup","WBS_dup","7q11.23dup_distal","8p23.1dup","10q11.21q11.23del",
    "10q11.21q11.23dup","10q23dup","13q12.12del","13q12.12dup","13q12del(CRYL1)",
    "13q12dup(CRYL1)","15q11.2del","15q11.2dup","PWS_dup","15q11q13del_BP3-BP4",
    "15q11q13dup_BP3-BP4","15q11q13dup_BP3-BP5","15q13.3del","15q13.3dup",
    "15q13.3del(CHRNA7)","15q13.3dup(CHRNA7)","15q24dup","16p11.2del","16p11.2dup",
    "16p11.2distal_del","16p11.2distal_dup","16p12.1del","16p12.1dup","16p13.11del",
    "16p13.11dup","17p12del(HNPP)","17p12dup(CMT1A)","Potocki_Lupski","17q11.2del(NF1)",
    "17q12del","17q12dup","22q11.2del","22q11.2dup","22q11.2distal_del","22q11.2distal_dup"
  ),
  carriers_ukb = c(
    75,436,113,177,163,31,29,53,71,2448,1976,41,59,9,5,14,24,6,57,41,7,
    85,236,379,10,1664,2041,19,16,53,9,42,240,10,3031,9,110,138,58,137,
    246,202,131,828,237,124,5,9,9,101,10,280,5,13
  ),
  freq_percent_ukb = c(
    0.018,0.103,0.027,0.042,0.039,0.007,0.007,0.013,0.017,0.581,0.469,0.010,0.014,
    0.002,0.001,0.003,0.006,0.001,0.014,0.010,0.002,0.020,0.056,0.090,0.002,0.395,
    0.484,0.005,0.004,0.013,0.002,0.010,0.057,0.002,0.719,0.002,0.026,0.033,0.014,
    0.033,0.058,0.048,0.031,0.197,0.056,0.029,0.001,0.002,0.002,0.024,0.002,0.066,
    0.001,0.003
  )
)

# Align naming for Potocki
ukb_complete[CNV_locus == "Potocki_Lupski", CNV_locus := "Potocki-Lupski"]

# -----------------------------
# Merge + write table (ALL 54)
# -----------------------------
df <- merge(ukb_complete, ckb_cmp, by = "CNV_locus", all.x = TRUE)

# If any are NA in CKB (shouldn't be, but safe), treat as 0 for frequency table display
df[is.na(ckb_carriers), ckb_carriers := 0L]
df[is.na(ckb_freq_percent), ckb_freq_percent := 0]

# -----------------------------
# Frequencies from carrier counts and the cohort denominators:
# fold-change = (CKB carriers / CKB_N) / (UKB carriers / UKB_N).
# Computed from counts because the UKB percentages in Crawford et al.
# (Supplementary Table 4) are rounded to 3 decimal places. The rounded
# percentages are used only for the frequency table and plot axes.
# -----------------------------
df[, ckb_freq_exact := ckb_carriers / CKB_N * 100]
df[, ukb_freq_exact := carriers_ukb / UKB_N * 100]

cat("--- Merge summary ---\n")
cat("UKB loci: ", nrow(ukb_complete), "\n")
cat("CKB loci present: ", nrow(ckb_cmp), "\n")
cat("CKB zeros (carriers=0): ", sum(df$ckb_carriers == 0), "\n\n")

out_cmp_table <- file.path(tables_main, "Table_pathogenic_ckb_vs_ukb.tsv")
fwrite(df, out_cmp_table, sep = "\t")
cat("Saved comparison table: ", out_cmp_table, "\n\n")

# -----------------------------
# Plot helper
# -----------------------------
save_plot_both <- function(p, out_prefix, width = 7.5, height = 6.5, dpi = 600) {
  ggsave(paste0(out_prefix, ".png"), p, width = width, height = height, dpi = dpi, units = "in")
  ggsave(paste0(out_prefix, ".pdf"), p, width = width, height = height, units = "in")
  cat("Saved:", paste0(out_prefix, ".(png|pdf)"), "\n")
}

# -----------------------------
# (A) Log-log scatter: NON-ZERO ONLY (no floor hacks)
# Axes plot raw measured frequency (rounded display columns are
# fine here; log2_fc below drives label selection only, and is
# computed from the exact columns).
# -----------------------------
df_scatter <- df[ckb_freq_percent > 0 & freq_percent_ukb > 0]
df_scatter[, log2_fc := log2(ckb_freq_exact / ukb_freq_exact)]
df_scatter[, abs_log2_fc := abs(log2_fc)]

# Label rule: only top N outliers among loci with decent frequency
label_min_freq <- 0.01   # only consider loci >=0.01% in either cohort for labels
label_top_n <- 10

cand <- df_scatter[pmax(ckb_freq_percent, freq_percent_ukb) >= label_min_freq]
lab_dt <- cand[order(-abs_log2_fc)][1:min(.N, label_top_n)]

p_scatter <- ggplot(df_scatter, aes(x = freq_percent_ukb, y = ckb_freq_percent)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.9) +
  geom_point(size = 2.7, alpha = 0.9) +
  geom_text_repel(
    data = lab_dt,
    aes(label = CNV_locus),
    size = 3.3,
    max.overlaps = Inf,
    box.padding = 0.35,
    min.segment.length = 0
  ) +
  scale_x_log10(labels = label_number(accuracy = 0.001)) +
  scale_y_log10(labels = label_number(accuracy = 0.001)) +
  coord_equal() +
  theme_classic(base_size = 12) +
  theme(
    axis.line = element_line(linewidth = 0.6)
  ) +
  labs(
    x = "UKB carrier frequency (%)",
    y = "CKB carrier frequency (%)"
  )

save_plot_both(p_scatter, file.path(fig_main, "Fig_pathogenic_ckb_vs_ukb_loglog_scatter_nonzero"),
               width = 7.8, height = 6.8)

# -----------------------------
# (B) log2 fold-change plot: NON-ZERO ONLY (same set as scatter).
# Labels eight loci with large differences; fold-changes are from
# exact carrier counts and cohort denominators.
# -----------------------------
df_fc <- copy(df_scatter)

main_text_loci <- c(
  "2q13dup(NPHP1)", "13q12del(CRYL1)", "13q12dup(CRYL1)",
  "15q13.3del(CHRNA7)", "15q11q13dup_BP3-BP5", "17q11.2del(NF1)",
  "15q13.3dup", "7q11.23dup_distal"
)
lab_fc <- df_fc[CNV_locus %in% main_text_loci]

missing_loci <- setdiff(main_text_loci, df_fc$CNV_locus)
if (length(missing_loci) > 0) {
  warning(
    "The following main-text loci were not found in df_fc (likely ",
    "zero frequency in one cohort, so excluded from this non-zero-only ",
    "plot) and will NOT be labelled: ", paste(missing_loci, collapse = ", ")
  )
}
cat("Labelled loci (n =", nrow(lab_fc), "of", length(main_text_loci), "expected):\n")
print(lab_fc[, .(CNV_locus, log2_fc)])

FC_BASE_SIZE <- 15

p_fc <- ggplot(df_fc, aes(x = reorder(CNV_locus, log2_fc), y = log2_fc)) +
  geom_hline(yintercept = 0, linewidth = 0.6) +
  geom_hline(yintercept = c(-1, 1), linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 2.8, alpha = 0.9) +
  geom_text_repel(
    data = lab_fc,
    aes(label = CNV_locus),
    size = 4.2,
    max.overlaps = Inf,
    box.padding = 0.25,
    min.segment.length = 0
  ) +
  coord_flip() +
  theme_classic(base_size = FC_BASE_SIZE) +
  theme(
    axis.line = element_line(linewidth = 0.6),
    axis.text.y = element_text(size = FC_BASE_SIZE),
    axis.text.x = element_text(size = FC_BASE_SIZE),
    axis.title.x = element_text(size = FC_BASE_SIZE + 1)
  ) +
  labs(
    x = NULL,
    y = expression("log"[2]*" fold-change")
  )

save_plot_both(p_fc, file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_log2_foldchange_nonzero"),
               width = 8.8, height = 11.5)

# -----------------------------
# (C) Bar chart: ALL 54 loci (includes zeros)
# Displays raw measured frequency, not fold-change; rounded
# columns are fine here.
# -----------------------------
bar_dt <- rbind(
  df[, .(CNV_locus, Cohort = "UKB", freq = freq_percent_ukb)],
  df[, .(CNV_locus, Cohort = "CKB", freq = ckb_freq_percent)]
)

ord <- df[, .(CNV_locus, maxfreq = pmax(freq_percent_ukb, ckb_freq_percent))][order(-maxfreq)]
bar_dt[, CNV_locus := factor(CNV_locus, levels = ord$CNV_locus)]

cohort_cols <- c("UKB" = "#2E86AB", "CKB" = "#E67E22")  # blue / orange

BAR_BASE_SIZE <- 14  # both bar charts

p_bar <- ggplot(bar_dt, aes(x = CNV_locus, y = freq, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  scale_fill_manual(values = cohort_cols) +
  theme_classic(base_size = BAR_BASE_SIZE) +
  theme(
    axis.text.x = element_text(size = BAR_BASE_SIZE - 2, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = BAR_BASE_SIZE - 1),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "top",
    legend.text = element_text(size = BAR_BASE_SIZE)
  ) +
  labs(
    x = NULL,
    y = "Carrier frequency (%)",
    fill = NULL
  )

save_plot_both(p_bar, file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_bar_all54"),
               width = 13, height = 7.5)

# -----------------------------
# Bar chart excluding 2q13dup(NPHP1) for better scale
# -----------------------------
bar_dt_no <- bar_dt[CNV_locus != "2q13dup(NPHP1)"]
ord_no <- df[CNV_locus != "2q13dup(NPHP1)",
             .(CNV_locus, maxfreq = pmax(freq_percent_ukb, ckb_freq_percent))][order(-maxfreq)]
bar_dt_no[, CNV_locus := factor(CNV_locus, levels = ord_no$CNV_locus)]

p_bar_no <- ggplot(bar_dt_no, aes(x = CNV_locus, y = freq, fill = Cohort)) +
  geom_col(position = position_dodge(width = 0.65), width = 0.55) +
  scale_fill_manual(values = cohort_cols) +
  theme_classic(base_size = BAR_BASE_SIZE) +
  theme(
    axis.text.x = element_text(size = BAR_BASE_SIZE - 2, angle = 60, hjust = 1, vjust = 1),
    axis.text.y = element_text(size = BAR_BASE_SIZE - 1),
    axis.line = element_line(linewidth = 0.6),
    legend.position = "top",
    legend.text = element_text(size = BAR_BASE_SIZE)
  ) +
  labs(
    x = NULL,
    y = "Carrier frequency (%)",
    fill = NULL
  )

save_plot_both(p_bar_no, file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_bar_all54_no_NPHP1"),
               width = 13, height = 7.5)

cat("\nSaved figures:\n")
cat("  ", file.path(fig_main, "Fig_pathogenic_ckb_vs_ukb_loglog_scatter_nonzero.(png|pdf)"), "\n")
cat("  ", file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_log2_foldchange_nonzero.(png|pdf)"), "\n")
cat("  ", file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_bar_all54.(png|pdf)"), "\n")
cat("  ", file.path(fig_supp, "SFig_pathogenic_ckb_vs_ukb_bar_all54_no_NPHP1.(png|pdf)"), "\n\n")
cat("Completed: ", format(Sys.time()), "\n")
sink()
