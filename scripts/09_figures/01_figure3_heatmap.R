#!/usr/bin/env Rscript
# Stage 9: Figure 3, summary of the associations by locus and trait. Each cell shows the model(s) significant
# for that trait at that locus (deletion, duplication or mirror only; deletion or duplication with mirror), with
# deletion and duplication carrier frequencies at the lead probe of each locus in participants passing CNV QC.
# 16p11.2 is split into BP2-BP3 and BP4-BP5 by lead position (29.3 Mb); the stage 5 label "1q21.2" is the 1q21.1
# distal CNV, so 1q21.1 and 1q21.2 are one locus.
# Output: results/09_figures/<version>/Fig3_cytoband_association_summary_unpruned_no14q.{png,pdf}
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(grid)
  library(data.table)
  library(Matrix)
})

source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
input           <- file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv")
mirror_mat_path <- file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")
qc_master_path  <- pa$qc_gwas_file
outdir          <- file.path(pa$results_root, "09_figures", ver)
prefix          <- "Fig3_cytoband_association_summary_unpruned_no14q"
effect_col     <- "lead_beta"

plot_width       <- 20
plot_height      <- 14
bar_height_frac  <- 0.38
label_size       <- 13

if (!file.exists(input))          stop(paste0("ERROR: input file not found: ", input))
if (!file.exists(mirror_mat_path)) stop(paste0("ERROR: mirror matrix not found: ", mirror_mat_path))
if (!file.exists(qc_master_path))  stop(paste0("ERROR: QC master not found: ", qc_master_path))
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ---------- helpers ----------
parse_num <- function(x) {
  x <- as.character(x)
  x <- str_trim(x)
  x <- str_replace_all(x, "\u2212", "-")
  x <- str_replace_all(x, ",", "")
  suppressWarnings(as.numeric(x))
}

first_non_na  <- function(x) dplyr::first(stats::na.omit(x))
parse_chr     <- function(cyt) suppressWarnings(as.integer(str_extract(cyt, "^[0-9]+")))
parse_arm     <- function(cyt) str_match(cyt, "^[0-9]+([pq])")[, 2]
parse_bandnum <- function(cyt) suppressWarnings(as.numeric(
  str_match(cyt, "^[0-9]+[pq]([0-9]+(?:\\.[0-9]+)?)")[, 2]
))

# ---------- read results file ----------
hdr <- names(read.csv(input, nrows = 1, check.names = FALSE))
cc  <- rep("character", length(hdr)); names(cc) <- hdr
df  <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE, colClasses = cc)

req     <- c("phenotype","cnv_type","locus_cytoband_final","pval","n_del","n_dup","total_n", effect_col)
missing <- setdiff(req, names(df))
if (length(missing) > 0) stop(paste0("ERROR: Missing required columns: ", paste(missing, collapse = ", ")))

# ---------- clean core table ----------
df2 <- df %>%
  transmute(
    phenotype = str_remove(str_trim(phenotype), "_rint$"),
    cnv_type  = tolower(str_trim(cnv_type)),
    cytoband  = str_trim(locus_cytoband_final),
    lead_pos  = suppressWarnings(as.numeric(str_trim(lead_pos))),
    lead_probe = str_trim(lead_probe),
    pval      = parse_num(pval),
    beta      = parse_num(.data[[effect_col]]),
    n_del     = parse_num(n_del),
    n_dup     = parse_num(n_dup),
    total_n   = parse_num(total_n)
  ) %>%
  filter(
    !is.na(phenotype), phenotype != "",
    !is.na(cytoband),  cytoband  != "",
    !is.na(cnv_type),  cnv_type  != ""
  ) %>%
  distinct()

# ---------- exclude mirror rows with zero carriers of either type ----------
df2 <- df2 %>%
  filter(!(cnv_type == "mirror" & (is.na(n_del) | n_del == 0 | is.na(n_dup) | n_dup == 0)))

# ---------- split loci sharing a cytoband ----------
aw_start_map <- df %>%
  transmute(
    phenotype  = str_remove(str_trim(phenotype), "_rint$"),
    cnv_type   = tolower(str_trim(cnv_type)),
    lead_pos   = suppressWarnings(as.numeric(str_trim(lead_pos))),
    lead_probe = str_trim(lead_probe),
    aw_start   = suppressWarnings(as.numeric(str_trim(auwerx_event_start)))
  ) %>%
  distinct()

q21_both <- any(df2$cytoband %in% c("1q21.1", "1q21.2") & df2$lead_pos < 145000000) &&
            any(df2$cytoband %in% c("1q21.1", "1q21.2") & df2$lead_pos > 145000000)
df2 <- df2 %>%
  mutate(cytoband = case_when(
    cytoband == "16p11.2" & lead_pos < 29300000 ~ "16p11.2 BP2-BP3",
    cytoband == "16p11.2" & lead_pos > 29300000 ~ "16p11.2 BP4-BP5",
    cytoband %in% c("1q21.1", "1q21.2") & q21_both & lead_pos < 145000000 ~ "1q21.1a",
    cytoband %in% c("1q21.1", "1q21.2") & q21_both & lead_pos > 145000000 ~ "1q21.1b",
    cytoband %in% c("1q21.1", "1q21.2") ~ "1q21.1",
    TRUE ~ cytoband
  ))

# ---------- phenotype display order + labels ----------
pheno_order <- c(
  "standing_height_mm", "sitting_height_mm",
  "weight_kg", "bmi_calc", "fat_percent",
  "waist_mm", "hip_mm",
  "dbp_mean", "heart_rate_mean", "random_glucose"
)

pheno_labels <- c(
  standing_height_mm = "Height (standing)",
  sitting_height_mm  = "Height (sitting)",
  weight_kg          = "Weight",
  bmi_calc           = "BMI",
  fat_percent        = "Body fat (%)",
  waist_mm           = "WC",
  hip_mm             = "HC",
  dbp_mean           = "DBP",
  heart_rate_mean    = "Resting heart rate",
  random_glucose     = "RPG"
)

df2 <- df2 %>% filter(phenotype %in% pheno_order, cytoband != "14q11.2")

# ---------- signal class of each trait x locus, from the stage 5 models ----------
signal_class_lookup <- df2 %>%
  group_by(phenotype, cytoband) %>%
  summarise(models = paste(sort(unique(cnv_type)), collapse = "+"), .groups = "drop") %>%
  mutate(signal_class = recode(models,
    "deletion" = "Deletion only", "duplication" = "Duplication only", "mirror" = "Mirror only",
    "deletion+mirror" = "Deletion + Mirror", "duplication+mirror" = "Duplication + Mirror", .default = NA_character_))
if (anyNA(signal_class_lookup$signal_class))
  stop("Model combination without a class/colour: ", paste(unique(signal_class_lookup$models[is.na(signal_class_lookup$signal_class)]), collapse = ", "))
cat("=== Signal classes (trait x locus):", nrow(signal_class_lookup), "cells,", n_distinct(signal_class_lookup$cytoband), "loci\n")
print(as.data.frame(signal_class_lookup %>% arrange(cytoband, phenotype)))
signal_class_lookup <- signal_class_lookup %>% select(-models)

# ---------- direction labels ----------
best_beta <- df2 %>%
  filter(!is.na(pval)) %>%
  group_by(phenotype, cytoband, cnv_type) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(phenotype, cytoband, cnv_type, beta)

dir_del <- best_beta %>%
  filter(cnv_type == "deletion") %>%
  select(phenotype, cytoband, beta_del = beta) %>%
  mutate(label_del = case_when(
    is.na(beta_del) ~ NA_character_,
    beta_del > 0    ~ "+",
    beta_del < 0    ~ "\u2212",
    TRUE            ~ "0"
  ))

dir_dup <- best_beta %>%
  filter(cnv_type == "duplication") %>%
  select(phenotype, cytoband, beta_dup = beta) %>%
  mutate(label_dup = case_when(
    is.na(beta_dup) ~ NA_character_,
    beta_dup > 0    ~ "+",
    beta_dup < 0    ~ "\u2212",
    TRUE            ~ "0"
  ))

dir_mirror <- best_beta %>%
  filter(cnv_type == "mirror") %>%
  select(phenotype, cytoband, beta_mirror = beta) %>%
  mutate(label_mirror = case_when(
    is.na(beta_mirror) ~ NA_character_,
    beta_mirror > 0    ~ "+",
    beta_mirror < 0    ~ "\u2212",
    TRUE               ~ "0"
  ))

plot_df <- signal_class_lookup %>%
  left_join(dir_del,    by = c("phenotype", "cytoband")) %>%
  left_join(dir_dup,    by = c("phenotype", "cytoband")) %>%
  left_join(dir_mirror, by = c("phenotype", "cytoband")) %>%
  mutate(
    direction_label = case_when(
      signal_class == "Deletion only"                        ~ label_del,
      signal_class == "Duplication only"                     ~ label_dup,
      signal_class == "Mirror only"                          ~ label_mirror,
      signal_class %in% c("Duplication + Mirror",
                          "Deletion + Mirror")               ~ label_mirror,
      TRUE ~ NA_character_
    )
  )

# ---------- order cytobands ----------
plot_df <- plot_df %>%
  mutate(
    cytoband_base = str_remove(cytoband, "[ab]$"),
    chr           = parse_chr(cytoband_base),
    arm           = parse_arm(cytoband_base),
    band          = parse_bandnum(cytoband_base),
    arm_rank      = case_when(arm == "p" ~ 0, arm == "q" ~ 1, TRUE ~ 2),
    ab_rank       = case_when(
      str_ends(cytoband, "a") ~ 0,
      str_ends(cytoband, "b") ~ 1,
      TRUE ~ 0
    )
  ) %>%
  arrange(chr, arm_rank, band, ab_rank, cytoband)

cyto_levels <- unique(plot_df$cytoband)

plot_df <- plot_df %>%
  mutate(
    cytoband  = factor(cytoband,  levels = cyto_levels),
    phenotype = factor(phenotype, levels = rev(pheno_order))
  )

# ---------- chromosome boundary lines ----------
cyto_chr_map <- tibble(
  cytoband = cyto_levels,
  chr      = parse_chr(str_remove(cyto_levels, "[ab]$"))
) %>%
  mutate(idx = row_number()) %>%
  group_by(chr) %>%
  summarise(last_idx = max(idx), .groups = "drop") %>%
  filter(!is.na(chr)) %>%
  arrange(chr)

vline_pos <- cyto_chr_map$last_idx + 0.5
if (length(vline_pos) > 0) vline_pos <- vline_pos[-length(vline_pos)]

# ============================================================
# CARRIER FREQUENCIES FROM MIRROR MATRIX
# Denominator: all samples in the matrices passing CNV QC (with or without phenotype data)
# ============================================================

# Lead probes per locus - most significant across all models
lead_probe_lookup <- df2 %>%
  filter(!is.na(pval)) %>%
  group_by(cytoband) %>%
  slice_min(order_by = pval, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cytoband, lead_probe)

cat("Loading QC master file...\n")
qc_master <- fread(qc_master_path)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_ids <- qc_master[qc_pass == TRUE]$GID
cat("QC-passing samples:", length(qc_pass_ids), "\n")

cat("Loading mirror matrix...\n")
mirror_mat <- readRDS(mirror_mat_path)
cat("Mirror matrix dimensions:", nrow(mirror_mat), "x", ncol(mirror_mat), "\n")

common_ids  <- intersect(rownames(mirror_mat), qc_pass_ids)
cat("Samples in mirror matrix passing QC:", length(common_ids), "\n\n")
mirror_qc   <- mirror_mat[common_ids, , drop = FALSE]
n_total     <- length(common_ids)

lead_by_cytoband <- lead_probe_lookup %>%
  rowwise() %>%
  mutate(
    del_freq_pct = if (lead_probe %in% colnames(mirror_qc)) {
      vec <- as.numeric(mirror_qc[, lead_probe])
      round(100 * sum(vec < 0, na.rm = TRUE) / n_total, 4)
    } else NA_real_,
    dup_freq_pct = if (lead_probe %in% colnames(mirror_qc)) {
      vec <- as.numeric(mirror_qc[, lead_probe])
      round(100 * sum(vec > 0, na.rm = TRUE) / n_total, 4)
    } else NA_real_
  ) %>%
  ungroup() %>%
  mutate(cytoband = factor(cytoband, levels = cyto_levels))

cat("=== Carrier frequencies per locus ===\n")
print(as.data.frame(lead_by_cytoband))

freq_long <- lead_by_cytoband %>%
  pivot_longer(
    cols      = c(del_freq_pct, dup_freq_pct),
    names_to  = "cnv",
    values_to = "freq_pct"
  ) %>%
  mutate(cnv = recode(cnv,
                      del_freq_pct = "Deletion",
                      dup_freq_pct = "Duplication")) %>%
  tidyr::complete(
    cytoband,
    cnv = c("Deletion", "Duplication"),
    fill = list(freq_pct = 0)
  )

# ---------- colours ----------
fill_vals <- c(
  "Deletion + Mirror"    = "#fc8d62",
  "Deletion only"        = "#d95f02",
  "Duplication + Mirror" = "#66c2a5",
  "Duplication only"     = "#1b9e77",
  "Mirror only"          = "#7570b3"
)

freq_cols <- c(
  "Deletion"    = fill_vals[["Deletion only"]],
  "Duplication" = fill_vals[["Duplication only"]]
)

# ---------- shared theme ----------
# Saved on a 20 x 14 inch canvas and shown at about a third of that width, so text, tick lengths, legend keys
# and line widths are all scaled up together.
base_theme <- theme_classic(base_size = 28) +
  theme(
    text                 = element_text(family = "Nimbus Sans"),
    axis.line.x          = element_line(color = "black", linewidth = 1.4),
    axis.line.y          = element_line(color = "black", linewidth = 1.4),
    axis.ticks           = element_line(color = "black", linewidth = 1.3),
    axis.ticks.length    = unit(6, "pt"),
    legend.position      = "right",
    legend.box           = "vertical",
    legend.justification = "top",
    legend.title         = element_text(size = 26, face = "bold", family = "Nimbus Sans"),
    legend.text          = element_text(size = 24, family = "Nimbus Sans"),
    legend.key.size      = unit(1.8, "lines"),
    legend.spacing.y     = unit(0.3, "cm"),
    legend.margin        = margin(t = 0, b = 0),
    legend.box.spacing   = unit(0.5, "cm")
  )

# ---------- top frequency bar panel ----------
p_top <- ggplot(freq_long, aes(x = cytoband, y = freq_pct, fill = cnv)) +
  geom_col(width = 0.82,
           position = position_dodge2(width = 0.9, preserve = "single")) +
  { if (length(vline_pos) > 0)
      geom_vline(xintercept = vline_pos, linewidth = 0.25, color = "grey80",
                 linetype = "solid") } +
  scale_fill_manual(
    values = freq_cols,
    name   = "Carrier frequency",
    guide  = guide_legend(nrow = 2, title.position = "top")
  ) +
  # square-root scale, so frequencies of 0.02-0.05% are visible next to the 3% at 12p13.31
  scale_y_sqrt(breaks = c(0.05, 0.25, 1, 2, 3), labels = c("0.05", "0.25", "1", "2", "3"), expand = expansion(mult = c(0, 0.06))) +
  base_theme +
  theme(
    axis.title.x         = element_blank(),
    axis.text.x          = element_blank(),
    axis.ticks.x         = element_blank(),
    legend.position      = "right",
    legend.justification = "centre",
    plot.margin          = margin(8, 10, 0, 10)
  ) +
  labs(y = "Frequency (%)\n(square-root scale)")

# ---------- heatmap panel ----------
p_heat <- ggplot(plot_df, aes(x = cytoband, y = phenotype, fill = signal_class)) +
  geom_tile(color = NA) +
  geom_text(
    aes(label = direction_label),
    color    = "black",
    size     = label_size,
    fontface = "bold",
    na.rm    = TRUE
  ) +
  { if (length(vline_pos) > 0)
      geom_vline(xintercept = vline_pos, linewidth = 0.25, color = "grey80",
                 linetype = "solid") } +
  scale_fill_manual(
    values = fill_vals,
    breaks = names(fill_vals),
    drop   = FALSE,
    name   = "Model",
    guide  = guide_legend(ncol = 1, title.position = "top")
  ) +
  scale_y_discrete(labels = pheno_labels) +
  base_theme +
  theme(
    axis.text.x          = element_text(angle = 45, hjust = 1, vjust = 1),
    legend.position      = "right",
    legend.justification = "centre",
    plot.margin          = margin(4, 10, 8, 10)
  ) +
  labs(x = "Cytoband", y = "Trait")

# ---------- combine ----------
final <- p_top / p_heat +
  plot_layout(heights = c(bar_height_frac, 1 - bar_height_frac))

# ---------- save ----------
png_out <- file.path(outdir, paste0(prefix, ".png"))
pdf_out <- file.path(outdir, paste0(prefix, ".pdf"))

ggsave(png_out, final, width = plot_width, height = plot_height, dpi = 600)
ggsave(pdf_out, final, width = plot_width, height = plot_height, device = cairo_pdf)

message("Saved:\n  ", png_out, "\n  ", pdf_out)
