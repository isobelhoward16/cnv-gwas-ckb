#!/usr/bin/env Rscript
# Stage 9: Figures S3-S5, from the stage 6 array-stratified checks: effect estimates by array and pooled
# (Figure S3), carrier counts at each lead probe by array (Figure S4), and CNV x array interaction P against the
# number of CKB1 carriers (Figure S5; Bonferroni line for the number of interaction tests). Also writes a
# scatter of CKB1 against CKB2 estimates (pooling_scatter), not used in the paper.
# Output: results/09_figures/<version>/array_checks/
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(ggrepel)
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version

# ============================================================
# INPUT
# ============================================================
csv_path <- file.path(pa$results_root, "06_followup", ver, "array_checks", "pooling_validation_unpruned_mincar1_related.csv")
if (!file.exists(csv_path)) stop("CSV not found: ", csv_path)
cat("Using:", csv_path, "\n")

df <- read.csv(csv_path, stringsAsFactors = FALSE)
df <- df[, !names(df) %in% "."]

if (!"locus_cytoband_final" %in% names(df)) {
  stop("Expected column 'locus_cytoband_final' not found in validation CSV")
}
df$cytoband <- df$locus_cytoband_final
df$cytoband[df$cytoband == "1q21.2"] <- "1q21.1"   # sitting height: same 1q21.1 distal CNV as standing height

# ============================================================
# LABELS
# ============================================================
pheno_map <- c(
  bmi_calc_rint           = "BMI",
  dbp_mean_rint           = "DBP",
  fat_percent_rint        = "Fat %",
  heart_rate_mean_rint    = "Heart rate",
  hip_mm_rint             = "Hip circ.",
  random_glucose_rint     = "Glucose",
  sitting_height_mm_rint  = "Sitting ht.",
  standing_height_mm_rint = "Standing ht.",
  waist_mm_rint           = "Waist circ.",
  weight_kg_rint          = "Weight"
)

type_map <- c(
  duplication = "Dup",
  deletion    = "Del",
  mirror      = "Mir"
)

df <- df %>%
  mutate(
    pheno_label = pheno_map[phenotype],
    type_label  = type_map[cnv_type],
    locus       = ifelse(!is.na(cytoband) & cytoband != "", cytoband,
                         sub("AX-", "", lead_probe)),
    label_base  = paste(pheno_label, type_label, locus),
    label       = paste(pheno_label, type_label, locus)
  )

# Disambiguate duplicate labels with probe ID suffix (temporary)
dup_labels <- df$label[duplicated(df$label) | duplicated(df$label, fromLast = TRUE)]
df$label[df$label %in% dup_labels] <- paste0(
  df$label[df$label %in% dup_labels], " (",
  sub("AX-", "", df$lead_probe[df$label %in% dup_labels]), ")"
)

# Fix 16p11.2 labels - replace probe ID suffix with BP2-BP3 or BP4-BP5
# Uses lead_pos cutpoint of 29.3Mb
df <- df %>%
  mutate(
    label = case_when(
      grepl("16p11.2", label) & lead_pos < 29300000 ~
        paste(pheno_label, type_label, "16p11.2 BP2-BP3"),
      grepl("16p11.2", label) & lead_pos >= 29300000 ~
        paste(pheno_label, type_label, "16p11.2 BP4-BP5"),
      TRUE ~ label
    )
  )

cat("Labels:\n")
print(df$label)

# ============================================================
# COLOURS / OUTPUT
# ============================================================
type_cols  <- c("duplication" = "#E07B39", "deletion" = "#4A90D9", "mirror" = "#6AAB5E")
phase_cols <- c("CKB1" = "#E07B39", "CKB2" = "#4A90D9", "Pooled" = "#1A1A1A")

out_dir <- file.path(pa$results_root, "09_figures", ver, "array_checks")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# PLOT 1: FOREST PLOT
# ============================================================
forest_long <- df %>%
  select(label, cnv_type, pool_beta, pool_se, p12_beta, p12_se, p3_beta, p3_se) %>%
  pivot_longer(
    cols = c(pool_beta, pool_se, p12_beta, p12_se, p3_beta, p3_se),
    names_to = c("phase", ".value"),
    names_pattern = "^(pool|p12|p3)_(beta|se)$"
  ) %>%
  mutate(
    phase_label = factor(
      case_when(
        phase == "pool" ~ "Pooled",
        phase == "p12"  ~ "CKB1",
        TRUE            ~ "CKB2"
      ),
      levels = c("CKB1", "CKB2", "Pooled")
    ),
    ci_lo    = beta - 1.96 * se,
    ci_hi    = beta + 1.96 * se,
    y_offset = case_when(
      phase == "p12" ~ -0.32,
      phase == "p3"  ~  0.32,
      TRUE           ~  0.0
    ),
    lw = if_else(phase == "pool", 1.2, 0.8)
  ) %>%
  filter(!is.na(beta))

# Row order: input row order reversed
df <- df %>% mutate(row_order = row_number())
label_order <- df %>% arrange(desc(row_order)) %>% pull(label)
label_order <- unique(label_order)

forest_long <- forest_long %>%
  mutate(label = factor(label, levels = label_order))

n <- nlevels(forest_long$label)
plot_height <- max(9, n * 0.55)

p_forest <- ggplot(forest_long,
                   aes(x = beta, y = as.numeric(label) + y_offset,
                       colour = phase_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60",
             linewidth = 0.5) +
  geom_errorbarh(aes(xmin = ci_lo, xmax = ci_hi, linewidth = I(lw)),
                 height = 0.12) +
  geom_point(aes(size = I(if_else(phase_label == "Pooled", 2.8, 2.0)))) +
  scale_colour_manual(values = phase_cols, name = NULL) +
  scale_y_continuous(
    breaks = seq_len(n),
    labels = levels(forest_long$label),
    expand = c(0.02, 0.02)
  ) +
  scale_x_continuous(breaks = seq(-1.5, 1.5, 0.5)) +
  labs(x = "Effect estimate (SD units, RINT)", y = NULL) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.y         = element_text(size = 12, colour = "black"),
    axis.text.x         = element_text(size = 12),
    axis.title.x        = element_text(size = 13),
    legend.position      = "bottom",
    legend.text          = element_text(size = 12),
    legend.key.size      = unit(0.8, "cm"),
    panel.grid.major.x   = element_line(colour = "grey93", linewidth = 0.3),
    plot.margin          = margin(6, 10, 6, 6, "mm")
  )

ggsave(file.path(out_dir, "pooling_forest.pdf"), p_forest,
       width = 10, height = plot_height)
ggsave(file.path(out_dir, "pooling_forest.png"), p_forest,
       width = 10, height = plot_height, dpi = 200)
cat("Forest plot saved\n")

# ============================================================
# PLOT 2: CKB1 vs CKB2 SCATTER
# Label top discordant points (|CKB1 beta - CKB2 beta| >= 0.8)
# Discordant signals concentrated in low-carrier-count 16p11.2
# dups - consistent with noise from small numbers in CKB1
# ============================================================
scatter_df <- df %>%
  filter(!is.na(p12_beta) & !is.na(p3_beta)) %>%
  mutate(
    ci_p12_lo       = p12_beta - 1.96 * p12_se,
    ci_p12_hi       = p12_beta + 1.96 * p12_se,
    ci_p3_lo        = p3_beta  - 1.96 * p3_se,
    ci_p3_hi        = p3_beta  + 1.96 * p3_se,
    cnv_type        = factor(cnv_type, levels = c("duplication", "deletion", "mirror")),
    concordance_gap = abs(p3_beta - p12_beta)
  )

concordance_label_threshold <- 0.8

lim <- max(
  abs(c(
    scatter_df$ci_p12_lo, scatter_df$ci_p12_hi,
    scatter_df$ci_p3_lo,  scatter_df$ci_p3_hi
  )),
  na.rm = TRUE
) * 1.05

scatter_label_df <- scatter_df %>%
  filter(concordance_gap >= concordance_label_threshold) %>%
  mutate(
    nudge_x     = ifelse(p12_beta >= 0,  0.22 * lim, -0.22 * lim),
    nudge_y     = ifelse(p3_beta  >= 0,  0.03 * lim, -0.03 * lim),
    label_hjust = ifelse(p12_beta >= 0, 0, 1)
  )

p_scatter <- ggplot(scatter_df,
                    aes(x = p12_beta, y = p3_beta,
                        colour = cnv_type, shape = cnv_type)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.4) +
  geom_errorbar(aes(ymin = ci_p3_lo, ymax = ci_p3_hi),
                width = 0, linewidth = 0.7, alpha = 0.5) +
  geom_errorbarh(aes(xmin = ci_p12_lo, xmax = ci_p12_hi),
                 height = 0, linewidth = 0.7, alpha = 0.5) +
  geom_point(size = 4.5) +
  geom_text_repel(
    data               = scatter_label_df,
    aes(label          = label),
    size               = 5,
    colour             = "black",
    min.segment.length = 0,
    segment.size       = 0.4,
    segment.colour     = "grey50",
    max.overlaps       = Inf,
    box.padding        = 0.4,
    point.padding      = 0.2,
    force              = 1.5,
    force_pull         = 0.05,
    nudge_x            = scatter_label_df$nudge_x,
    nudge_y            = scatter_label_df$nudge_y,
    hjust              = scatter_label_df$label_hjust,
    direction          = "y",
    show.legend        = FALSE
  ) +
  scale_colour_manual(
    values = type_cols,
    name   = "CNV type",
    labels = c("Duplication", "Deletion", "Mirror")
  ) +
  scale_shape_manual(
    values = c(duplication = 16, deletion = 15, mirror = 17),
    name   = "CNV type",
    labels = c("Duplication", "Deletion", "Mirror")
  ) +
  coord_fixed(
    xlim = c(-lim * 1.28, lim * 1.28),
    ylim = c(-lim, lim),
    clip = "off"
  ) +
  labs(x = "CKB1 effect estimate (RINT)",
       y = "CKB2 effect estimate (RINT)") +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18),
    plot.margin     = margin(8, 28, 8, 28, "mm")
  )

ggsave(file.path(out_dir, "pooling_scatter.pdf"), p_scatter,
       width = 13, height = 13)
ggsave(file.path(out_dir, "pooling_scatter.png"), p_scatter,
       width = 13, height = 13, dpi = 200)
cat("Scatter plot saved\n")

# ============================================================
# PLOT 3: CKB1 CARRIER COUNT vs INTERACTION P
# Label only top 4 most significant interactions to avoid clutter
# Point: nominally significant interactions concentrated in
# signals with low CKB1 carrier counts - noise not heterogeneity
# ============================================================
int_df <- df %>%
  filter(!is.na(int_p)) %>%
  mutate(
    log10p   = -log10(int_p),
    cnv_type = factor(cnv_type, levels = c("duplication", "deletion", "mirror"))
  )

# Bonferroni threshold for the number of interaction tests
bonf_n_tests  <- nrow(int_df)
bonf_threshold <- -log10(0.05 / bonf_n_tests)
noml_threshold <- -log10(0.05)

# Label only top 4 by interaction p-value
int_label_df <- int_df %>%
  slice_min(int_p, n = 4)
# invisible anchor points along the two threshold lines, so the labels are pushed off the lines;
# seed fixed so the label layout is reproducible
if (any(int_df$carriers_phase12_after_mask < 1)) stop("A lead has no CKB1 carriers: the log-scale x axis needs >= 1")
x_rng <- range(int_df$carriers_phase12_after_mask)
line_anchor_df <- expand.grid(carriers_phase12_after_mask = 10^seq(log10(x_rng[1]), log10(x_rng[2]), length.out = 200),
                              log10p = c(noml_threshold, bonf_threshold))
line_anchor_df$label <- ""
int_repel_df <- bind_rows(int_label_df %>% select(carriers_phase12_after_mask, log10p, label), line_anchor_df)

p_int <- ggplot(int_df,
                aes(x = carriers_phase12_after_mask, y = log10p,
                    colour = cnv_type, shape = cnv_type)) +
  geom_hline(yintercept = noml_threshold, linetype = "dashed",
             colour = "grey40", linewidth = 0.8) +
  geom_hline(yintercept = bonf_threshold, linetype = "dotted",
             colour = "black", linewidth = 0.8) +
  geom_point(size = 4.5, alpha = 0.85) +
  geom_text_repel(
    data               = int_repel_df,
    inherit.aes        = FALSE,
    mapping            = aes(x = carriers_phase12_after_mask, y = log10p, label = label),
    seed               = 42,
    nudge_x            = 0.12,     # log10 units; labels to the right of each point
    hjust              = 0,
    direction          = "y",
    size               = 5,
    colour             = "black",
    min.segment.length = 0,
    segment.size       = 0.4,
    segment.colour     = "grey50",
    max.overlaps       = Inf,
    box.padding        = 0.7,
    point.padding      = 0.5,
    force              = 2,
    force_pull         = 0.2,
    show.legend        = FALSE
  ) +
  annotate(
    "text",
    x = Inf, y = noml_threshold + 0.06,
    hjust = 1.05, vjust = 0,
    size = 5.5, colour = "grey40",
    label = "p = 0.05"
  ) +
  annotate(
    "text",
    x = Inf, y = bonf_threshold + 0.06,
    hjust = 1.05, vjust = 0,
    size = 5.5, colour = "black",
    label = paste0("Bonferroni (n=", bonf_n_tests, ")")
  ) +
  scale_colour_manual(
    values = type_cols,
    name   = "CNV type",
    labels = c("Duplication", "Deletion", "Mirror")
  ) +
  scale_shape_manual(
    values = c(duplication = 16, deletion = 15, mirror = 17),
    name   = "CNV type",
    labels = c("Duplication", "Deletion", "Mirror")
  ) +
  scale_x_log10(breaks = c(3, 10, 30, 100, 300, 1000), labels = c("3", "10", "30", "100", "300", "1,000")) +   # log scale (12p13.31 has over 1,000 CKB1 carriers)
  labs(
    x = "CKB1 carrier count at lead probe (log scale)",
    y = expression(-log[10](italic(p)[interaction]))
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18)
  )

ggsave(file.path(out_dir, "pooling_interaction.pdf"), p_int,
       width = 14, height = 10)
ggsave(file.path(out_dir, "pooling_interaction.png"), p_int,
       width = 14, height = 10, dpi = 200)
cat("Interaction plot saved\n")

# ============================================================
# PLOT 4: CARRIER COUNTS BY ARRAY FOR SIGNIFICANT SIGNALS
# ============================================================
sig_pool_threshold <- 0.05 / nrow(df)

carrier_bar_df <- df %>%
  filter(!is.na(pool_p) & pool_p < sig_pool_threshold) %>%
  select(label, cnv_type,
         carriers_phase12_after_mask, carriers_phase3_after_mask) %>%
  pivot_longer(
    cols      = c(carriers_phase12_after_mask, carriers_phase3_after_mask),
    names_to  = "phase",
    values_to = "carriers"
  ) %>%
  mutate(
    phase = recode(
      phase,
      carriers_phase12_after_mask = "CKB1",
      carriers_phase3_after_mask  = "CKB2"
    ),
    phase = factor(phase, levels = c("CKB1", "CKB2"))
  )

if (nrow(carrier_bar_df) > 0) {
  carrier_order <- carrier_bar_df %>%
    group_by(label) %>%
    summarise(total_carriers = sum(carriers, na.rm = TRUE), .groups = "drop") %>%
    arrange(total_carriers) %>%
    pull(label)

  carrier_bar_df <- carrier_bar_df %>%
    mutate(label = factor(label, levels = carrier_order))

  n_bars     <- length(carrier_order)
  bar_height <- max(7, n_bars * 0.6)

  p_carriers <- ggplot(carrier_bar_df,
                       aes(x = carriers, y = label, fill = phase)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.62,
             colour = NA) +
    scale_fill_manual(values = phase_cols, name = NULL) +
    labs(x = "Carrier count at lead probe after masking", y = NULL) +
    theme_classic(base_size = 13) +
    theme(
      legend.position = "bottom",
      legend.text     = element_text(size = 12),
      legend.key.size = unit(0.8, "cm"),
      axis.text.y     = element_text(size = 12, colour = "black"),
      axis.text.x     = element_text(size = 12),
      axis.title.x    = element_text(size = 13),
      plot.margin     = margin(6, 10, 6, 6, "mm")
    )

  ggsave(file.path(out_dir, "pooling_sig_carrier_counts.pdf"), p_carriers,
         width = 10, height = bar_height)
  ggsave(file.path(out_dir, "pooling_sig_carrier_counts.png"), p_carriers,
         width = 10, height = bar_height, dpi = 200)
  cat("Carrier-count bar plot saved\n")
} else {
  cat("No Bonferroni-significant pooled associations found; skipped\n")
}

cat("All plots written to:", out_dir, "\n")
