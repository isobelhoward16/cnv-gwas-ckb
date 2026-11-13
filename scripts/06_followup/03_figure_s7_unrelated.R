#!/usr/bin/env Rscript
# Stage 6, step 3: scatter plot of effect estimates at the lead probes, all participants against
# unrelated participants.
# Output: results/06_followup/<version>/unrelated_leads/sensitivity_unrelated_scatter.{pdf,png}
suppressPackageStartupMessages({ library(ggplot2); library(dplyr) })
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version
od <- file.path(pa$results_root, "06_followup", ver, "unrelated_leads")
df <- fread(file.path(od, "table_s9_unrelated_sensitivity.csv"))
df[, cnv_type := model]

type_cols <- c("duplication" = "#E07B39", "deletion" = "#4A90D9", "mirror" = "#6AAB5E")

df <- df %>%
  mutate(cnv_type = factor(cnv_type, levels = c("duplication", "deletion", "mirror")))

lim <- max(abs(c(df$beta_all, df$beta_unrel)), na.rm = TRUE) * 1.15

p <- ggplot(df, aes(x = beta_all, y = beta_unrel,
                    colour = cnv_type, shape = cnv_type)) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", colour = "grey40", linewidth = 0.7) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
  geom_point(size = 4.5, alpha = 0.9) +
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
  coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(
    x = "Effect estimate, full sample",
    y = "Effect estimate, unrelated subsample"
  ) +
  theme_classic(base_size = 18) +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size = 16),
    legend.key.size = unit(1, "cm"),
    axis.text       = element_text(size = 16),
    axis.title      = element_text(size = 18),
    plot.margin     = margin(8, 12, 8, 12, "mm")
  )

ggsave(file.path(od, "sensitivity_unrelated_scatter.pdf"), p, width = 10, height = 10)
ggsave(file.path(od, "sensitivity_unrelated_scatter.png"), p, width = 10, height = 10, dpi = 200)
cat("Done:", file.path(od, "sensitivity_unrelated_scatter.pdf"), "and .png\n")
