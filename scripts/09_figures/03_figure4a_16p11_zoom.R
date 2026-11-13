#!/usr/bin/env Rscript
# Stage 9: Figure 4A, CNV calls across 16p11.2 (BP2-BP3 and BP4-BP5) in participants passing CNV QC, packed
# into rows (duplications above the baseline, deletions below). Dashed lines mark the stage 5 CNVRs: BP2-BP3 is
# the union of the 16p11.2 CNVRs with lead probe < 29.3 Mb, BP4-BP5 the union of those with lead probe > 29.3 Mb.
# Diagnostics print counts only. Output: results/09_figures/<version>/zooms/
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(grid)
})
source("R/config.R"); suppressPackageStartupMessages(library(Matrix))
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version

# =============================================================================
# CONFIG
# =============================================================================

cnv_dir <- file.path(pa$results_root, "01_callset", ver, "final", "calls")

chr_target <- 16L
win_start  <- 28000000
win_end    <- 30500000

plot_type <- "both"

# Proximal / distal anchors used for row ordering + filtering "other" calls
prox_start <- 28450000
prox_end   <- 29450000

dist_start <- 29550000
dist_end   <- 30350000

# Plot appearance
packed_row_half_height <- 0.38
base_size <- 12
show_legend <- TRUE
show_lead_hits <- FALSE
show_region_labels <- FALSE

drop_other_broad <- TRUE

lead_hits <- data.table(
  label = c("BMI lead", "Body fat% lead"),
  pos   = c(28865000, 30174024)
)

save_plot <- TRUE
out_file  <- file.path(pa$results_root, "09_figures", ver, "zooms", "16p11_stacked_pubstyle_twocolour.png")
plot_width  <- 14
plot_height <- 8
plot_dpi    <- 300

# =============================================================================
# LOAD CNV FILES
# =============================================================================

cnv_files <- list.files(cnv_dir, pattern = "\\.txt$", full.names = TRUE)
if (!length(cnv_files)) stop("No CNV files found in: ", cnv_dir)

cat("Loading CNV files:", length(cnv_files), "\n")
cnv <- rbindlist(lapply(cnv_files, fread), fill = TRUE)

if (ncol(cnv) < 17) {
  stop("Unexpected CNV table format. Expected >=17 columns, got: ", ncol(cnv))
}

setnames(cnv, c(
  "region","length","n_probes","source_file",
  "start_probe","end_probe","conf","state","sample_qc",
  "chr","start","end","type","col14","ckid","gene","col17"
))

cnv[, chr   := as.integer(gsub("^chr", "", as.character(chr)))]
cnv[, start := suppressWarnings(as.numeric(start))]
cnv[, end   := suppressWarnings(as.numeric(end))]
cnv[, ckid  := as.character(ckid)]
cnv[, type  := tolower(as.character(type))]
cnv[, source_file := as.character(source_file)]

cnv <- cnv[
  !is.na(chr) & !is.na(start) & !is.na(end) &
    is.finite(start) & is.finite(end)
]

# participants passing CNV QC (samples in the matrices that pass the QC file)
mat_ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))
qc_ids  <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
n_before <- nrow(cnv); cnv <- cnv[ckid %in% intersect(mat_ids, qc_ids)]
cat("Calls restricted to participants passing CNV QC:", n_before, "->", nrow(cnv), "\n")

# dashed lines: stage 5 CNVRs (event boundaries where available, as in Table 1)
x5 <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
x5 <- x5[chr == 16 & locus_cytoband_final == "16p11.2"]
x5[, `:=`(s = fifelse(is.na(auwerx_event_start), event_start, auwerx_event_start), e = fifelse(is.na(auwerx_event_end), event_end, auwerx_event_end))]
cat("Stage 5 16p11.2 rows:\n"); print(x5[order(lead_pos), .(phenotype, cnv_type, lead_pos, s, e)])
if (!nrow(x5[lead_pos < 29.3e6]) || !nrow(x5[lead_pos > 29.3e6])) stop("Expected 16p11.2 signals at both BP2-BP3 and BP4-BP5")
bp23_start <- x5[lead_pos < 29.3e6, min(s)]; bp23_end <- x5[lead_pos < 29.3e6, max(e)]
bp45_start <- x5[lead_pos > 29.3e6, min(s)]; bp45_end <- x5[lead_pos > 29.3e6, max(e)]
cat("BP2-BP3 lines:", bp23_start, "-", bp23_end, "| BP4-BP5 lines:", bp45_start, "-", bp45_end, "\n")

idx_swap <- which(cnv$start > cnv$end)
if (length(idx_swap)) {
  tmp <- cnv$start[idx_swap]
  cnv$start[idx_swap] <- cnv$end[idx_swap]
  cnv$end[idx_swap] <- tmp
  cat("Swapped start/end for", length(idx_swap), "rows with inverted coordinates.\n")
}

# =============================================================================
# SUBSET TO BROAD 16p11.2 WINDOW
# =============================================================================

dt <- cnv[
  chr == chr_target &
    start <= win_end &
    end >= win_start
]

if (plot_type == "duplication") {
  dt <- dt[type %in% c("dup", "duplication")]
} else if (plot_type == "deletion") {
  dt <- dt[type %in% c("del", "deletion")]
} else if (plot_type == "both") {
  dt <- dt[type %in% c("dup", "duplication", "del", "deletion")]
} else {
  stop("plot_type must be one of: duplication, deletion, both")
}

if (!nrow(dt)) {
  stop("No CNVs found for chr", chr_target, " in broad window ",
       win_start, "-", win_end, " for plot_type=", plot_type)
}

dt <- unique(dt, by = c("ckid", "start", "end", "type", "source_file"))
dt[, width_bp := end - start + 1]

dt[, cnv_type := fifelse(type %in% c("dup", "duplication"), "duplication",
                  fifelse(type %in% c("del", "deletion"), "deletion", NA_character_))]
dt <- dt[!is.na(cnv_type)]

if (!nrow(dt)) stop("No rows remaining after cnv_type standardization.")

# =============================================================================
# CLASSIFY BREAKPOINT PATTERN
# =============================================================================

overlaps_interval <- function(s, e, a, b) {
  (s <= b) & (e >= a)
}

dt[, overlaps_prox := overlaps_interval(start, end, prox_start, prox_end)]
dt[, overlaps_dist := overlaps_interval(start, end, dist_start, dist_end)]

dt[, class := fifelse(overlaps_prox & overlaps_dist, "bridges_both",
               fifelse(overlaps_prox & !overlaps_dist, "proximal_only",
               fifelse(!overlaps_prox & overlaps_dist, "distal_only",
                       "other_broad_16p11")))]

class_order <- c("proximal_only", "bridges_both", "distal_only", "other_broad_16p11")
dt[, class := factor(class, levels = class_order)]

if (isTRUE(drop_other_broad)) {
  n_before <- nrow(dt)
  dt <- dt[class != "other_broad_16p11"]
  cat("Dropped ", n_before - nrow(dt), " CNVs in class 'other_broad_16p11'.\n", sep = "")
}

if (!nrow(dt)) stop("No CNVs left after filtering; try drop_other_broad <- FALSE")

# =============================================================================
# PACK ROWS: DUPS ABOVE, DELS BELOW
# =============================================================================

dt_dup <- copy(dt[cnv_type == "duplication"])
dt_del <- copy(dt[cnv_type == "deletion"])

if (nrow(dt_dup)) setorder(dt_dup, class, start, end, ckid)
if (nrow(dt_del)) setorder(dt_del, class, start, end, ckid)

if (nrow(dt_dup)) {
  dt_dup[, y := seq_len(.N)]
  dt_dup[, `:=`(ymin = y - packed_row_half_height,
                ymax = y + packed_row_half_height)]
}

if (nrow(dt_del)) {
  dt_del[, y := -seq_len(.N)]
  dt_del[, `:=`(ymin = y - packed_row_half_height,
                ymax = y + packed_row_half_height)]
}

dt_plot <- rbindlist(list(dt_dup, dt_del), use.names = TRUE, fill = TRUE)
if (!nrow(dt_plot)) stop("No plotted rows after packing.")

# =============================================================================
# PLOT RANGE
# =============================================================================

x_min <- floor(min(dt_plot$start, na.rm = TRUE))
x_max <- ceiling(max(dt_plot$end,   na.rm = TRUE))
x_pad <- max(1000, round((x_max - x_min) * 0.01))
x_lim_min <- x_min - x_pad
x_lim_max <- x_max + x_pad

y_top <- if (nrow(dt_dup)) max(dt_dup$ymax, na.rm = TRUE) else 0.5
y_bot <- if (nrow(dt_del)) min(dt_del$ymin, na.rm = TRUE) else -0.5

pad_top <- 0.6
pad_bot <- 0.35
y_lim <- c(y_bot - pad_bot, y_top + pad_top)
y_abs <- max(abs(y_top), abs(y_bot))

if (show_lead_hits && nrow(lead_hits)) {
  lead_hits <- copy(lead_hits)
  lead_hits <- lead_hits[!is.na(pos)]
  lead_hits[, label := as.character(label)]
  lead_hits[, pos := as.numeric(pos)]
  lead_hits <- lead_hits[pos >= x_lim_min & pos <= x_lim_max]
  if (nrow(lead_hits)) {
    lead_hits[, track := seq_len(.N)]
    lead_hits[, y := y_abs + 0.55 + (track - 1) * 0.32]
  }
}

# =============================================================================
# DIAGNOSTICS
# =============================================================================

cat("\nSelected CNVs for plotting:", nrow(dt_plot), "\n")
cat("Unique carriers:", uniqueN(dt_plot$ckid), "\n")
cat("Observed plot range:", format(x_min, big.mark=","), "-", format(x_max, big.mark=","), "\n")
cat("\nEmpirical boundaries used for vertical lines:\n")
cat("  BP2-BP3:", format(bp23_start, big.mark=","), "-", format(bp23_end, big.mark=","),
    "(stage 5 CNVR union)\n")
cat("  BP4-BP5:", format(bp45_start, big.mark=","), "-", format(bp45_end, big.mark=","),
    "(stage 5 CNVR union)\n")

cat("\nCounts by CNV type:\n")
print(dt_plot[, .N, by = cnv_type][order(cnv_type)])

cat("\nCounts by breakpoint class:\n")
print(dt_plot[, .N, by = class][order(class)])

# =============================================================================
# BUILD PLOT
# =============================================================================

mb_breaks <- seq(
  floor(x_lim_min / 5e5) * 5e5,
  ceiling(x_lim_max / 5e5) * 5e5,
  by = 5e5
)

fill_vals_type <- c(
  duplication = "#2C7FB8",
  deletion    = "#D95F0E"
)

vline_lwd    <- 0.6
baseline_lwd <- 0.3

p <- ggplot(dt_plot) +
  geom_hline(yintercept = 0, linewidth = baseline_lwd, color = "black") +

  # CNVR boundaries as vertical dashed lines
  geom_vline(
    xintercept = c(bp23_start, bp23_end),
    linetype = "dashed", linewidth = vline_lwd, color = "#1f78b4", alpha = 0.9
  ) +
  geom_vline(
    xintercept = c(bp45_start, bp45_end),
    linetype = "dashed", linewidth = vline_lwd, color = "#e31a1c", alpha = 0.9
  ) +

  geom_rect(
    aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = cnv_type),
    color = NA, alpha = 0.90
  ) +

  {
    if (show_lead_hits && exists("lead_hits") && nrow(lead_hits)) {
      list(
        geom_segment(
          data = lead_hits,
          aes(x = pos, xend = pos, y = y - 0.12, yend = y + 0.12),
          inherit.aes = FALSE,
          linewidth = 0.35,
          color = "black"
        ),
        geom_text(
          data = lead_hits,
          aes(x = pos, y = y + 0.16, label = label),
          inherit.aes = FALSE,
          size = 3,
          vjust = 0,
          hjust = 0.5
        )
      )
    } else NULL
  } +

  scale_fill_manual(
    values = fill_vals_type,
    breaks = c("deletion", "duplication"),
    labels = c("DEL", "DUP"),
    drop = FALSE
  ) +
  scale_x_continuous(
    limits = c(x_lim_min, x_lim_max),
    breaks = mb_breaks,
    labels = function(x) sprintf("%.1f", x / 1e6),
    expand = c(0, 0)
  ) +
  coord_cartesian(ylim = y_lim, clip = "off") +

  labs(
    x    = "Position on chr16 (Mb)",
    y    = "CNV calls",
    fill = if (show_legend) "CNV type" else NULL
  ) +

  theme_void(base_size = base_size) +
  theme(
    axis.line.x          = element_line(color = "black", linewidth = 0.3),
    axis.ticks.x         = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length.x  = unit(2, "pt"),
    axis.text.x          = element_text(size = 17, color = "black", margin = margin(t = 4)),
    axis.title.x         = element_text(size = 18, face = "bold", color = "black",
                                        margin = margin(t = 10)),
    axis.title.y         = element_text(size = 18, face = "bold", color = "black", angle = 90, margin = margin(r = 10)),
    legend.position      = if (show_legend) "right" else "none",
    legend.title         = element_text(size = 15),
    legend.text          = element_text(size = 14),
    legend.key.height    = unit(20, "pt"),
    legend.key.width     = unit(24, "pt"),
    plot.margin          = margin(8, 12, 8, 8)
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(alpha = 1),
      keywidth  = unit(24, "pt"),
      keyheight = unit(20, "pt")
    )
  )


# =============================================================================
# SAVE
# =============================================================================

if (isTRUE(save_plot)) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, plot = p,
         width = plot_width, height = plot_height,
         dpi = plot_dpi, bg = "white")
  cat("\nSaved plot to:\n", out_file, "\n")
}

# =============================================================================
# DIAGNOSTICS
# =============================================================================

carrier_counts <- dt_plot[, .N, by = ckid][order(-N)]
cat("\nPeople with >1 selected CNV in plot window:", carrier_counts[N > 1, .N], "\n")
