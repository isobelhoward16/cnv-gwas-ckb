#!/usr/bin/env Rscript
# Stage 9: CNV calls around one stage 5 signal, packed into rows (Figures 5A and 5C, and the zooms of the other
# loci), in participants passing CNV QC, with protein-coding genes below.
# Usage: Rscript scripts/09_figures/05_figure5_zoom.R <chr> <model> <phenotype> <plot_type> <out_name> [lead_gene]
#   model     = stage 5 model whose CNVR (dashed lines) and lead probe (dotted line) are drawn
#   plot_type = deletion | duplication | both (which calls are drawn)
#   lead_gene = optional gene to highlight in red
# The CNVR uses the event boundaries where available, as in Table 1 (a single-probe CNVR is widened by 5 kb each
# side). Diagnostics print counts only. Output: results/09_figures/<version>/zooms/
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

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) stop("Usage: Rscript scripts/09_figures/05_figure5_zoom.R <chr> <model> <phenotype> <plot_type> <out_name> [lead_gene]")
chr_target <- as.integer(args[1]); model <- args[2]; phenotype_arg <- sub("_rint$", "", args[3]); plot_type <- args[4]
out_name <- args[5]; lead_gene <- if (length(args) >= 6) args[6] else ""
stopifnot(model %in% c("deletion", "duplication", "mirror"), plot_type %in% c("deletion", "duplication", "both"))


packed_row_half_height <- 0.16
base_size    <- 12
axis_text_size  <- 22
axis_title_size <- 30
show_legend  <- FALSE   # set to TRUE below when deletions and duplications are both drawn
show_lead_hits <- TRUE


gtf_path <- pa$gtf_gencode_v19

save_plot   <- TRUE
out_file    <- file.path(pa$results_root, "09_figures", ver, "zooms", paste0(out_name, ".png"))
plot_width  <- 14
plot_height <- 8.5
plot_dpi    <- 300

# =============================================================================
# LOAD CNV FILES
# =============================================================================

cnv_files <- list.files(cnv_dir, pattern = "\\.txt$", full.names = TRUE)
if (!length(cnv_files)) stop("No CNV files found in: ", cnv_dir)

cat("Loading CNV files:", length(cnv_files), "\n")
cnv <- rbindlist(lapply(cnv_files, fread), fill = TRUE)

if (ncol(cnv) < 17) stop("Unexpected CNV table format. Expected >=17 columns, got: ", ncol(cnv))

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

cnv <- cnv[!is.na(chr) & !is.na(start) & !is.na(end) & is.finite(start) & is.finite(end)]

# participants passing CNV QC (samples in the matrices that pass the master QC file)
mat_ids <- rownames(readRDS(file.path(pa$results_root, "02_matrices", ver, "combined_mirror_matrix_shared_probes_fixed.rds")))
qc_ids  <- fread(pa$qc_gwas_file)[LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100, GID]
n_before <- nrow(cnv); cnv <- cnv[ckid %in% intersect(mat_ids, qc_ids)]
cat("Calls restricted to participants passing CNV QC:", n_before, "->", nrow(cnv), "\n")

# CNVR and lead from stage 5 (model and phenotype from the command line)
x5 <- fread(file.path(pa$results_root, "05_signals", ver, "all", "cnvr_outputs", "cnvr_by_model_supplemental_unpruned_related.csv"))
r5 <- x5[chr == chr_target & cnv_type == model & phenotype == paste0(phenotype_arg, "_rint")]
if (nrow(r5) != 1) stop("Expected one stage 5 row for chr", chr_target, " ", model, " ", phenotype_arg, ", found ", nrow(r5))
s5 <- if (is.na(r5$auwerx_event_start)) r5$event_start else r5$auwerx_event_start
e5 <- if (is.na(r5$auwerx_event_end))   r5$event_end   else r5$auwerx_event_end
if (e5 - s5 < 1000) { s5 <- s5 - 5000; e5 <- e5 + 5000 }   # single-probe CNVR: widen so the dashed lines are visible
cat("Signal:", r5$locus_cytoband_final, phenotype_arg, model, "| CNVR (dashed):", s5, "-", e5, "| lead:", r5$lead_pos, r5$lead_probe,
    "| P", signif(r5$pval, 3), "| carriers del/dup", r5$n_del, "/", r5$n_dup, "\n")
focal_start <- s5; focal_end <- e5; lead_probe_pos <- r5$lead_pos
lead_hits <- data.table(label = paste(phenotype_arg, "lead"), pos = lead_probe_pos)

idx_swap <- which(cnv$start > cnv$end)
if (length(idx_swap)) {
  tmp <- cnv$start[idx_swap]
  cnv$start[idx_swap] <- cnv$end[idx_swap]
  cnv$end[idx_swap] <- tmp
  cat("Swapped start/end for", length(idx_swap), "rows with inverted coordinates.\n")
}

# =============================================================================
# SUBSET TO LOCUS WINDOW
# =============================================================================

dt <- cnv[chr == chr_target & start <= focal_end & end >= focal_start]

if (plot_type == "duplication") {
  dt <- dt[type %in% c("dup", "duplication")]
} else if (plot_type == "deletion") {
  dt <- dt[type %in% c("del", "deletion")]
} else if (plot_type == "both") {
  dt <- dt[type %in% c("dup", "duplication", "del", "deletion")]
} else {
  stop("plot_type must be one of: duplication, deletion, both")
}

if (!nrow(dt)) stop("No CNVs found for chr", chr_target, " overlapping focal interval ", focal_start, "-", focal_end)

dt <- unique(dt, by = c("ckid", "start", "end", "type", "source_file"))
dt[, width_bp := end - start + 1]

dt[, cnv_type := fifelse(type %in% c("dup", "duplication"), "duplication",
                  fifelse(type %in% c("del", "deletion"), "deletion", NA_character_))]
dt <- dt[!is.na(cnv_type)]
if (!nrow(dt)) stop("No rows remaining after cnv_type standardization.")

# =============================================================================
# CLASSIFICATION
# =============================================================================

overlaps_interval <- function(s, e, a, b) (s <= b) & (e >= a)

dt[, overlaps_focal := overlaps_interval(start, end, focal_start, focal_end)]
dt[, class := fifelse(overlaps_focal & start < focal_start & end > focal_end, "spans_focal",
               fifelse(overlaps_focal, "overlaps_focal", "flanking"))]
dt[, class := factor(class, levels = c("overlaps_focal", "spans_focal", "flanking"))]

# =============================================================================
# PACK ROWS
# =============================================================================

dt_dup <- copy(dt[cnv_type == "duplication"])
dt_del <- copy(dt[cnv_type == "deletion"])

row_spacing <- 0.5  # try 0.8, 0.7, 0.6

if (nrow(dt_dup)) setorder(dt_dup, class, start, end, ckid)
if (nrow(dt_del)) setorder(dt_del, class, start, end, ckid)

if (nrow(dt_dup)) {
  dt_dup[, y := seq_len(.N) * row_spacing]
  dt_dup[, `:=`(ymin = y - packed_row_half_height, ymax = y + packed_row_half_height)]
}
if (nrow(dt_del)) {
  dt_del[, y := -seq_len(.N) * row_spacing]
  dt_del[, `:=`(ymin = y - packed_row_half_height, ymax = y + packed_row_half_height)]
}

dt_plot <- rbindlist(list(dt_dup, dt_del), use.names = TRUE, fill = TRUE)
if (!nrow(dt_plot)) stop("No plotted rows after packing.")

# =============================================================================
# PLOT RANGE
# =============================================================================

x_min <- floor(min(dt_plot$start, na.rm = TRUE))
x_max <- ceiling(max(dt_plot$end, na.rm = TRUE))
x_pad <- max(1000, round((x_max - x_min) * 0.01))
x_lim_min <- x_min - x_pad
x_lim_max <- x_max + x_pad
# cap the window at the CNVR +/- max(250 kb, CNVR width), so a few very long calls cannot stretch it; calls are
# drawn clipped to the window (their full extent is in the callset)
win_pad <- max(250e3, focal_end - focal_start)
x_lim_min <- max(x_lim_min, focal_start - win_pad); x_lim_max <- min(x_lim_max, focal_end + win_pad)
n_clip <- dt_plot[start < x_lim_min | end > x_lim_max, .N]
dt_plot[, `:=`(start = pmax(start, x_lim_min), end = pmin(end, x_lim_max))]
cat("Window:", x_lim_min, "-", x_lim_max, "| calls clipped at the window edge:", n_clip, "\n")
if (plot_type == "both") show_legend <- TRUE

y_top <- if (nrow(dt_dup)) max(dt_dup$ymax, na.rm = TRUE) else 0.5
y_bot <- if (nrow(dt_del)) min(dt_del$ymin, na.rm = TRUE) else -0.5
pad_top <- 0.6
pad_bot <- 0.35
y_abs <- max(abs(y_top), abs(y_bot))

if (show_lead_hits && nrow(lead_hits)) {
  lead_hits <- copy(lead_hits)
  lead_hits <- lead_hits[!is.na(pos)]
  lead_hits[, pos := as.numeric(pos)]
  lead_hits <- lead_hits[pos >= x_lim_min & pos <= x_lim_max]
}

# =============================================================================
# GENE TRACK
# =============================================================================

gtf_raw <- fread(
  cmd = paste("zcat", gtf_path, "| grep -v '^#'"),
  sep = "\t", header = FALSE, quote = ""
)
setnames(gtf_raw, c("seqname","source","feature","start","end","score","strand","frame","attribute"))

gtf_raw <- gtf_raw[feature == "gene" & grepl("protein_coding", attribute)]
gtf_raw[, seqname := gsub("^chr", "", seqname)]
gtf_raw <- gtf_raw[seqname == as.character(chr_target)]
gtf_raw[, start := as.numeric(start)]
gtf_raw[, end   := as.numeric(end)]
gtf_raw[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attribute)]

genes <- gtf_raw[(start + end) / 2 >= x_lim_min & (start + end) / 2 <= x_lim_max, .(gene_name, start, end, strand)]   # gene midpoint in window
genes <- unique(genes, by = "gene_name")
genes <- genes[toupper(gene_name) != "LYRM5"]
setorder(genes, start)

cat("\nProtein-coding genes in window:", nrow(genes), "\n")
if (nrow(genes)) print(genes[, .(gene_name, start, end)])

# pack genes by the width of their italic labels (about 2% of the window per character at this text size),
# so labels never overlap; at most 4 label rows, genes that do not fit are drawn without a label
char_bp <- (x_lim_max - x_lim_min) * 0.02
pack_genes <- function(dt) {
  dt[, `:=`(row = NA_integer_, show_label = TRUE)]
  if (!nrow(dt)) { dt[, xlab := numeric(0)]; return(dt) }
  # label centre kept inside the window, so no label runs past the panel edge into the axis title
  dt[, xlab := pmin(pmax((start + end) / 2, x_lim_min + nchar(gene_name) * char_bp / 2), x_lim_max - nchar(gene_name) * char_bp / 2)]
  dt[, `:=`(l = pmin(start, xlab - nchar(gene_name) * char_bp / 2), r = pmax(end, xlab + nchar(gene_name) * char_bp / 2))]
  row_ends <- numeric(0)
  for (i in seq_len(nrow(dt))) {
    placed <- FALSE
    for (k in seq_along(row_ends)) if (dt$l[i] > row_ends[k]) { dt$row[i] <- k; row_ends[k] <- dt$r[i]; placed <- TRUE; break }
    if (!placed) {
      if (length(row_ends) < 4) { row_ends <- c(row_ends, dt$r[i]); dt$row[i] <- length(row_ends) }
      else { dt$row[i] <- 4L; dt$show_label[i] <- FALSE }
    }
  }
  dt[, c("l", "r") := NULL]
  dt
}
genes <- pack_genes(genes)
cat("Genes labelled:", sum(genes$show_label), "of", nrow(genes), "\n")

gs <- max(1, (y_top - y_bot) / 15)   # gene track scaled to the stack height
gene_row_height <- 2.8 * gs
gene_bar_height <- 1.2 * gs
gene_sep_y <- y_bot - (pad_bot + 0.3) * gs

genes[, ymid      := gene_sep_y - 0.6 * gs - (row - 1) * gene_row_height]
genes[, ymin_gene := ymid - gene_bar_height / 2]
genes[, ymax_gene := ymid + gene_bar_height / 2]
genes[, xmid      := (start + end) / 2]
genes[, is_lead_gene := nzchar(lead_gene) & toupper(gene_name) == toupper(lead_gene)]

gene_y_bot <- if (nrow(genes)) min(genes$ymin_gene, na.rm = TRUE) - 0.6 * gs else gene_sep_y - gs
y_lim <- c(gene_y_bot, y_top + pad_top * gs)

# =============================================================================
# DIAGNOSTICS
# =============================================================================

cat("\nSelected CNVs for plotting:", nrow(dt_plot), "\n")
cat("Unique carriers:", uniqueN(dt_plot$ckid), "\n")
cat("Observed plot range:", format(x_min, big.mark=","), "-", format(x_max, big.mark=","), "\n")
cat("\nCounts by CNV type:\n"); print(dt_plot[, .N, by = cnv_type][order(cnv_type)])
cat("\nCounts by ordering class:\n"); print(dt_plot[, .N, by = class][order(class)])

# =============================================================================
# BUILD PLOT
# =============================================================================

# tick step from the window width; labels with enough decimals to be exact
w_mb <- (x_lim_max - x_lim_min) / 1e6
step <- if (w_mb <= 1) 1e5 else if (w_mb <= 2.5) 2.5e5 else 5e5
mb_breaks <- seq(ceiling(x_lim_min / step) * step, floor(x_lim_max / step) * step, by = step)
mb_digits <- if (step %% 1e5 == 0) 1 else 2


fill_vals_type <- c(duplication = "#2C7FB8", deletion = "#D95F0E")
vline_lwd    <- 0.6
baseline_lwd <- 0.3

p <- ggplot(dt_plot) +
  geom_hline(yintercept = 0, linewidth = baseline_lwd, color = "black") +
  geom_vline(xintercept = c(focal_start, focal_end),
             linetype = "dashed", linewidth = vline_lwd, color = "black", alpha = 0.9) +
  geom_rect(aes(xmin = start, xmax = end, ymin = ymin, ymax = ymax, fill = cnv_type),
            color = NA, alpha = 0.90) +
  geom_hline(yintercept = gene_sep_y, linewidth = baseline_lwd, color = "black") +
  geom_rect(data = genes[is_lead_gene == FALSE],
            aes(xmin = start, xmax = end, ymin = ymin_gene, ymax = ymax_gene),
            inherit.aes = FALSE, fill = "grey30", color = "grey10", linewidth = 0.25, alpha = 0.85) +
  geom_rect(data = genes[is_lead_gene == TRUE],
            aes(xmin = start, xmax = end, ymin = ymin_gene, ymax = ymax_gene),
            inherit.aes = FALSE, fill = "#D7301F", color = "grey10", linewidth = 0.25, alpha = 0.55) +
  geom_text(data = genes[is_lead_gene == FALSE & show_label == TRUE],
            aes(x = xlab, y = ymin_gene - 0.12 * gs, label = gene_name),
            inherit.aes = FALSE, size = 8.4, fontface = "italic", vjust = 1, color = "black") +
  geom_text(data = genes[is_lead_gene == TRUE & show_label == TRUE],
            aes(x = xlab, y = ymin_gene - 0.12 * gs, label = gene_name),
            inherit.aes = FALSE, size = 8.4, fontface = "italic", vjust = 1, color = "#D7301F") +

  {
    if (show_lead_hits && exists("lead_hits") && nrow(lead_hits)) {
      geom_vline(
        data = lead_hits,
        aes(xintercept = pos),
        inherit.aes = FALSE,
        linetype = "dotted",
        linewidth = 0.8,
        color = "#D7301F",
        alpha = 0.95
      )
    } else NULL
  } +
  scale_fill_manual(values = fill_vals_type,
                    breaks = c("deletion", "duplication"),
                    labels = c("DEL", "DUP"), drop = FALSE) +
  scale_x_continuous(limits = c(x_lim_min, x_lim_max), breaks = mb_breaks,
                     labels = function(x) sprintf(paste0("%.", mb_digits, "f"), x / 1e6), expand = c(0, 0)) +
  coord_cartesian(ylim = y_lim, clip = "off") +
  labs(x = sprintf("Position on chr%d (Mb)", chr_target),
       y = c(duplication = "Duplication calls", deletion = "Deletion calls", both = "CNV calls")[[plot_type]],
       fill = if (show_legend) "CNV type" else NULL) +
  theme_void(base_size = base_size) +
  theme(
    axis.line.x  = element_line(color = "black", linewidth = 0.3),
    axis.ticks.x = element_line(color = "black", linewidth = 0.3),
    axis.ticks.length.x = unit(2, "pt"),
    axis.text.x  = element_text(size = axis_text_size, color = "black", margin = margin(t = 4)),
    axis.title.x = element_text(size = axis_title_size, face = "bold", color = "black", margin = margin(t = 10)),
    axis.title.y = element_text(size = axis_title_size, face = "bold", color = "black", angle = 90, margin = margin(r = 10)),
    legend.position  = if (show_legend) "right" else "none",
    legend.title     = element_text(size = 15),
    legend.text      = element_text(size = 14),
    legend.key.height = unit(20, "pt"),
    legend.key.width  = unit(24, "pt"),
    plot.margin = margin(8, 12, 8, 8)
  ) +
  guides(fill = guide_legend(override.aes = list(alpha = 1),
                             keywidth = unit(24, "pt"), keyheight = unit(20, "pt")))


# =============================================================================
# SAVE
# =============================================================================

if (isTRUE(save_plot)) {
  dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_file, plot = p, width = plot_width, height = plot_height, dpi = plot_dpi, bg = "white")
  cat("\nSaved plot to:\n", out_file, "\n")
}

# =============================================================================
# OPTIONAL DIAGNOSTICS
# =============================================================================

carrier_counts <- dt_plot[, .N, by = ckid][order(-N)]
cat("\nCalls plotted:", nrow(dt_plot), "| deletions", dt_plot[cnv_type == "deletion", .N], "| duplications", dt_plot[cnv_type == "duplication", .N], "\n")
cat("People with >1 selected CNV in plot window:", carrier_counts[N > 1, .N], "\n")   # counts only, no IDs
