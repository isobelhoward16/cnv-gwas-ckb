#!/usr/bin/env Rscript
# Stage 8, step 3: Figure 2, genome-wide CNV density (deletions and duplications per bin; bin size from
# config/pipeline.yaml, 500 kb in the paper), with labels for the highest-frequency regions.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(cowplot)
  library(stringr)
  library(yaml)
  library(scales)
})

# -----------------------------
# Helpers
# -----------------------------
timestamp_string <- function() format(Sys.time(), "%Y-%m-%d_%H%M%S")
dir_create <- function(x) if (!dir.exists(x)) dir.create(x, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf(...), "\n")
save_plot_both <- function(p, out_prefix, width = 12, height = 4, dpi = 320) {
  ggsave(paste0(out_prefix, ".png"), plot = p, width = width, height = height, dpi = dpi, units = "in")
  ggsave(paste0(out_prefix, ".pdf"), plot = p, width = width, height = height, units = "in")
}
safe_range <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return("NA")
  sprintf("[%s, %s]", format(min(x), scientific = FALSE), format(max(x), scientific = FALSE))
}
normalize_cnv_type <- function(x) {
  x <- tolower(trimws(x))
  ifelse(x %in% c("del","deletion"), "del",
         ifelse(x %in% c("dup","duplication"), "dup", NA_character_))
}
parse_chr_int <- function(x) {
  x <- as.character(x)
  x <- str_replace(x, "^chr", "")
  suppressWarnings(as.integer(x))
}

# -----------------------------
# Config
# -----------------------------
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)
bin_size <- if (!is.null(paths$bin_size)) as.numeric(paths$bin_size) else 1e6
bin_lab  <- if (bin_size >= 1e6) paste0(bin_size / 1e6, "Mb") else paste0(bin_size / 1e3, "kb")

out_root             <- paths$out_root
gene_cleaned_merged_dir <- paths$gene_cleaned_merged_dir
qc_master_file          <- paths$qc_master_file
manifest_path <- file.path(out_root, "data/derived/manifests/postQC_denominator_samples.tsv")

tables_dir   <- file.path(out_root, "data/derived/tables")
fig_main_dir <- file.path(out_root, "results/figures/main")
fig_supp_dir <- file.path(out_root, "results/figures/supp")
logs_dir     <- file.path(out_root, "logs")

dir_create(tables_dir); dir_create(fig_main_dir); dir_create(fig_supp_dir); dir_create(logs_dir)

run_ts   <- timestamp_string()
log_file <- file.path(logs_dir, paste0("landscape_bin_freq_", bin_lab, "_", run_ts, ".log"))
zz <- file(log_file, open = "wt")
sink(zz, type = "output", split = TRUE)

msg("=== Paper1: %s bin CNV landscape (carrier-based) ===", bin_lab)
msg("paths.yaml: %s", paths_yaml)
msg("manifest:   %s", manifest_path)
msg("cnv_dir:    %s", gene_cleaned_merged_dir)
msg("qc_master:  %s", qc_master_file)
msg("log:        %s", log_file)

# -----------------------------
# Denominators
# -----------------------------
denom_dt <- fread(manifest_path)
if (!("phase" %in% names(denom_dt))) stop("Manifest missing 'phase' column: ", manifest_path)

msg("\n--- Denominators (pre-QC filter) ---")
msg("Manifest n = %d", nrow(denom_dt))
print(denom_dt[, .N, by = phase])

# ---- Apply QC filter from master QC file ----
# Mirrors inclusion criteria used in all GWAS analyses:
#   LRR_SD <= 0.35, |WF| <= 0.05, BAF_drift <= 0.01, NumCNV <= 100
msg("\nLoading QC master file: %s", qc_master_file)
qc_master <- fread(qc_master_file)
qc_master[, qc_pass := LRR_SD <= 0.35 & abs(WF) <= 0.05 & BAF_drift <= 0.01 & NumCNV <= 100]
qc_pass_ids <- qc_master[qc_pass == TRUE]$GID
msg("QC-passing samples in master file: %d", length(qc_pass_ids))
msg("QC-failing samples excluded: %d", sum(!qc_master$qc_pass))

before_n <- nrow(denom_dt)
denom_dt <- denom_dt[ckid %in% qc_pass_ids]
msg("Removed %d QC-failing samples from manifest denominator", before_n - nrow(denom_dt))
msg("Denominator N after QC filter: %d", nrow(denom_dt))

denom_by_phase    <- denom_dt[, .N, by = phase]
setnames(denom_by_phase, "N", "denom")
overall_denom     <- nrow(denom_dt)  # derived from QC-filtered manifest, not hardcoded

msg("\n--- Denominators after QC filter ---")
msg("Overall denominator: %d", overall_denom)
print(denom_by_phase)
# ---------------------------------------------

# -----------------------------
# CNV ingest (HARD-CODED columns)
# -----------------------------
CHR_COL   <- 10L
START_COL <- 11L
END_COL   <- 12L
TYPE_COL  <- 13L
CKID_COL  <- 15L

cnv_files <- list.files(gene_cleaned_merged_dir, pattern = "\\.txt$", full.names = TRUE)
msg("\nFound CNV files: %d", length(cnv_files))
if (length(cnv_files) == 0) stop("No CNV files found in gene_cleaned_merged_dir.")

peek <- fread(cnv_files[1], sep = "\t", header = FALSE, nrows = 5, fill = TRUE, quote = "")
msg("CNV file peek ncol (first file): %d", ncol(peek))
if (ncol(peek) < CKID_COL) stop("CNV file has fewer columns than expected (need at least 15).")

read_cnv_file <- function(f) {
  dt <- fread(f, sep = "\t", header = FALSE, fill = TRUE, quote = "", showProgress = FALSE)

  phase <- if (str_detect(basename(f), "Phase12")) "Phase12" else if (str_detect(basename(f), "Phase3")) "Phase3" else NA_character_
  region_phase <- str_replace(basename(f), "^merged_gene_mapped_", "")
  region_phase <- str_replace(region_phase, "\\.txt$", "")

  data.table(
    file = basename(f),
    region_phase = region_phase,
    phase = phase,
    chr = parse_chr_int(dt[[CHR_COL]]),
    start = suppressWarnings(as.integer(dt[[START_COL]])),
    end = suppressWarnings(as.integer(dt[[END_COL]])),
    cnv_type = normalize_cnv_type(as.character(dt[[TYPE_COL]])),
    ckid = as.character(dt[[CKID_COL]])
  )
}

msg("\n--- Loading CNVs ---")
cnv_dt <- rbindlist(lapply(cnv_files, read_cnv_file), use.names = TRUE, fill = TRUE)

# Restrict CNV records to QC-passing samples
before_cnv_n <- nrow(cnv_dt)
cnv_dt <- cnv_dt[ckid %in% qc_pass_ids]
msg("CNV records after restricting to QC-passing samples: %d (removed %d)",
    nrow(cnv_dt), before_cnv_n - nrow(cnv_dt))

# -----------------------------
# Sanity checks
# -----------------------------
msg("CNV records loaded (post-QC): %d", nrow(cnv_dt))
msg("Unique CKIDs in CNV records:  %d", uniqueN(cnv_dt$ckid))

msg("\n--- Sanity checks (parsed fields) ---")
msg("Phases present: %s",   paste(sort(unique(cnv_dt$phase)),    collapse = ", "))
msg("CNV types present: %s", paste(sort(unique(cnv_dt$cnv_type)), collapse = ", "))
msg("Chr range: %s",   safe_range(cnv_dt$chr))
msg("Start range: %s", safe_range(cnv_dt$start))
msg("End range:   %s", safe_range(cnv_dt$end))
msg("Records with start > end: %d", cnv_dt[start > end, .N])
msg("Records with NA chr: %d",      cnv_dt[is.na(chr), .N])
msg("Records with NA cnv_type: %d", cnv_dt[is.na(cnv_type), .N])

if (quantile(cnv_dt$end, 0.95, na.rm = TRUE) < 1e6) {
  stop("end positions look too small to be bp coordinates (95th percentile < 1e6). Check hard-coded column mapping.")
}

cnv_dt_filt <- cnv_dt[
  !is.na(chr) & chr >= 1 & chr <= 22 &
    !is.na(start) & !is.na(end) & start <= end &
    !is.na(cnv_type) &
    !is.na(ckid) & ckid != "" &
    !is.na(phase)
]

msg("\nRecords remaining after chr1-22 + valid coords/type: %d", nrow(cnv_dt_filt))
msg("Unique CKIDs after filtering: %d", uniqueN(cnv_dt_filt$ckid))

msg("\nCNV record counts by phase x type (records):")
print(cnv_dt_filt[, .N, by = .(phase, cnv_type)][order(phase, cnv_type)])

# -----------------------------
# bins (hg19 autosomes); size set above from config
# -----------------------------
chr_lengths_hg19 <- c(
  `1`=249250621, `2`=243199373, `3`=198022430, `4`=191154276, `5`=180915260,
  `6`=171115067, `7`=159138663, `8`=146364022, `9`=141213431, `10`=135534747,
  `11`=135006516, `12`=133851895, `13`=115169878, `14`=107349540, `15`=102531392,
  `16`=90354753, `17`=81195210, `18`=78077248, `19`=59128983, `20`=63025520,
  `21`=48129895, `22`=51304566
)


bins <- rbindlist(lapply(1:22, function(cc) {
  L      <- chr_lengths_hg19[as.character(cc)]
  starts <- seq(1, L, by = bin_size)
  ends   <- pmin(starts + bin_size - 1, L)
  data.table(chr = cc, bin_start = as.integer(starts), bin_end = as.integer(ends))
}))
bins[, bin_id := sprintf("chr%d:%d-%d", chr, bin_start, bin_end)]

chr_offsets <- data.table(chr = 1:22)
chr_offsets[, chr_len := as.numeric(chr_lengths_hg19[as.character(chr)])]
chr_offsets[, offset  := c(0, cumsum(chr_len)[-.N])]

bins <- merge(bins, chr_offsets[, .(chr, offset, chr_len)], by = "chr", all.x = TRUE)
bins[, genome_mid := offset + (bin_start + bin_end)/2]
bins[, plot_pos := .I]
setkey(bins, chr, bin_start, bin_end)

chr_ticks      <- chr_offsets[, .(chr, tick = offset + chr_len/2)]
chr_boundaries <- chr_offsets[, .(boundary = offset)]
chr_midpoints  <- bins[, .(plot_mid = mean(plot_pos)), by = chr][order(chr)]

# -----------------------------
# Overlap + carrier frequencies
# -----------------------------
cnv_int <- cnv_dt_filt[, .(phase, cnv_type, chr, start, end, ckid)]
setnames(cnv_int, c("start", "end"), c("cnv_start", "cnv_end"))
setkey(cnv_int, chr, cnv_start, cnv_end)

ov <- foverlaps(
  cnv_int, bins,
  by.x = c("chr", "cnv_start", "cnv_end"),
  by.y = c("chr", "bin_start", "bin_end"),
  type = "any", nomatch = 0L
)

ov_unique <- unique(ov[, .(phase, cnv_type, bin_id, ckid)])

msg("\n--- Sanity checks (overlaps) ---")
msg("Total CNV-bin overlaps (raw): %d", nrow(ov))
msg("Unique (phase,type,bin,CKID) overlaps: %d", nrow(ov_unique))

counts_phase   <- ov_unique[, .(carriers = .N), by = .(phase, cnv_type, bin_id)]
counts_overall <- ov_unique[, .(carriers = uniqueN(ckid)), by = .(cnv_type, bin_id)]
counts_overall[, phase := "Overall"]

counts_all <- rbindlist(list(
  counts_overall[, .(phase, cnv_type, bin_id, carriers)],
  counts_phase
), use.names = TRUE)

denom_map <- rbind(
  data.table(phase = "Overall", denom = overall_denom),
  denom_by_phase[, .(phase, denom)]
)

counts_all <- merge(counts_all, denom_map, by = "phase", all.x = TRUE)
counts_all[, freq := carriers / denom]

grid <- CJ(phase = unique(counts_all$phase), cnv_type = c("del","dup"), bin_id = bins$bin_id, unique = TRUE)
counts_all <- merge(grid, counts_all, by = c("phase","cnv_type","bin_id"), all.x = TRUE)
counts_all[is.na(carriers), carriers := 0L]
counts_all[is.na(freq),     freq     := 0]
counts_all <- merge(counts_all, bins[, .(bin_id, chr, bin_start, bin_end, genome_mid, plot_pos)],
                    by = "bin_id", all.x = TRUE)

msg("\n--- Sanity checks (top bins) ---")
top_overall <- counts_all[phase == "Overall"][order(-freq)][1:10]
msg("Top 10 bins by frequency (overall) - bin_id + freq only:")
print(top_overall[, .(bin_id, cnv_type, freq, carriers, denom)])

# -----------------------------
# Save tables
# -----------------------------
out_table_stamp <- file.path(tables_dir, paste0("cnv_bin_freq_", bin_lab, "_", run_ts, ".tsv"))
fwrite(counts_all, out_table_stamp, sep = "\t")
fwrite(counts_all[phase == "Overall"],  file.path(tables_dir, paste0("cnv_bin_freq_", bin_lab, "_overall.tsv")),   sep = "\t")
fwrite(counts_all[phase != "Overall"],  file.path(tables_dir, paste0("cnv_bin_freq_", bin_lab, "_by_phase.tsv")),  sep = "\t")

sanity_out <- file.path(tables_dir, paste0("cnv_landscape_sanity_summary_", run_ts, ".tsv"))
sanity_lines <- rbindlist(list(
  data.table(metric = "cnv_records_loaded",                  value = nrow(cnv_dt)),
  data.table(metric = "unique_ckids_pre_filter",             value = uniqueN(cnv_dt$ckid)),
  data.table(metric = "cnv_records_post_filter",             value = nrow(cnv_dt_filt)),
  data.table(metric = "unique_ckids_post_filter",            value = uniqueN(cnv_dt_filt$ckid)),
  data.table(metric = "overlaps_raw",                        value = nrow(ov)),
  data.table(metric = "overlaps_unique_phase_type_bin_ckid", value = nrow(ov_unique)),
  data.table(metric = "overall_denom_used",                  value = overall_denom),
  data.table(metric = "qc_pass_samples_in_master_file",      value = length(qc_pass_ids))
))
fwrite(sanity_lines, sanity_out, sep = "\t")

msg("\nSaved sanity summary: %s", sanity_out)
msg("Saved bin frequency table (timestamped): %s", out_table_stamp)

# -----------------------------
# Plots
# -----------------------------
make_bins_for_plot <- function(dt) {
  del <- dt[cnv_type == "del", .(
    Chromosome = as.character(chr),
    Bin_Start  = bin_start,
    Bin_End    = bin_end,
    plot_pos   = plot_pos,
    percentage = freq * 100,
    Phase      = phase
  )]
  dup <- dt[cnv_type == "dup", .(
    Chromosome = as.character(chr),
    Bin_Start  = bin_start,
    Bin_End    = bin_end,
    plot_pos   = plot_pos,
    percentage = freq * 100,
    Phase      = phase
  )]
  list(del = del, dup = dup)
}

chr_midpoints_style <- bins[, .(midpoint = mean(plot_pos)), by = chr][order(chr)]
chr_midpoints_style[, Chromosome := as.character(chr)]

color_blue_dark    <- "#2E4057"
color_blue_light   <- "#7A9CC6"
color_orange_dark  <- "#B8621B"
color_orange_light <- "#E6A968"

axis_theme_bigger <- theme(
  axis.title = element_text(size = 13),
  axis.text  = element_text(size = 13)
)

cytoband_file <- paths$cytoband_file
find_cytoband <- NULL
if (file.exists(cytoband_file)) {
  cytoband <- fread(cytoband_file, col.names = c("chr", "start", "end", "band", "stain"))
  cytoband[, chr_clean := suppressWarnings(as.integer(gsub("chr", "", chr)))]
  cytoband <- cytoband[!is.na(chr_clean) & chr_clean %in% 1:22]
  cytoband[, chr := chr_clean]
  cytoband[, chr_clean := NULL]

  find_cytoband <- function(chromosome, bin_start, bin_end) {
    band_info <- cytoband[chr == as.integer(chromosome) & start < bin_end & end > bin_start]
    if (nrow(band_info) > 0) {
      band_info[, overlap := pmin(end, bin_end) - pmax(start, bin_start)]
      return(paste0(chromosome, band_info$band[which.max(band_info$overlap)]))
    } else {
      return(paste0(chromosome, "q"))
    }
  }
} else {
  msg("Cytoband file not found (%s); proceeding without labels.", cytoband_file)
  find_cytoband <- function(chr, start, end) paste0(chr, "q")
}

avoid_label_overlap <- function(label_data, min_distance = 100 * 1e6 / bin_size) {   # 100 Mb, in bins
  if (nrow(label_data) == 0) return(label_data)
  label_data <- label_data[order(plot_pos)]
  label_data[, show_label := TRUE]
  if (nrow(label_data) > 1) {
    for (i in 2:nrow(label_data)) {
      prev_shown_idx <- which(label_data$show_label[1:(i-1)])
      if (length(prev_shown_idx) > 0) {
        prev_shown_idx <- max(prev_shown_idx)
        distance <- label_data$plot_pos[i] - label_data$plot_pos[prev_shown_idx]
        if (distance < min_distance) {
          if (label_data$max_freq[i] > label_data$max_freq[prev_shown_idx]) {
            label_data$show_label[prev_shown_idx] <- FALSE
          } else {
            label_data$show_label[i] <- FALSE
          }
        }
      }
    }
  }
  label_data[show_label == TRUE]
}

# (A) Overall mirror density (main)
plot_overall  <- counts_all[phase == "Overall"]
bins_overall  <- make_bins_for_plot(plot_overall)
del_bins      <- bins_overall$del
dup_bins      <- bins_overall$dup

del_bins[, chr_group := as.factor(as.numeric(Chromosome) %% 2)]
dup_bins[, chr_group := as.factor(as.numeric(Chromosome) %% 2)]

top_n_labels_overall <- 12L
high_del <- del_bins[percentage > 0][order(-percentage)][1:min(.N, top_n_labels_overall)]
if (nrow(high_del) > 0) high_del[, cytoband := mapply(find_cytoband, Chromosome, Bin_Start, Bin_End)]
high_dup <- dup_bins[percentage > 0][order(-percentage)][1:min(.N, top_n_labels_overall)]
if (nrow(high_dup) > 0) high_dup[, cytoband := mapply(find_cytoband, Chromosome, Bin_Start, Bin_End)]

del_label_data <- if (nrow(high_del) > 0) high_del[, .(
  plot_pos = mean(plot_pos), max_freq = max(percentage), bin_count = .N
), by = cytoband] else data.table(cytoband = character(0))

dup_label_data <- if (nrow(high_dup) > 0) high_dup[, .(
  plot_pos = mean(plot_pos), max_freq = max(percentage), bin_count = .N
), by = cytoband] else data.table(cytoband = character(0))

del_label_data_clean <- avoid_label_overlap(del_label_data)
dup_label_data_clean <- avoid_label_overlap(dup_label_data)

p_del <- ggplot(del_bins, aes(x = plot_pos, y = percentage, color = chr_group)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(color_blue_dark, color_blue_light)) +
  scale_x_continuous(breaks = chr_midpoints_style$midpoint, labels = chr_midpoints_style$Chromosome,
                     expand = c(0.02, 0.01)) +
  scale_y_continuous(name = "Deletion Density (%)", expand = expansion(mult = c(0.02, 0.1))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 20, r = 20, b = 5, l = 40),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black")
  ) + axis_theme_bigger

if (nrow(del_label_data_clean) > 0) {
  del_label_data_adjusted <- copy(del_label_data_clean)
  plot_range     <- range(del_bins$plot_pos)
  left_boundary  <- plot_range[1] + (plot_range[2] - plot_range[1]) * 0.02
  right_boundary <- plot_range[2] - (plot_range[2] - plot_range[1]) * 0.02
  del_label_data_adjusted[plot_pos < left_boundary,  plot_pos := left_boundary]
  del_label_data_adjusted[plot_pos > right_boundary, plot_pos := right_boundary]
  p_del <- p_del + geom_text(
    data = del_label_data_adjusted,
    aes(x = plot_pos, y = max_freq + max(del_bins$percentage) * 0.05, label = cytoband),
    color = "black", size = 4, hjust = 0.5, vjust = 0, inherit.aes = FALSE
  )
}

p_dup <- ggplot(dup_bins, aes(x = plot_pos, y = percentage, color = chr_group)) +
  geom_point(size = 0.5) +
  scale_color_manual(values = c(color_orange_dark, color_orange_light)) +
  scale_x_continuous(breaks = chr_midpoints_style$midpoint, labels = chr_midpoints_style$Chromosome,
                     expand = c(0.02, 0.01), position = "top") +
  scale_y_continuous(name = "Duplication Density (%)", trans = "reverse",
                     expand = expansion(mult = c(0.1, 0.02))) +
  theme_minimal() +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    panel.grid = element_blank(),
    plot.margin = margin(t = 5, r = 20, b = 20, l = 40),
    axis.line = element_line(color = "black", linewidth = 0.5),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(color = "black")
  ) + axis_theme_bigger

if (nrow(dup_label_data_clean) > 0) {
  dup_label_data_adjusted <- copy(dup_label_data_clean)
  plot_range     <- range(dup_bins$plot_pos)
  left_boundary  <- plot_range[1] + (plot_range[2] - plot_range[1]) * 0.02
  right_boundary <- plot_range[2] - (plot_range[2] - plot_range[1]) * 0.02
  dup_label_data_adjusted[plot_pos < left_boundary,  plot_pos := left_boundary]
  dup_label_data_adjusted[plot_pos > right_boundary, plot_pos := right_boundary]
  p_dup <- p_dup + geom_text(
    data = dup_label_data_adjusted,
    aes(x = plot_pos, y = max_freq + max(dup_bins$percentage) * 0.05, label = cytoband),
    color = "black", size = 4, hjust = 0.5, vjust = 1, inherit.aes = FALSE
  )
}

p_chr <- ggplot(data.frame(x = chr_midpoints_style$midpoint, y = 0, chr = chr_midpoints_style$Chromosome)) +
  geom_text(aes(x = x, y = y, label = chr), size = 3.5, hjust = 0.5, vjust = 0.5) +
  scale_x_continuous(limits = range(del_bins$plot_pos), expand = c(0.02, 0.01)) +
  scale_y_continuous(limits = c(-0.5, 0.5)) +
  theme_void() +
  theme(plot.margin = margin(t = 0, r = 20, b = 0, l = 40),
        panel.background = element_rect(fill = "white", color = NA)) +
  labs(x = "Chromosome")

combined_overall <- plot_grid(
  p_del,
  p_chr + theme(plot.margin = margin(t = -5, b = -5)),
  p_dup,
  ncol = 1,
  rel_heights = c(1, 0.04, 1),
  align = "v"
) + theme(plot.background = element_rect(fill = "white", color = NA))
# x-axis title (the chromosome strip above uses theme_void(), which hides its title)
combined_overall <- ggdraw(add_sub(combined_overall, "Chromosomal position (hg19)", size = 13)) +
  theme(plot.background = element_rect(fill = "white", color = NA))

out_prefix_overall <- file.path(fig_main_dir, paste0("Fig_cnv_density_mirrorplot_", bin_lab, "_no_grid_overall"))
ggsave(paste0(out_prefix_overall, ".pdf"), plot = combined_overall, width = 12, height = 6, dpi = 300)
ggsave(paste0(out_prefix_overall, ".png"), plot = combined_overall, width = 12, height = 6, dpi = 300)

msg("Saved overall mirror-density style plot:")
msg("  %s.(pdf|png)", out_prefix_overall)

# (B) Supplementary: by-phase landscape (del +, dup -)
mirror_phase <- counts_all[phase %in% c("Phase12","Phase3")]
mirror_phase[, freq_signed := ifelse(cnv_type == "dup", -freq, freq)]

mirror_lim <- max(abs(mirror_phase$freq_signed), na.rm = TRUE)
if (!is.finite(mirror_lim) || mirror_lim <= 0) mirror_lim <- 1e-6
mirror_pad <- mirror_lim * 0.05

theme_pub <- theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background  = element_rect(fill = "white", color = NA),
    axis.line   = element_line(color = "black", linewidth = 0.4),
    axis.ticks  = element_line(color = "black", linewidth = 0.4),
    axis.text   = element_text(color = "black"),
    axis.title  = element_text(color = "black"),
    strip.background = element_blank(),
    strip.text  = element_text(face = "bold"),
    plot.title  = element_text(face = "bold"),
    legend.position = "none"
  )

p_mirror_phase <- ggplot(mirror_phase, aes(x = genome_mid, y = freq_signed, color = cnv_type)) +
  geom_segment(aes(xend = genome_mid, y = 0, yend = freq_signed),
               linewidth = 0.2, alpha = 0.9) +
  geom_vline(data = chr_boundaries, aes(xintercept = boundary),
             linewidth = 0.2, alpha = 0.4, color = "grey50", inherit.aes = FALSE) +
  scale_x_continuous(breaks = chr_ticks$tick, labels = chr_ticks$chr) +
  scale_color_manual(values = c(del = "#2E4057", dup = "#B8621B")) +
  coord_cartesian(ylim = c(-(mirror_lim + mirror_pad), mirror_lim + mirror_pad), clip = "off") +
  scale_y_continuous(labels = scales::label_number(accuracy = 0.01)) +
  facet_wrap(~phase, ncol = 1) +
  theme_pub +
  labs(
    x = sprintf("Chromosome (hg19 autosomes; %s bins)", bin_lab),
    y = "Carrier frequency (del +, dup -)",
    title = sprintf("CNV landscape by phase (carrier-based, %s bins)", bin_lab)
  )

save_plot_both(
  p_mirror_phase,
  file.path(fig_supp_dir, paste0("SFig_landscape_mirror_by_phase_", bin_lab)),
  width = 13, height = 7
)

msg("Saved by-phase landscape plot:")
msg("  %s.(pdf|png)", file.path(fig_supp_dir, paste0("SFig_landscape_mirror_by_phase_", bin_lab)))

# -----------------------------
# Histogram + ECDF panels of bin carrier frequencies
# -----------------------------
plot_bins_all     <- counts_all[freq > 0]
plot_bins_overall <- plot_bins_all[phase == "Overall"]
plot_bins_overall[, log10_freq := log10(freq)]

p_hist_overall <- ggplot(plot_bins_overall, aes(x = log10_freq)) +
  geom_histogram(bins = 80) +
  facet_wrap(~cnv_type, ncol = 1, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(x = bquote(log[10]*"("*.(bin_lab)*" bin carrier frequency)"),
       y = sprintf("Number of %s bins", bin_lab),
       title = sprintf("Distribution of non-zero %s bin carrier frequencies (Overall)", bin_lab))

save_plot_both(p_hist_overall,
               file.path(fig_supp_dir, "SFig_landscape_binfreq_hist_log10_overall"),
               width = 9, height = 6)

p_ecdf_overall <- ggplot(plot_bins_overall, aes(x = freq)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  scale_x_log10() +
  facet_wrap(~cnv_type, ncol = 1) +
  theme_classic(base_size = 12) +
  labs(x = sprintf("%s bin carrier frequency (log10 scale)", bin_lab),
       y = "ECDF",
       title = sprintf("ECDF of non-zero %s bin carrier frequencies (Overall)", bin_lab))

save_plot_both(p_ecdf_overall,
               file.path(fig_supp_dir, "SFig_landscape_binfreq_ecdf_log10x_overall"),
               width = 9, height = 6)

plot_bins_phase <- plot_bins_all[phase %in% c("Phase12","Phase3")]
plot_bins_phase[, log10_freq := log10(freq)]

p_hist_byphase <- ggplot(plot_bins_phase, aes(x = log10_freq)) +
  geom_histogram(bins = 80) +
  facet_grid(phase ~ cnv_type, scales = "free_y") +
  theme_classic(base_size = 12) +
  labs(x = bquote(log[10]*"("*.(bin_lab)*" bin carrier frequency)"),
       y = sprintf("Number of %s bins", bin_lab),
       title = sprintf("Distribution of non-zero %s bin carrier frequencies (by phase)", bin_lab))

save_plot_both(p_hist_byphase,
               file.path(fig_supp_dir, "SFig_landscape_binfreq_hist_log10_byphase"),
               width = 10, height = 6)

p_ecdf_byphase <- ggplot(plot_bins_phase, aes(x = freq)) +
  stat_ecdf(geom = "step", linewidth = 0.7) +
  scale_x_log10() +
  facet_grid(phase ~ cnv_type) +
  theme_classic(base_size = 12) +
  labs(x = sprintf("%s bin carrier frequency (log10 scale)", bin_lab),
       y = "ECDF",
       title = sprintf("ECDF of non-zero %s bin carrier frequencies (by phase)", bin_lab))

save_plot_both(p_ecdf_byphase,
               file.path(fig_supp_dir, "SFig_landscape_binfreq_ecdf_log10x_byphase"),
               width = 10, height = 6)

msg("Saved hist + ECDF panels (overall + by phase) to: %s", fig_supp_dir)

sink(type = "output")
close(zz)
