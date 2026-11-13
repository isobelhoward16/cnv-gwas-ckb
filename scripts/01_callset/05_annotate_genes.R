#!/usr/bin/env Rscript
# Stage 1, step 5: genes overlapping each call. PennCNV scan_region.pl is run with -refgene and
# -reflink on the call coordinates; gene identifiers are then converted to symbols with the reflink
# file (column 1 -> column 2). Calls with no gene get NOT_FOUND. The conversion is skipped if the
# reflink file is not a two-column identifier -> symbol table.
#
# Input:  04_sample_filtered/calls_sample_filtered.rds; paths$refgene, paths$reflink, paths$penncnv_dir
# Output: 05_annotated/calls_annotated.rds (adds Genes and Gene_distance)
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths
STEP <- "05_annotate_genes"
start_step(cfg, STEP)
out_dir <- callset_dir(cfg, "05_annotated")

x <- readRDS(callset_dir(cfg, "04_sample_filtered", "calls_sample_filtered.rds"))
record(cfg, STEP, "calls in", nrow(x))

# query file: region and row number only, so output rows join back exactly
qf <- file.path(out_dir, "scan_region_query.tmp"); of <- file.path(out_dir, "scan_region_out.tmp")
ef <- file.path(out_dir, "scan_region_stderr.log")
fwrite(data.table(paste0(x$Chromosome, ":", x$StartPos, "-", x$EndPos), seq_len(nrow(x))), qf, sep = "\t", col.names = FALSE)
st <- system2("perl", c(file.path(pa$penncnv_dir, "scan_region.pl"), qf, pa$refgene, "-refgene", "-reflink", pa$reflink),
              stdout = of, stderr = ef)
if (!identical(st, 0L)) stop("scan_region.pl failed (exit ", st, "); see ", ef)
o <- fread(of, sep = "\t", header = FALSE, colClasses = "character", quote = "", fill = TRUE)
if (ncol(o) < 4) stop("Unexpected scan_region.pl output: ", ncol(o), " columns")
o <- o[, .(row = as.integer(V2), genes_raw = V3, Gene_distance = V4)]
if (anyDuplicated(o$row) || !setequal(o$row, seq_len(nrow(x)))) stop("scan_region.pl output does not match the input rows one to one")
setorder(o, row)
record(cfg, STEP, "calls annotated by scan_region.pl", nrow(o))

# identifier -> symbol table (column 1 -> column 2; later lines win)
rl <- fread(pa$reflink, sep = "\t", header = FALSE, colClasses = "character", quote = "", fill = TRUE, select = 1:2)
map <- rl[!duplicated(V1, fromLast = TRUE)]; map <- setNames(map$V2, map$V1)
toks <- strsplit(gsub("[[:space:]]+", "", o$genes_raw), ",", fixed = TRUE)
all_tok <- unlist(toks); all_tok <- all_tok[all_tok != "NOT_FOUND"]
p_mapped <- if (length(all_tok)) mean(all_tok %in% names(map)) else 0
mapped_vals <- map[intersect(unique(all_tok), names(map))]
p_space <- if (length(mapped_vals)) mean(grepl("[[:space:]]", mapped_vals)) else 0
use_map <- p_mapped > 0.5 && p_space < 0.01
record(cfg, STEP, "share of scan_region identifiers found in reflink column 1", round(p_mapped, 4))
record(cfg, STEP, "share of mapped values containing spaces", round(p_space, 4))
record(cfg, STEP, "identifier conversion", if (use_map) "applied (as clean_gene_symbols.sh)" else "skipped (scan_region output already symbols)")

clean <- function(t) {
  g <- if (use_map) unname(map[t]) else t
  g <- unique(g[!is.na(g) & g != "" & g != "NOT_FOUND"])
  if (length(g)) paste(sort(g), collapse = ",") else "NOT_FOUND"
}
x[, Genes := vapply(toks, clean, "")]
x[, Gene_distance := o$Gene_distance]
record(cfg, STEP, "calls overlapping at least one gene", sum(x$Genes != "NOT_FOUND"))
cat("\n   Examples (scan_region output -> Genes):\n")
ex <- which(o$genes_raw != "NOT_FOUND")[1:min(5, sum(o$genes_raw != "NOT_FOUND"))]
print(data.table(region = paste0(x$Chromosome, ":", x$StartPos, "-", x$EndPos)[ex], scan_region = o$genes_raw[ex], Genes = x$Genes[ex]))

saveRDS(x, file.path(out_dir, "calls_annotated.rds"))
invisible(file.remove(qf, of))
