#!/usr/bin/env Rscript
# Stage 8, step 4: exonic and genic (gene body +/-10 kb) annotation of each call against RefSeq genes (hg19),
# used for Table S4.

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(GenomicRanges)
  library(parallel)
  library(HGNChelper)
  library(yaml)
})

# ==== Config ====
paths_yaml <- Sys.getenv("DESCRIPTIVE_PATHS")   # written by scripts/08_descriptive/00_setup_inputs.R
if (!nzchar(paths_yaml) || !file.exists(paths_yaml)) stop("Set DESCRIPTIVE_PATHS to results/08_descriptive/<version>/descriptive_paths.yaml (see slurm/08_descriptive.sbatch)")
paths <- yaml::read_yaml(paths_yaml)
cnv_dir <- paths$gene_cleaned_merged_dir
if (!dir.exists(cnv_dir)) stop("CNV directory not found: ", cnv_dir)

refgene_file <- paths$refgene_annotation
if (!file.exists(refgene_file)) stop("RefGene file not found: ", refgene_file)

# Derive run_id and write to a run-specific output dir to avoid overwrites
run_id <- paths$run_name
output_dir <- paths$annotated_cnv_dir

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

expand_kb <- 10  # gene body expansion (kb)
n_cores <- max(1, min(detectCores() - 1, 8))

logfile <- file.path(
  output_dir,
  paste0("exonic_annotation_", run_id, "_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log")
)
sink(logfile, split = TRUE)

cat("======================================\n")
cat("CNV Gene Annotation (Exonic + Genic +/-", expand_kb, "kb)\n")
cat("Run:", run_id, "\n")
cat("CNV dir:", cnv_dir, "\n")
cat("RefGene:", refgene_file, "\n")
cat("Output dir:", output_dir, "\n")
cat("Using", n_cores, "cores\n")
cat("======================================\n\n")

# ==== Load and normalise RefGene ====
cat("Loading RefGene...\n")
refgene <- fread(refgene_file, header = FALSE, sep = "\t")
setnames(refgene, c(
  "bin", "name", "chrom", "strand", "txStart", "txEnd",
  "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds",
  "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames"
))
cat("Total genes in refGene:", nrow(refgene), "\n")

cat("Normalising gene symbols with HGNChelper...\n")
hgnc_map <- checkGeneSymbols(refgene$name2, unmapped.as.na = FALSE)
refgene$name2 <- ifelse(!is.na(hgnc_map$Suggested.Symbol), hgnc_map$Suggested.Symbol, refgene$name2)

# ==== Exon ranges ====
cat("Extracting exon ranges...\n")
extract_exons <- function(row) {
  starts <- as.numeric(strsplit(row$exonStarts, ",")[[1]])
  ends   <- as.numeric(strsplit(row$exonEnds, ",")[[1]])
  if (length(starts) != length(ends) || length(starts) == 0) return(NULL)
  data.table(
    chrom = row$chrom,
    start = starts,
    end   = ends,
    gene  = row$name2
  )
}
exon_list <- mclapply(1:nrow(refgene), function(i) extract_exons(refgene[i]), mc.cores = n_cores)
exons <- rbindlist(exon_list[!sapply(exon_list, is.null)])
exons_gr <- makeGRangesFromDataFrame(exons, keep.extra.columns = TRUE)
cat("Exons extracted:", length(exons_gr), "\n")

# ==== Full gene body +/- kb ====
cat("Building gene body ranges +/-", expand_kb, "kb...\n")
genes_gr <- makeGRangesFromDataFrame(
  data.table(
    chrom = refgene$chrom,
    start = pmax(1, refgene$txStart - expand_kb * 1000),
    end   = refgene$txEnd + expand_kb * 1000,
    gene  = refgene$name2
  ),
  keep.extra.columns = TRUE
)
cat("Gene body ranges created:", length(genes_gr), "\n")

# ==== Process CNVs ====
cnv_files <- list.files(cnv_dir, pattern = "\\.txt$", full.names = TRUE)
if (length(cnv_files) == 0) stop("No CNV files in ", cnv_dir)

cat("Found", length(cnv_files), "CNV files\n")

process_cnv <- function(cnv_file) {
  cat("Processing:", basename(cnv_file), "\n")
  cnv <- tryCatch(fread(cnv_file), error = function(e) { cat("  Failed to read:", e$message, "\n"); return(NULL) })
  if (is.null(cnv) || nrow(cnv) == 0) return(NULL)

  # Harmonise columns if present (first 16 expected)
  expected_cols <- c("coord","numsnp","length","source_path","startsnp","endsnp",
                     "conf","state","cn","chr","start","end",
                     "cnv_type","n_merged","ckid","genes")
  if (ncol(cnv) >= length(expected_cols)) setnames(cnv, 1:length(expected_cols), expected_cols)

  # Build GRanges for CNVs
  cnv_gr <- makeGRangesFromDataFrame(cnv, keep.extra.columns = TRUE)

  # Overlaps: exons
  ol_exon <- findOverlaps(cnv_gr, exons_gr)
  cnv$exonic_genes <- ""
  if (length(ol_exon) > 0) {
    exon_map <- data.table(
      idx  = queryHits(ol_exon),
      gene = mcols(exons_gr)$gene[subjectHits(ol_exon)]
    )
    exon_map <- exon_map[, .(exonic_genes = paste(unique(gene), collapse = ",")), by = idx]
    cnv$exonic_genes[exon_map$idx] <- exon_map$exonic_genes
  }

  # Overlaps: gene body +/- kb
  ol_gene <- findOverlaps(cnv_gr, genes_gr)
  cnv$genic_genes <- ""
  if (length(ol_gene) > 0) {
    gene_map <- data.table(
      idx  = queryHits(ol_gene),
      gene = mcols(genes_gr)$gene[subjectHits(ol_gene)]
    )
    gene_map <- gene_map[, .(genic_genes = paste(unique(gene), collapse = ",")), by = idx]
    cnv$genic_genes[gene_map$idx] <- gene_map$genic_genes
  }

  # Write output (prefix with run_id to avoid any collisions)
  out_file <- file.path(output_dir, paste0(run_id, "_annotated_", basename(cnv_file)))
  fwrite(cnv, out_file, sep = "\t")

  cat("  CNVs:", nrow(cnv),
      "| Exonic:", sum(cnv$exonic_genes != ""),
      "| Genic:",  sum(cnv$genic_genes  != ""), "\n")

  return(list(file = basename(cnv_file),
              total = nrow(cnv),
              exonic = sum(cnv$exonic_genes != ""),
              genic  = sum(cnv$genic_genes  != "")))
}

cat("Annotating CNVs...\n")
results <- mclapply(cnv_files, process_cnv, mc.cores = n_cores)
results <- results[!sapply(results, is.null)]

# ==== Summary ====
cat("\n=== Summary ===\n")
if (length(results) > 0) {
  total      <- sum(sapply(results, `[[`, "total"))
  total_exon <- sum(sapply(results, `[[`, "exonic"))
  total_gene <- sum(sapply(results, `[[`, "genic"))
  cat("Total CNVs:", total, "\n")
  cat("With exonic overlap:", total_exon, "(", ifelse(total>0, round(100 * total_exon / total, 2), 0), "% )\n")
  cat("With genic  overlap:", total_gene, "(", ifelse(total>0, round(100 * total_gene / total, 2), 0), "% )\n")
} else {
  cat("No CNVs processed.\n")
}
sink()
