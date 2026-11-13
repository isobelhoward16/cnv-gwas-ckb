#!/usr/bin/env Rscript
# Stage 8 set-up: writes the small config file read by the descriptive scripts (stage 8), pointing at the
# stage 1 callset and its list of samples with no CNV call, with all outputs under results/08_descriptive/<version>/.
# Output: results/08_descriptive/<version>/descriptive_paths.yaml, qc_thresholds.yaml,
#         inputs/qc_passed_zero_cnv_samples.txt
source("R/config.R")
cfg <- load_config(); pa <- cfg$paths; ver <- cfg$pipeline$callset_version; q <- cfg$pipeline$sample_qc
root <- file.path(pa$results_root, "08_descriptive", ver)
fin  <- file.path(pa$results_root, "01_callset", ver, "final")
dir.create(file.path(root, "inputs"), recursive = TRUE, showWarnings = FALSE)
for (p in c(pa$qc_gwas_file, pa$phenotype_csv, pa$refgene_annotation, pa$cytoband_landscape))
  if (is.null(p) || !file.exists(p)) stop("Missing input in config/paths.yaml: ", p)

zero <- readLines(file.path(fin, "samples", "zero_cnv_ckids.txt"))
zf <- file.path(root, "inputs", "qc_passed_zero_cnv_samples.txt")
fwrite(data.table(GID = zero), zf)            # one column, GID

y <- list(run_name = paste0("callset_", ver),
          gene_cleaned_merged_dir = file.path(fin, "calls"),
          qc_master_file = pa$qc_gwas_file,
          analysis_dataframe_csv = pa$phenotype_csv,
          qc_passed_zero_cnv_file = zf,
          out_root = root,
          cytoband_file = pa$cytoband_landscape,
          refgene_annotation = pa$refgene_annotation,
          annotated_cnv_dir = file.path(root, "inputs", "refgene_annotated"),
          bin_size = cfg$pipeline$descriptive$density_bin_size)
write_yaml(y, file.path(root, "descriptive_paths.yaml"))
write_yaml(list(LRR_SD_max = q$max_lrr_sd, WF_abs_max = q$max_wf, BAF_drift_max = q$max_baf_drift),
           file.path(root, "qc_thresholds.yaml"))
cat("Zero-CNV samples:", length(zero), "\nWritten:", file.path(root, "descriptive_paths.yaml"), "\n")
cat(readLines(file.path(root, "descriptive_paths.yaml")), sep = "\n")
