#!/usr/bin/env Rscript
# Stage 2: deletion, duplication and mirror probe matrices from the frozen call set.
#
# Each call is assigned to every probe of its array (CKB1 = phase12, CKB2 = phase3) that lies
# within the call (chromosome and position overlap). The two arrays' matrices are restricted to
# the probes present on both arrays and stacked. In the deletion matrix a probe covered by a
# duplication in the same sample is set to NA, and vice versa; the mirror matrix codes deletions
# as -1 and duplications as +1. Samples with no call are included as all-zero rows.
#
# Input:  results/01_callset/<version>/final/ (callset.rds, samples/zero_cnv_ckids.txt);
#         paths$array1_key, paths$array2_key, paths$probe_map_phase12, paths$probe_map_phase3
# Output: results/02_matrices/<version>/combined_{deletion,duplication,mirror}_matrix_shared_probes_fixed.rds,
#         combined_probe_map_shared_probes.rds, combined_sample_phase_info_fixed.rds

library(data.table)
library(Matrix)

source("R/config.R")
cfg <- load_config(); pa <- cfg$paths
callset_final_dir <- file.path(pa$results_root, "01_callset", cfg$pipeline$callset_version, "final")
output_dir <- file.path(pa$results_root, "02_matrices", cfg$pipeline$callset_version)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
if (length(list.files(output_dir, pattern = "\\.rds$")) && Sys.getenv("FORCE") != "1")
  stop("Matrices already exist in ", output_dir, "; set FORCE=1 only if you mean to rebuild them")

cat("=== CNV matrix creation ===\n")
cat("Callset:", callset_final_dir, "\n")
cat("Output :", output_dir, "\n\n")

# ============================================================================
# STEP 1: Load probe maps (full, per phase) and build join-ready copies
# ============================================================================
cat("Step 1: Loading probe maps...\n")

phase12_map <- fread(pa$probe_map_phase12,
                     col.names = c("chr", "probe_id", "genetic_dist", "pos"))
phase3_map <- fread(pa$probe_map_phase3,
                    col.names = c("chr", "probe_id", "genetic_dist", "pos"))

cat("Phase12 probes (full):", nrow(phase12_map), "\n")
cat("Phase3 probes (full):", nrow(phase3_map), "\n")

# probe_idx is the row in the original probe map (used for the sparse matrix columns).
# chr is character in the join copies only, to match the call coordinates.
phase12_map_idx <- copy(phase12_map); phase12_map_idx[, probe_idx := .I]; phase12_map_idx[, chr := as.character(chr)]
phase3_map_idx  <- copy(phase3_map);  phase3_map_idx[,  probe_idx := .I]; phase3_map_idx[,  chr := as.character(chr)]
setkey(phase12_map_idx, chr, pos)
setkey(phase3_map_idx,  chr, pos)

# Probes on both arrays, used to subset the finished phase matrices
shared_probes <- intersect(phase12_map$probe_id, phase3_map$probe_id)
shared_probe_map <- phase12_map[probe_id %in% shared_probes]
setorder(shared_probe_map, chr, pos)
shared_order <- shared_probe_map$probe_id

cat("Shared probes (intersection):", length(shared_probes), "\n\n")

# ============================================================================
# Samples with no call (one list, so the deletion and duplication matrices have identical rows)
# ============================================================================
cat("Loading sample info...\n")

array1_samples <- fread(pa$array1_key)$IID
array2_samples <- fread(pa$array2_key)$IID
zero_cnv_samples <- readLines(file.path(callset_final_dir, "samples", "zero_cnv_ckids.txt"))

phase12_zero_samples <- zero_cnv_samples[zero_cnv_samples %in% array1_samples]
phase3_zero_samples  <- zero_cnv_samples[zero_cnv_samples %in% array2_samples]

cat("Zero-CNV samples (phase12):", length(phase12_zero_samples), "\n")
cat("Zero-CNV samples (phase3):", length(phase3_zero_samples), "\n\n")

# ============================================================================
# STEP 2: Load raw CNV calls
# ============================================================================
cat("Step 2: Loading raw CNV calls...\n")

cs <- readRDS(file.path(callset_final_dir, "callset.rds"))
# column names used in the rest of this script
all_cnvs <- cs[, .(region = Region, n_probes = NumSNP, length = Length, source_file = File,
                   start_probe = StartSNP, end_probe = EndSNP, score = Conf, state = State,
                   cn = CN, chr = Chromosome, start_pos = StartPos, end_pos = EndPos, type = CN_type,
                   conf = n_merged, sample_id = CKID, genes = Genes, overlap = Gene_distance)]
rm(cs)

cat("CHECKPOINT - Total raw CNVs loaded:", nrow(all_cnvs), "\n")
cat("  Deletions:", sum(all_cnvs$type == "del"), "| Duplications:", sum(all_cnvs$type == "dup"), "\n\n")

# ============================================================================
# STEP 3: Assign each call to its array (phase)
# ============================================================================
cat("Step 3: Assigning phases...\n")

all_cnvs[, phase := ifelse(sample_id %in% array1_samples, "phase12",
                           ifelse(sample_id %in% array2_samples, "phase3", NA))]
n_unassigned <- sum(is.na(all_cnvs$phase))
if(n_unassigned > 0) cat("  WARNING:", n_unassigned, "CNVs dropped (sample not in array1/array2)\n")
all_cnvs <- all_cnvs[!is.na(phase)]
all_cnvs[, chr_clean := gsub("^chr", "", chr)]

phase12_cnv_samples <- unique(all_cnvs[phase == "phase12"]$sample_id)
phase3_cnv_samples  <- unique(all_cnvs[phase == "phase3"]$sample_id)

cat("CHECKPOINT - CNVs after phase assignment:", nrow(all_cnvs), "\n")
cat("  Phase12:", sum(all_cnvs$phase=="phase12"), "CNVs /", length(phase12_cnv_samples), "samples\n")
cat("  Phase3: ", sum(all_cnvs$phase=="phase3"),  "CNVs /", length(phase3_cnv_samples), "samples\n\n")

# ============================================================================
# STEP 4: Complete sample sets (CNV samples + zero-CNV samples)
# ============================================================================
complete_phase12_samples <- sort(unique(c(phase12_cnv_samples, phase12_zero_samples)))
complete_phase3_samples  <- sort(unique(c(phase3_cnv_samples,  phase3_zero_samples)))

cat("Step 4: Complete sample sets - Phase12:", length(complete_phase12_samples),
    "| Phase3:", length(complete_phase3_samples),
    "| Total:", length(complete_phase12_samples) + length(complete_phase3_samples), "\n\n")

# ============================================================================
# STEP 5/6: Build the phase-specific matrices (chromosome + position overlap join)
# ============================================================================
cat("Step 5/6: Building phase-specific matrices (vectorized overlap join)...\n")

build_matrix <- function(cnv_subset, sample_set, probe_map_idx, probe_map_full, label) {
  cnv_subset <- copy(cnv_subset)
  cnv_subset[, sample_idx := match(sample_id, sample_set)]
  n_dropped <- sum(is.na(cnv_subset$sample_idx))
  if(n_dropped > 0) cat("    WARNING (", label, "):", n_dropped, "CNVs - sample not in complete set (dropped)\n")
  cnv_valid <- cnv_subset[!is.na(sample_idx)]

  if(nrow(cnv_valid) == 0) {
    return(data.table(sample_idx = integer(), probe_idx = integer(), value = integer()))
  }

  matches <- probe_map_idx[
    cnv_valid,
    on = .(chr = chr_clean, pos >= start_pos, pos <= end_pos),
    nomatch = 0L,
    allow.cartesian = TRUE
  ]

  triplets <- unique(data.table(sample_idx = matches$sample_idx, probe_idx = matches$probe_idx, value = 1L))
  cat("    (", label, ") CNVs:", nrow(cnv_valid), "-> triplets:", nrow(triplets), "\n")
  return(triplets)
}

del12_t <- build_matrix(all_cnvs[phase=="phase12" & type=="del"], complete_phase12_samples, phase12_map_idx, phase12_map, "phase12/del")
dup12_t <- build_matrix(all_cnvs[phase=="phase12" & type=="dup"], complete_phase12_samples, phase12_map_idx, phase12_map, "phase12/dup")
del3_t  <- build_matrix(all_cnvs[phase=="phase3"  & type=="del"], complete_phase3_samples,  phase3_map_idx,  phase3_map,  "phase3/del")
dup3_t  <- build_matrix(all_cnvs[phase=="phase3"  & type=="dup"], complete_phase3_samples,  phase3_map_idx,  phase3_map,  "phase3/dup")
cat("\n")

del_phase12 <- sparseMatrix(i=del12_t$sample_idx, j=del12_t$probe_idx, x=del12_t$value,
                             dims=c(length(complete_phase12_samples), nrow(phase12_map)),
                             dimnames=list(complete_phase12_samples, phase12_map$probe_id))
dup_phase12 <- sparseMatrix(i=dup12_t$sample_idx, j=dup12_t$probe_idx, x=dup12_t$value,
                             dims=c(length(complete_phase12_samples), nrow(phase12_map)),
                             dimnames=list(complete_phase12_samples, phase12_map$probe_id))
del_phase3  <- sparseMatrix(i=del3_t$sample_idx,  j=del3_t$probe_idx,  x=del3_t$value,
                             dims=c(length(complete_phase3_samples), nrow(phase3_map)),
                             dimnames=list(complete_phase3_samples, phase3_map$probe_id))
dup_phase3  <- sparseMatrix(i=dup3_t$sample_idx,  j=dup3_t$probe_idx,  x=dup3_t$value,
                             dims=c(length(complete_phase3_samples), nrow(phase3_map)),
                             dimnames=list(complete_phase3_samples, phase3_map$probe_id))

cat("CHECKPOINT - Full phase-specific matrix nnz (before shared-probe subset):\n")
cat("  del_phase12:", nnzero(del_phase12), "| dup_phase12:", nnzero(dup_phase12), "\n")
cat("  del_phase3: ", nnzero(del_phase3),  "| dup_phase3: ", nnzero(dup_phase3), "\n\n")

# Free the triplet tables and join scratch data now that matrices are built
rm(del12_t, dup12_t, del3_t, dup3_t)
gc()

# ============================================================================
# STEP 7: Subset to shared probes, stack phases
# ============================================================================
cat("Step 7: Subsetting to shared probes and stacking...\n")

del_phase12_s <- del_phase12[, shared_order]
del_phase3_s  <- del_phase3[, shared_order]
combined_del  <- rbind(del_phase12_s, del_phase3_s)
rm(del_phase12, del_phase3, del_phase12_s, del_phase3_s)
gc()

dup_phase12_s <- dup_phase12[, shared_order]
dup_phase3_s  <- dup_phase3[, shared_order]
combined_dup  <- rbind(dup_phase12_s, dup_phase3_s)
rm(dup_phase12, dup_phase3, dup_phase12_s, dup_phase3_s)
gc()

cat("CHECKPOINT - Combined nnz (before masking):\n")
cat("  combined_del:", nnzero(combined_del), "| combined_dup:", nnzero(combined_dup), "\n")
cat("  Dimensions:", dim(combined_del), "\n\n")

# ============================================================================
# STEP 8: Masking (both position sets taken from the unmasked matrices)
# ============================================================================
cat("Step 8: Applying masking...\n")

dup_positions <- which(combined_dup == 1, arr.ind = TRUE, useNames = FALSE)
del_positions <- which(combined_del == 1, arr.ind = TRUE, useNames = FALSE)

if(nrow(dup_positions) > 0) combined_del[dup_positions] <- NA_real_
if(nrow(del_positions) > 0) combined_dup[del_positions] <- NA_real_

cat("CHECKPOINT - After masking - del:", nnzero(combined_del), "| dup:", nnzero(combined_dup), "\n\n")

# ============================================================================
# STEP 9: Mirror matrix
# ============================================================================
cat("Step 9: Creating mirror matrix...\n")
combined_mirror <- combined_dup
if(nrow(del_positions) > 0) combined_mirror[del_positions] <- -1
cat("CHECKPOINT - Mirror nnz:", nnzero(combined_mirror), "\n\n")

# ============================================================================
# STEP 10/11: Save
# ============================================================================
sample_phase_info <- data.table(
  sample_id = c(complete_phase12_samples, complete_phase3_samples),
  phase = c(rep("phase12", length(complete_phase12_samples)), rep("phase3", length(complete_phase3_samples)))
)

saveRDS(combined_del,    file.path(output_dir, "combined_deletion_matrix_shared_probes_fixed.rds"))
saveRDS(combined_dup,    file.path(output_dir, "combined_duplication_matrix_shared_probes_fixed.rds"))
saveRDS(combined_mirror, file.path(output_dir, "combined_mirror_matrix_shared_probes_fixed.rds"))
saveRDS(shared_probe_map, file.path(output_dir, "combined_probe_map_shared_probes.rds"))
saveRDS(sample_phase_info, file.path(output_dir, "combined_sample_phase_info_fixed.rds"))

cat("=== FINAL SUMMARY ===\n")
cat("Total samples:", nrow(combined_del), "| Shared probes:", ncol(combined_del), "\n")
cat("Deletions:", nnzero(combined_del), "| Duplications:", nnzero(combined_dup), "| Mirror:", nnzero(combined_mirror), "\n")
cat("=== DONE ===\n")
