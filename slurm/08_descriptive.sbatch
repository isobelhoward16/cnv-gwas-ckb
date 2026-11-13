#!/bin/bash
#SBATCH --job-name=descriptive
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 8: sample manifest, Figure 1, Figure 2, exonic annotation, Table S4, Figure S6.
# Needs only the stage 1 callset. Submit from the repository root: sbatch slurm/08_descriptive.sbatch
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
Rscript scripts/08_descriptive/00_setup_inputs.R
export DESCRIPTIVE_PATHS="$PWD/results/08_descriptive/$VER/descriptive_paths.yaml"
for s in 01_build_manifest 02_figure1_callset 03_figure2_density 04_exon_annotation 05_table_s4_genomic_disorders 06_figure_s6_ckb_vs_ukb; do
  echo "==================== $s | $(date)"
  Rscript "scripts/08_descriptive/$s.R"
done
echo "Finished $(date)"
