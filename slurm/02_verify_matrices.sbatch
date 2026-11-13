#!/bin/bash
#SBATCH --job-name=verify_matrices
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 2 check: independent rebuild of every matrix entry from the callset (scripts/02_matrices/03_verify_matrices.R).
# Runs alongside Meff and the GWAS (they do not wait for it); read its log before stage 5 results are used.
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
Rscript scripts/02_matrices/03_verify_matrices.R
echo "Finished $(date)"
