#!/bin/bash
#SBATCH --job-name=matrices
#SBATCH --ntasks=1
#SBATCH --mem=60G
#SBATCH --time=1-00:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 2: CNV matrices (samples x probes) from the stage 1 callset.
# Submit from the repository root:  sbatch slurm/02_matrices.sbatch
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
Rscript scripts/02_matrices/01_build_matrices.R
echo "Finished $(date)"
