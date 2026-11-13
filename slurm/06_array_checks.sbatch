#!/bin/bash
#SBATCH --job-name=array_checks
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 6: array-stratified validation at the lead probes (Figures S3-S5).
# Submitted by scripts/run_post_stage5.sh.
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Callset version: $(awk '/^callset_version:/{print $2}' config/pipeline.yaml) | commit: $(git rev-parse --short HEAD) | $(date)"
Rscript scripts/06_followup/04_array_validation.R
echo "Finished $(date)"
