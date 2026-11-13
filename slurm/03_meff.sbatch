#!/bin/bash
#SBATCH --job-name=meff
#SBATCH --ntasks=1
#SBATCH --mem=80G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --array=1-78%10
#SBATCH --output=logs/%x_%A_%a.log
# Stage 3: one Meff task per model x phenotype x sample set (3 x 13 x 2 = 78).
# Submit from the repository root:  sbatch slurm/03_meff.sbatch
# Afterwards:                       Rscript scripts/03_meff/02_meff_summary.R
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
read -r MODEL PHENO REL < <(Rscript R/tasks.R 2>/dev/null | sed -n "${SLURM_ARRAY_TASK_ID}p")
echo "Task $SLURM_ARRAY_TASK_ID: $MODEL $PHENO $REL | $(date)"
Rscript scripts/03_meff/01_meff_task.R "$MODEL" "$PHENO" "$REL"
echo "Finished $(date)"
