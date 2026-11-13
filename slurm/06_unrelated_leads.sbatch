#!/bin/bash
#SBATCH --job-name=unrel_leads
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%A_%a.log
# Stage 6: unrelated sensitivity at the lead probes of the main signals. Same GWAS code and model as the main
# analysis (scripts/04_gwas/01_gwas_task.R), unrelated participants, lead probes only, no minimum-carrier filter
# (so leads with < 15 unrelated carriers still get an estimate).
# Submit with scripts/06_followup/run_unrelated_leads.sh (it sets --array to the number of tasks).
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
read -r MODEL PHENO LIST < <(sed -n "${SLURM_ARRAY_TASK_ID}p" results/06_followup/$VER/unrelated_leads/tasks.txt)
echo "Task $SLURM_ARRAY_TASK_ID: $MODEL $PHENO | $(date)"
PROBE_LIST="$LIST" OUT_TAG=leads MIN_CARRIERS=1 Rscript scripts/04_gwas/01_gwas_task.R "$MODEL" "$PHENO" unrelated
echo "Finished $(date)"
