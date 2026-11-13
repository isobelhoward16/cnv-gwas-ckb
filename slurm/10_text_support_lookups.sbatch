#!/bin/bash
#SBATCH --job-name=lead16a
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --partition=general
#SBATCH --array=1-9
#SBATCH --output=logs/%x_%A_%a.log
# Glucose, SBP and DBP at the 16p11.2 BP2-BP3 lead probe (AX-31870881), deletion, duplication and mirror models
# (reported in the 16p11.2 Results paragraph). Same GWAS code and model as the main analysis, all
# participants, this probe only, no minimum-carrier filter (the deletion has fewer than 15 carriers at this probe).
# Output: results/04_gwas/<version>/all/sensitivity/*_lead16a.csv, read by scripts/10_paper/05_text_support.R.
# Before submitting, write the probe list:  echo AX-31870881 > results/10_paper/<version>/lead16a_probe.txt
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | commit: $(git rev-parse --short HEAD) | $(date)"
MODELS=(deletion duplication mirror); PHENOS=(random_glucose_rint sbp_mean_rint dbp_mean_rint)
i=$((SLURM_ARRAY_TASK_ID - 1)); MODEL=${MODELS[$((i / 3))]}; PHENO=${PHENOS[$((i % 3))]}
L=results/10_paper/$VER/lead16a_probe.txt
echo "Task $SLURM_ARRAY_TASK_ID: $MODEL $PHENO | probe list $L: $(cat "$L")"
PROBE_LIST="$L" OUT_TAG=lead16a MIN_CARRIERS=1 Rscript scripts/04_gwas/01_gwas_task.R "$MODEL" "$PHENO" all
echo "Finished $(date)"
