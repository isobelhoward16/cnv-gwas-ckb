#!/bin/bash
#SBATCH --job-name=batch_affyhmm
#SBATCH --output=logs/%x_%j.log
#SBATCH --time=1:00:00
#SBATCH -p general
# Submits extract_and_submit_cnv_calling.sh for each of the 20 region x array strata, at most 6 at a time.
# Usage (repository root):  mkdir -p logs; sbatch calling/batch_affyhmm_launcher.sh

# Define region-phase combinations
REGION_PHASE_LIST=(
    "Gansu_Rural Phase12"
    "Gansu_Rural Phase3"
    "Harbin_Urban Phase12"
    "Harbin_Urban Phase3"
    "Haikou_Urban Phase12"
    "Haikou_Urban Phase3"
    "Henan_Rural Phase12"
    "Henan_Rural Phase3"
    "Hunan_Rural Phase12"
    "Hunan_Rural Phase3"
    "Liuzhou_Urban Phase12"
    "Liuzhou_Urban Phase3"
    "Qingdao_Urban Phase12"
    "Qingdao_Urban Phase3"
    "Sichuan_Rural Phase12"
    "Sichuan_Rural Phase3"
    "Suzhou_Urban Phase12"
    "Suzhou_Urban Phase3"
    "Zhejiang_Rural Phase12"
    "Zhejiang_Rural Phase3"
)

EXTRACT_SCRIPT="calling/extract_and_submit_cnv_calling.sh"
MAX_PARALLEL=6

# Submit the first batch
JOB_IDS=()
for i in $(seq 0 $((MAX_PARALLEL - 1))); do
    read REGION PHASE <<< "${REGION_PHASE_LIST[i]}"
    JOB_ID=$(sbatch --parsable "$EXTRACT_SCRIPT" "$REGION" "$PHASE")
    echo "Submitted job for $REGION $PHASE with Job ID: $JOB_ID"
    JOB_IDS+=("$JOB_ID")
done

# Submit remaining jobs with dependencies
for (( i=MAX_PARALLEL; i<${#REGION_PHASE_LIST[@]}; i++ )); do
    read REGION PHASE <<< "${REGION_PHASE_LIST[i]}"
    DEP_LIST=$(IFS=:; echo "${JOB_IDS[*]}")
    JOB_ID=$(sbatch --parsable --dependency=afterany:${DEP_LIST} "$EXTRACT_SCRIPT" "$REGION" "$PHASE")
    echo "Submitted job for $REGION $PHASE with Job ID: $JOB_ID (after ${JOB_IDS[*]})"
    JOB_IDS=("${JOB_IDS[@]:1}" "$JOB_ID")
done

echo "All jobs submitted with up to $MAX_PARALLEL concurrent runs."
