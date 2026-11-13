#!/bin/bash
#SBATCH --job-name=external
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 7: UK Biobank comparison (Tables S7, S8), GWAS Catalog counts per CNVR, segmental-duplication context.
# Each step runs even if an earlier one fails; the job fails at the end if any step did. Submitted by scripts/run_post_stage5.sh.
# Set up R here (packages listed in the README), e.g. module load R
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Callset version: $(awk '/^callset_version:/{print $2}' config/pipeline.yaml) | commit: $(git rev-parse --short HEAD) | $(date)"
FAILED=""
for s in 02_table_s7_cnvr_freq_ukb 03_table_s8_ukb_lookup 04_catalog_counts 05_segdup_overlap; do
  echo; echo "==================== $s | $(date)"
  Rscript scripts/07_external/$s.R || FAILED="$FAILED $s"
done
echo; echo "Finished $(date) | failed:${FAILED:- none}"
[ -z "$FAILED" ]
