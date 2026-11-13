#!/bin/bash
#SBATCH --job-name=meff_gwas_checks
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=06:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# After all Meff and GWAS tasks: Meff thresholds (needed by stage 5), the GWAS completeness check, and Table S2
# (its trait N must equal the GWAS N). Stage 5 starts only if this job succeeds.
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
echo "==================== Meff summary"; Rscript scripts/03_meff/02_meff_summary.R
echo "==================== GWAS check";   Rscript scripts/04_gwas/02_check_gwas.R
echo "==================== Table S2";     Rscript scripts/08_descriptive/07_table_s2_population.R || echo "Table S2 failed: look at the log (stage 5 continues)"
echo "Finished $(date)"
