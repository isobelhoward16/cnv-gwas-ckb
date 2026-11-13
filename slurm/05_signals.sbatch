#!/bin/bash
#SBATCH --job-name=signals
#SBATCH --ntasks=1
#SBATCH --mem=100G
#SBATCH --time=3-00:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 5: significant probes, CNVRs and signal relationships, for all participants and for unrelated participants.
# Run after the Meff summary (scripts/03_meff/02_meff_summary.R) and GWAS check have passed.
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
for REL in ${RELS:-all unrelated}; do
  echo "==================== $REL | $(date)"
  Rscript scripts/05_signals/01_sig_extraction.R "$REL"
  Rscript scripts/05_signals/02_cnvr_annotation.R "$REL"
  Rscript scripts/05_signals/03_signal_relationships.R "$REL"
  Rscript scripts/05_signals/05_check_signals.R "$REL" || echo "05_check_signals.R failed for $REL: look at the log"
done
echo "Finished $(date)"
