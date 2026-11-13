#!/bin/bash
#SBATCH --job-name=callset
#SBATCH --ntasks=1
#SBATCH --mem=48G
#SBATCH --time=12:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 1: build the callset from the raw PennCNV calls. Submitted by scripts/run_main_chain.sh.
# Submit from the repository root:  sbatch slurm/01_callset.sbatch
# Set up R here (packages listed in the README), e.g. module load R
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | repository: $PWD | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none) | $(date)"
for s in 01_parse_merge_calls 02_sample_qc 03_filter_calls 04_filter_high_count 05_annotate_genes 06_freeze_callset; do
  Rscript "scripts/01_callset/$s.R"
done
echo "Finished $(date)"
