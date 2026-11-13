#!/bin/bash
# Makes the lead lists, submits one unrelated GWAS task per model x trait, then the Table S9 / Figure S7 job.
# Run from the repository root:  bash scripts/06_followup/run_unrelated_leads.sh
set -euo pipefail
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml); mkdir -p logs
mkdir -p results/06_followup/$VER/unrelated_leads/logs
Rscript scripts/06_followup/01_unrelated_lead_lists.R
N=$(wc -l < results/06_followup/$VER/unrelated_leads/tasks.txt)
J1=$(sbatch --parsable --array=1-"$N" slurm/06_unrelated_leads.sbatch)
# the --wrap job needs R on its PATH (as in the slurm/ scripts)
J2=$(sbatch --parsable --dependency=afterany:"$J1" --job-name=unrel_tables --mem=16G --time=01:00:00 --partition=general \
     --output=logs/unrel_tables_%j.log \
     --wrap "cd $PWD; Rscript scripts/06_followup/02_table_s9_unrelated.R && Rscript scripts/06_followup/03_figure_s7_unrelated.R")
echo "Submitted $N unrelated lead tasks (array $J1) and the Table S9 / Figure S7 job ($J2, starts when they finish)."
echo "Then:  cat logs/unrel_tables_$J2.log"
