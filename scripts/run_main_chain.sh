#!/bin/bash
# Submits the main analysis for the callset version in config/pipeline.yaml, each job waiting for the one it needs:
#   stage 1 callset -> stage 8 descriptive (Figures 1-2, Table S4, Figure S6)
#                   -> stage 2 matrices -> matrix verification (in parallel with the next step)
#                                       -> stage 3 Meff (78 tasks) + stage 4 GWAS (78 tasks)
#                                          -> Meff thresholds, GWAS check, Table S2
#                                             -> stage 5 signals, all participants and unrelated (in parallel)
# Stages 6, 7, 9 and 10 are run after the stage 5 signals have been checked (scripts/06_followup/run_unrelated_leads.sh,
# slurm/06_fasting_rpg.sbatch, then scripts/run_post_stage5.sh).
# Run from the repository root:  bash scripts/run_main_chain.sh
set -euo pipefail
[ -f config/pipeline.yaml ] || { echo "Run from the repository root"; exit 1; }
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
if [ -n "$(ls -A "results/01_callset/$VER/final" 2>/dev/null)" ]; then
  echo "results/01_callset/$VER/final already exists: callset $VER has been built. Nothing submitted."; exit 1
fi
mkdir -p logs
echo "Callset version: $VER | commit: $(git rev-parse --short HEAD 2>/dev/null || echo none)"
Rscript scripts/01_callset/00_check_inputs.R
J1=$(sbatch --parsable slurm/01_callset.sbatch)
J8=$(sbatch --parsable --dependency=afterok:"$J1" slurm/08_descriptive.sbatch)
J2=$(sbatch --parsable --dependency=afterok:"$J1" slurm/02_matrices.sbatch)
J2V=$(sbatch --parsable --dependency=afterok:"$J2" slurm/02_verify_matrices.sbatch)
J3=$(sbatch --parsable --dependency=afterok:"$J2" slurm/03_meff.sbatch)
J4=$(sbatch --parsable --dependency=afterok:"$J2" slurm/04_gwas.sbatch)
J5=$(sbatch --parsable --dependency=afterok:"$J3":"$J4" slurm/03_04_checks.sbatch)
J6=$(sbatch --parsable --dependency=afterok:"$J5" --export=ALL,RELS=all --job-name=signals_all slurm/05_signals.sbatch)
J7=$(sbatch --parsable --dependency=afterok:"$J5" --export=ALL,RELS=unrelated --job-name=signals_unrel slurm/05_signals.sbatch)
{
  echo "Main chain for callset $VER submitted $(date) from commit $(git rev-parse --short HEAD 2>/dev/null || echo none)"
  printf "%-22s %s\n" "1 callset" "$J1" "8 descriptive" "$J8" "2 matrices" "$J2" "2 verify matrices" "$J2V" \
         "3 Meff (array)" "$J3" "4 GWAS (array)" "$J4" "3-4 checks" "$J5" "5 signals all" "$J6" "5 signals unrelated" "$J7"
  echo "Logs: logs/<job name>_<job id>.log (arrays: logs/<job name>_<array id>_<task>.log)"
} | tee "logs/chain_${VER}.txt"
