#!/bin/bash
# After stage 5: array checks (stage 6), external comparisons (stage 7) and figures (stage 9) in parallel, then the
# paper job (stage 10) once those and the already-submitted unrelated-leads and fasting-glucose jobs have finished.
# Usage (repository root):  bash scripts/run_post_stage5.sh <unrelated tables job id> <fasting job id>
set -eo pipefail
[ $# -eq 2 ] || { echo "Usage: bash scripts/run_post_stage5.sh <unrelated tables job id> <fasting job id>"; exit 1; }
mkdir -p logs
A=$(sbatch --parsable slurm/06_array_checks.sbatch)
E=$(sbatch --parsable slurm/07_external.sbatch)
F=$(sbatch --parsable slurm/09_figures.sbatch)
P=$(sbatch --parsable --dependency=afterany:"$A":"$E":"$F":"$1":"$2" slurm/10_paper.sbatch)
echo "array_checks $A | external $E | figures $F | paper $P (after $A $E $F $1 $2)" | tee -a logs/chain_post_stage5.txt
echo "When the queue is empty:  tail -60 logs/array_checks_$A.log; cat logs/external_$E.log | tail -150; tail -80 logs/figures_$F.log; cat logs/paper_$P.log"
