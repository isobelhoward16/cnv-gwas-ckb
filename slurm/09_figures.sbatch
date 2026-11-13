#!/bin/bash
#SBATCH --job-name=figures
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=08:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# Stage 9: Figure 3 (heatmap), Figure 4A (16p11.2 zoom), boxplots for every stage 5 signal (Figures 4B-C, 5B, 5D
# and the supplement), Figure 5A and 5C zooms (1q21.1 standing height, mirror; 8p23.1 random glucose, deletion), and zooms
# of the other loci (7p21.1, 8p22, 12p13.31, 15q11.2).
# Each step runs even if an earlier one fails; the job fails at the end if any step did. Submitted by scripts/run_post_stage5.sh.
# Set up R here (packages listed in the README), e.g. module load R
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
echo "Callset version: $(awk '/^callset_version:/{print $2}' config/pipeline.yaml) | commit: $(git rev-parse --short HEAD) | $(date)"
FAILED=""
run() { echo; echo "==================== $* | $(date)"; Rscript "$@" || FAILED="$FAILED $(basename "$1")"; }
run scripts/09_figures/01_figure3_heatmap.R
run scripts/09_figures/03_figure4a_16p11_zoom.R
run scripts/09_figures/05_figure5_zoom.R 1 mirror standing_height_mm both 1q21_1_standing_height_mirror_stacked_pubstyle
run scripts/09_figures/05_figure5_zoom.R 8 deletion random_glucose deletion 8p23_1_random_glucose_del_stacked_pubstyle
# zooms of the other loci in the same style; model = the Table 1 model
run scripts/09_figures/05_figure5_zoom.R 7 deletion heart_rate_mean deletion 7p21_1_heart_rate_del_stacked_pubstyle
run scripts/09_figures/05_figure5_zoom.R 8 mirror sitting_height_mm both 8p22_sitting_height_mirror_stacked_pubstyle
run scripts/09_figures/05_figure5_zoom.R 12 duplication standing_height_mm duplication 12p13_31_standing_height_dup_stacked_pubstyle
run scripts/09_figures/05_figure5_zoom.R 15 deletion standing_height_mm deletion 15q11_2_standing_height_del_stacked_pubstyle
run scripts/09_figures/04_boxplots.R
echo; echo "Finished $(date) | failed:${FAILED:- none}"
[ -z "$FAILED" ]
