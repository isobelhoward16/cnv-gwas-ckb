#!/bin/bash
#SBATCH --job-name=paper
#SBATCH --ntasks=1
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --partition=general
#SBATCH --output=logs/%x_%j.log
# After stages 6, 7 and 9: Figures S3-S5 (array checks) and S7 (unrelated), Table 1, Tables S1/S5/S6, the manuscript
# numbers, then the collection of every final file in results/paper_final/<version>/ (MANIFEST.tsv) and a copy of all
# signal boxplots in results/paper_final/<version>/boxplots_all/ and of all zooms in zooms_all/. Each step runs even if an earlier one fails.
# Submitted by scripts/run_post_stage5.sh (afterany: it reports what is missing rather than not running).
# Set up R here (packages listed in the README), e.g. module load R
set -uo pipefail
cd "$SLURM_SUBMIT_DIR"
VER=$(awk '/^callset_version:/{print $2}' config/pipeline.yaml)
echo "Callset version: $VER | commit: $(git rev-parse --short HEAD) | $(date)"
FAILED=""
run() { echo; echo "==================== $* | $(date)"; Rscript "$@" || FAILED="$FAILED $(basename "$1")"; }
run scripts/09_figures/02_figures_s3_s5_array.R
run scripts/06_followup/03_figure_s7_unrelated.R
run scripts/10_paper/00_table1_carriers_both_types.R
run scripts/10_paper/01_table1.R
run scripts/10_paper/02_tables_s1_s5_s6.R
run scripts/10_paper/03_manuscript_numbers.R
run scripts/10_paper/05_text_support.R
run scripts/10_paper/06_landscape_numbers.R
run scripts/10_paper/09_collect_paper_final.R
B=results/paper_final/$VER/boxplots_all; rm -rf "$B"; mkdir -p "$B"
cp results/09_figures/$VER/boxplots/unpruned_related_withpoints/boxplot_*.png "$B"/ 2>/dev/null || FAILED="$FAILED boxplots_copy"
echo; echo "Boxplots copied to $B: $(ls "$B" | wc -l) (expected 23, one per stage 5 signal)"
Z=results/paper_final/$VER/zooms_all; rm -rf "$Z"; mkdir -p "$Z"
cp results/09_figures/$VER/zooms/*.png "$Z"/ 2>/dev/null || FAILED="$FAILED zooms_copy"
echo "Zooms copied to $Z: $(ls "$Z" | wc -l) (16p11.2 + one per other locus)"
echo "Finished $(date) | failed:${FAILED:- none}"
[ -z "$FAILED" ]
