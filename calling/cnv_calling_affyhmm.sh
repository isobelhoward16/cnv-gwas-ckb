#!/bin/bash
#SBATCH --job-name=cnv_calling_affyhmm
#SBATCH --ntasks=2
#SBATCH --mem=50G
#SBATCH --time=30:00:00
#SBATCH --output=logs/%x_%j.log
#SBATCH -p general
# PennCNV-Affy calling for one batch of samples from one region and array, with the region- and array-specific
# PFB file and the GC model of the array. Submitted by extract_and_submit_cnv_calling.sh.
# Usage: sbatch calling/cnv_calling_affyhmm.sh <Region> <Phase12|Phase3> <list of signal files> <output dir>
source calling/config.sh

REGION="$1"
PHASE="$2"
EXTRACTED_LIST="$3"
OUTPUT_DIR="$4"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$LOG_DIR"

PFBFILE="${PFB_DIR}/${REGION}_${PHASE}_stage2.pfb"
GCFILE="${GC_MODEL_DIR}/full_gc_output_array1.txt"
[[ "$PHASE" == "Phase3" ]] && GCFILE="${GC_MODEL_DIR}/full_gc_output_array2.txt"

perl "$PENNCNV_DIR"/detect_cnv.pl --test \
  --hmmfile "$PENNCNV_DIR/affy/libgw6/affygw6.hmm" \
  --pfbfile "$PFBFILE" \
  --gcmodel "$GCFILE" \
  --listfile "$EXTRACTED_LIST" \
  --output "$OUTPUT_DIR/${REGION}_${PHASE}_cnv_calls_batch_${SLURM_JOB_ID}.txt" \
  --logfile "$LOG_DIR/${REGION}_${PHASE}_affyhmm_cnv.log" \
  --conf

echo "CNV calling complete for $REGION $PHASE."
