#!/bin/bash
#SBATCH --job-name=extract_cnv_affyhmm
#SBATCH --ntasks=2
#SBATCH --cpus-per-task=4
#SBATCH --mem=65G
#SBATCH --time=60:00:00
#SBATCH --output=logs/%x_%j.log
#SBATCH -p general
# Calls CNVs for one region and array in batches of 1,000 samples: unzips the signal files of a batch, submits
# cnv_calling_affyhmm.sh for it and waits for that job before starting the next batch.
# Usage: sbatch calling/extract_and_submit_cnv_calling.sh <Region> <Phase12|Phase3>
source calling/config.sh

REGION="$1"
PHASE="$2"
BATCH_SIZE=1000
TEMP_DIR="${TEMP_ROOT}/${REGION}_${PHASE}"
CNV_SCRIPT="calling/cnv_calling_affyhmm.sh"
OUTPUT_DIR="${CALL_DIR}/${REGION}_${PHASE}"
LOG_DIR="${OUTPUT_DIR}/logs"

mkdir -p "$TEMP_DIR" "$OUTPUT_DIR" "$LOG_DIR"

LISTFILE="${FILE_LIST_DIR}/${REGION}_${PHASE}_file_list.txt"
TOTAL_FILES=$(wc -l < "$LISTFILE")
PROCESSED=0

while [ "$PROCESSED" -lt "$TOTAL_FILES" ]; do
    echo "Processing batch starting at line $PROCESSED for $REGION $PHASE"
    rm -rf "$TEMP_DIR"/*
    mkdir -p "$TEMP_DIR"

    BATCH_FILES=$(sed -n "$((PROCESSED + 1)),$((PROCESSED + BATCH_SIZE))p" "$LISTFILE")
    [[ -z "$BATCH_FILES" ]] && break

    EXTRACTED_LIST="$TEMP_DIR/extracted_file_list.txt"
    rm -f "$EXTRACTED_LIST"

    for FILEPATH in $BATCH_FILES; do
        ZIP_PATH=$(echo "$FILEPATH" | sed -E 's|^(.*\.zip)/.*$|\1|')
        FILE_IN_ZIP=$(echo "$FILEPATH" | sed -E 's|^.*/([^/]+)$|\1|')

        if [[ -f "$ZIP_PATH" ]]; then
            unzip -j "$ZIP_PATH" "$FILE_IN_ZIP" -d "$TEMP_DIR"
            echo "$TEMP_DIR/$FILE_IN_ZIP" >> "$EXTRACTED_LIST"
        else
            echo "Warning: ZIP file $ZIP_PATH not found, skipping $FILEPATH" | tee -a "$LOG_DIR/extract_and_submit_affyhmm_${SLURM_JOB_ID}.log"
        fi
    done

    echo "Batch extraction completed for $REGION $PHASE."

    if [[ -s "$EXTRACTED_LIST" ]]; then
        CNV_JOB_ID=$(sbatch --parsable "$CNV_SCRIPT" "$REGION" "$PHASE" "$EXTRACTED_LIST" "$OUTPUT_DIR")
        echo "Submitted CNV calling job $CNV_JOB_ID for $REGION $PHASE"
        while squeue -j "$CNV_JOB_ID" | grep -q "$CNV_JOB_ID"; do sleep 30; done
        echo "CNV calling job $CNV_JOB_ID completed for $REGION $PHASE."
    fi

    PROCESSED=$((PROCESSED + BATCH_SIZE))
done

echo "All CNV calling batches for ${REGION} ${PHASE} completed. Results in $OUTPUT_DIR."
