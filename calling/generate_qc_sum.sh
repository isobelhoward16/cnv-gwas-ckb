#!/bin/bash
#SBATCH --job-name=generate_qcsum_affyhmm
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --output=logs/%x_%j.log
#SBATCH -p general
# Per-sample QC summary (LRR SD, BAF drift, waviness factor) for each region and array, from the PennCNV log
# (filter_cnv.pl -qclogfile). Output: <region_phase>/logs/<region_phase>_qcsum.txt
# Usage (repository root):  sbatch calling/generate_qc_sum.sh
source calling/config.sh

echo "Generating QC summaries..."

for dir in "$CALL_DIR"/*_Phase*; do
  [ -d "$dir" ] || continue
  rp=$(basename "$dir")
  log_path="$dir/logs/${rp}_affyhmm_cnv.log"
  out_qcsum="$dir/logs/${rp}_qcsum.txt"

  # filter_cnv.pl needs a CNV file as input; any batch file of this stratum is used
  cnv_file=$(find "$dir" -maxdepth 1 -name "${rp}_cnv_calls_batch_*.txt" | head -n 1)

  if [[ ! -f "$log_path" || ! -f "$cnv_file" ]]; then
    echo "Skipping $rp: missing log or CNV file"
    continue
  fi

  echo "Processing: $rp"

  "$PENNCNV_DIR"/filter_cnv.pl "$cnv_file" \
    -qclogfile "$log_path" \
    -qcsumout "$out_qcsum" \
    -qcpassout /dev/null \
    -out /dev/null
done

echo "Done generating QC summary files."
