# Paths for the CNV calling scripts. Sourced by every script in calling/ (submit them from the repository root).
PENNCNV_DIR=/path/to/PennCNV-1.0.5
PFB_DIR=/path/to/pfb                    # one PFB file per region and array: <Region>_<Phase>_stage2.pfb
GC_MODEL_DIR=/path/to/gc_model_files    # full_gc_output_array1.txt (CKB1), full_gc_output_array2.txt (CKB2)
FILE_LIST_DIR=/path/to/file_lists       # <Region>_<Phase>_file_list.txt: one signal file per line, as <archive>.zip/<file>
CALL_DIR=/path/to/penncnv_calls         # output: one folder per region and array
TEMP_ROOT=/path/to/tmp                  # signal files are unzipped here, one batch at a time
