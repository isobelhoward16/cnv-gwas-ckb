# CNV association analysis in the China Kadoorie Biobank

Code for the paper *Copy number variant association analysis in 94,725 Chinese adults reveals loci influencing
anthropometric and cardiometabolic traits* (Howard et al.).

CNVs were called with PennCNV-Affy from Axiom CKB1 and CKB2 genotyping array data and tested for association with
13 anthropometric and cardiometabolic traits under deletion-only, duplication-only and mirror models.

No individual-level data are included. Paths to the data are set in `config/paths.yaml` (not committed; copy
`config/paths.example.yaml`). Analysis settings are in `config/pipeline.yaml`. All outputs are written to `results/`.

## Contents

| Folder | Content |
|---|---|
| `calling/` | PennCNV-Affy calling, one job per region and array, in batches of 1,000 samples; per-sample QC summaries |
| `R/` | Shared functions (configuration, task list, call merging) |
| `scripts/01_callset` | Merging of fragmented calls, sample QC, call filters, removal of samples with more than 100 CNVs, gene annotation |
| `scripts/02_matrices` | Deletion, duplication and mirror probe matrices on the probes shared by both arrays |
| `scripts/03_meff` | Effective number of independent tests (simpleM) and significance thresholds |
| `scripts/04_gwas` | Probe-level association tests (linear regression on rank-inverse-normal transformed traits) |
| `scripts/05_signals` | Significant probes, CNV regions (CNVRs), annotation and relationships between signals |
| `scripts/06_followup` | Lead probes in unrelated participants, array-stratified checks |
| `scripts/07_external` | UK Biobank comparison, GWAS Catalog look-ups, segmental duplications |
| `scripts/08_descriptive` | CNV landscape, genomic-disorder CNVs, analysis population, flowchart counts |
| `scripts/09_figures` | Figures 3-5 and Figures S3-S5 |
| `scripts/10_paper` | Table 1, Tables S1, S5 and S6, numbers quoted in the text |
| `slurm/` | SLURM job scripts for each stage |

`config/paper_outputs.tsv` lists which output file is each figure and table in the paper.

## Running

All commands are run from the repository root, on a SLURM cluster. Set up R in the job scripts as needed
(each has a line marking where).

1. CNV calling. Inputs are the per-sample LRR/BAF signal files, a population B allele frequency (PFB) file for each
   region and array, and a GC model file for each array (see Methods). Set the paths in `calling/config.sh`, then
   `mkdir -p logs; sbatch calling/batch_affyhmm_launcher.sh`, and when all jobs have finished,
   `sbatch calling/generate_qc_sum.sh`.
2. Stages 1-5 and 8: `bash scripts/run_main_chain.sh` (jobs are submitted with dependencies).
3. Stage 6: `bash scripts/06_followup/run_unrelated_leads.sh` and `sbatch slurm/06_fasting_rpg.sbatch`.
4. Stages 6 (array checks), 7, 9 and 10: `bash scripts/run_post_stage5.sh <unrelated tables job id> <fasting job id>`.

## Requirements

R (>= 4.1) with data.table, Matrix, yaml, ggplot2, dplyr, tidyr, stringr, scales, cowplot, patchwork, gridExtra,
gtable, ggrepel, ggbeeswarm, HGNChelper, and the Bioconductor packages GenomicRanges, IRanges and rtracklayer.
Perl and PennCNV 1.0.5. External files: RefSeq genes, GENCODE v19 and cytobands (hg19), the UCSC GWAS Catalog
table, UCSC genomicSuperDups (hg19), and the supplementary data of Auwerx et al. 2022 (UK Biobank CNV-GWAS).

## Data availability

China Kadoorie Biobank data are available to bona fide researchers; see https://www.ckbiobank.org/data-access.
