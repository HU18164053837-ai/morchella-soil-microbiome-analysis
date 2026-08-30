#!/usr/bin/env bash
#SBATCH --job-name=mor_16s_A34qc
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_16s_A3_A4_quality_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_16s_A3_A4_quality_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"
SCRIPT="$PROJECT/scripts/09_dada2_16s_A3_A4_quality_diagnostic.R"
OUTDIR="$PROJECT/results/diagnostic_16s_A3_A4_quality_v1"

module load apps/apptainer/1.3.4
cd "$PROJECT"
test -s "$CONTAINER"
test -s "$SCRIPT"
test -s "$PROJECT/input/trimmed/16S/A3-1_16S_trimmed.fastq.gz"
test -s "$PROJECT/input/trimmed/16S/A4-1_16S_trimmed.fastq.gz"
test -s "$PROJECT/results/diagnostic_16s_balanced_A_v1/filtered/A3-1_16S_filtered.fastq.gz"
test -s "$PROJECT/results/diagnostic_16s_balanced_A_v1/filtered/A4-1_16S_filtered.fastq.gz"
test ! -e "$OUTDIR"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/09_dada2_16s_A3_A4_quality_diagnostic.R \
    /project "$SLURM_CPUS_PER_TASK"

