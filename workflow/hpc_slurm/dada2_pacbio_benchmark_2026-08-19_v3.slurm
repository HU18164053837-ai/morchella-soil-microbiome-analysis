#!/usr/bin/env bash
#SBATCH --job-name=mor_dada2_bench2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_benchmark_v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_benchmark_v2_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"

module load apps/apptainer/1.3.4
cd "$PROJECT"

test -s "$CONTAINER"
test -s "$PROJECT/input/trimmed/16S/A2-2_16S_trimmed.fastq.gz"
test -s "$PROJECT/input/trimmed/ITS/A2-2_ITS_trimmed.fastq.gz"
test -s "$PROJECT/scripts/05_dada2_pacbio_benchmark_v2.R"
test ! -e "$PROJECT/results/benchmark_a2_2_full_v2"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/05_dada2_pacbio_benchmark_v2.R /project


