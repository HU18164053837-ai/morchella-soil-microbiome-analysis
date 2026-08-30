#!/usr/bin/env bash
#SBATCH --job-name=mor_16s_band32
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_16s_band32_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_16s_band32_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"

module load apps/apptainer/1.3.4
cd "$PROJECT"
test -s "$CONTAINER"
test -s "$PROJECT/scripts/06_dada2_16s_band32_diagnostic.R"
test -s "$PROJECT/results/benchmark_a2_2_full_v2/16S/error_model_benchmark.rds"
test ! -e "$PROJECT/results/diagnostic_16s_band32_v1"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/06_dada2_16s_band32_diagnostic.R /project "$SLURM_CPUS_PER_TASK"

