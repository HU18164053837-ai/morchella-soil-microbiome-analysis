#!/usr/bin/env bash
#SBATCH --job-name=mor_16s_balA
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_16s_balanced_A_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_16s_balanced_A_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"

module load apps/apptainer/1.3.4
cd "$PROJECT"
test -s "$CONTAINER"
test -s "$PROJECT/scripts/07_dada2_16s_balanced_A_diagnostic.R"
for sample in A2-1 A3-1 A4-1; do
  test -s "$PROJECT/input/trimmed/16S/${sample}_16S_trimmed.fastq.gz"
done
test ! -e "$PROJECT/results/diagnostic_16s_balanced_A_v1"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/07_dada2_16s_balanced_A_diagnostic.R \
    /project "$SLURM_CPUS_PER_TASK"

