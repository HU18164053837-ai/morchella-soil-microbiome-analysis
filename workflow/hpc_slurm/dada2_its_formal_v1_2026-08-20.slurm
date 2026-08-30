#!/usr/bin/env bash
#SBATCH --job-name=mor_its_formal1
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=36:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_its_formal_v1_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_its_formal_v1_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"
SCRIPT="$PROJECT/scripts/11_dada2_its_formal_v1.R"
WORKDIR="$PROJECT/work/dada2_ITS_v1"
OUTDIR="$PROJECT/results/dada2_ITS_v1"
CHECKPOINT="$PROJECT/checkpoints/dada2_ITS_v1"

module load apps/apptainer/1.3.4
cd "$PROJECT"
test -s "$CONTAINER"
test -s "$SCRIPT"
test "$(find "$PROJECT/input/trimmed/ITS" -maxdepth 1 -type f -name '*_ITS_trimmed.fastq.gz' | wc -l)" -eq 24
test ! -e "$WORKDIR"
test ! -e "$OUTDIR"
test ! -e "$CHECKPOINT"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/11_dada2_its_formal_v1.R \
    /project "$SLURM_CPUS_PER_TASK"

