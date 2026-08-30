#!/usr/bin/env bash
#SBATCH --job-name=mor_16s_balA2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/dada2_16s_balanced_A_v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/dada2_16s_balanced_A_v2_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"
SOURCE="$PROJECT/results/diagnostic_16s_balanced_A_v1"

module load apps/apptainer/1.3.4
cd "$PROJECT"
test -s "$CONTAINER"
test -s "$PROJECT/scripts/08_dada2_16s_balanced_A_resume_v2.R"
test -s "$SOURCE/error_model_A_balanced.rds"
test -s "$SOURCE/dada_A_independent.rds"
test -s "$SOURCE/seqtab_A_independent.rds"
test -s "$SOURCE/filter_stats.tsv"
test ! -e "$PROJECT/results/diagnostic_16s_balanced_A_v2"

apptainer exec --cleanenv \
  --bind "$PROJECT:/project" \
  "$CONTAINER" \
  Rscript /project/scripts/08_dada2_16s_balanced_A_resume_v2.R \
    /project "$SLURM_CPUS_PER_TASK"

