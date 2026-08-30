#!/bin/bash
#SBATCH --job-name=mor_pv_16s_dada2v2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/path/to/morchella_microbiome/public_validation/logs/mor_pv_16s_dada2v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/public_validation/logs/mor_pv_16s_dada2v2_%j.err
set -euo pipefail
module load apps/apptainer/1.3.4
ROOT=/path/to/morchella_microbiome
PV=$ROOT/public_validation
IMG=$ROOT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif
OUT=$PV/projects/PRJNA993383/dada2/dada2_single_end_v1
test -s "$IMG"
test -s "$PV/scripts/06_dada2_single_end_public_v1.R"
test "$(find "$PV/projects/PRJNA993383/raw/fastq" -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)" -eq 24
test ! -e "$OUT"
apptainer exec --cleanenv --bind "$ROOT:$ROOT" "$IMG" Rscript \
  "$PV/scripts/06_dada2_single_end_public_v1.R" \
  PRJNA993383 16S "$PV/projects/PRJNA993383/raw/fastq" "$OUT" 240 270 5


