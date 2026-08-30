#!/bin/bash
#SBATCH --job-name=mor_pv_its_dada2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/path/to/morchella_microbiome/public_validation/logs/mor_pv_its_dada2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/public_validation/logs/mor_pv_its_dada2_%j.err
set -euo pipefail
ROOT=/path/to/morchella_microbiome
PV=$ROOT/public_validation
IMG=$ROOT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif
OUT=$PV/projects/PRJNA935967/dada2/dada2_single_end_v1
test ! -e "$OUT"
singularity exec --bind "$ROOT:$ROOT" "$IMG" Rscript \
  "$PV/scripts/06_dada2_single_end_public_v1.R" \
  PRJNA935967 ITS "$PV/projects/PRJNA935967/raw/fastq" "$OUT" 100 400 5


