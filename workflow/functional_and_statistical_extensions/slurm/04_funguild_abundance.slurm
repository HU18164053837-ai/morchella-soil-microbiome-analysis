#!/usr/bin/env bash
#SBATCH --job-name=mor_fg_abund
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_fg_abund_%j.out
#SBATCH --error=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_fg_abund_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome; EXT=$BASE/extension_analysis_20260822; IMG=$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif; OUT=$EXT/results/funguild_abundance_v1
module load apps/apptainer/1.3.4
mkdir -p "$OUT"
apptainer exec --cleanenv --bind "$BASE:$BASE" "$IMG" python "$EXT/scripts/05_funguild_abundance_paired.py" "$BASE/results/dada2_ITS_v1/asv_count_table.tsv" "$EXT/results/funguild_v3/ITS_ASV_FUNGuild_offline.tsv" "$OUT"
test -s "$OUT/FUNGuild_paired_function_statistics.tsv"
echo MORCHELLA_FUNGUILD_ABUNDANCE_COMPLETE

