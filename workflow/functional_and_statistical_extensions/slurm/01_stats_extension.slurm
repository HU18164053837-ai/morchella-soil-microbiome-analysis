#!/usr/bin/env bash
#SBATCH --job-name=mor_stats_ext
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_stats_ext_%j.out
#SBATCH --error=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_stats_ext_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome; EXT=$BASE/extension_analysis_20260822; IMG=$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif
module load apps/apptainer/1.3.4
mkdir -p "$EXT/results/stats_extension_v2"
apptainer exec --cleanenv --bind "$BASE:$BASE" "$IMG" Rscript "$EXT/scripts/01_stats_extension.R" "$BASE" "$EXT/results/stats_extension_v2"
test -s "$EXT/results/stats_extension_v2/16S_exact_paired_beta.tsv"
test -s "$EXT/results/stats_extension_v2/ITS_candidate_clr_sensitivity.tsv"
echo MORCHELLA_STATS_EXTENSION_COMPLETE

