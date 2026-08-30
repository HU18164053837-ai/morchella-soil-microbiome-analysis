#!/usr/bin/env bash
#SBATCH --job-name=mor_fungal_func
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=04:00:00
#SBATCH --output=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_fungal_func_%j.out
#SBATCH --error=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_fungal_func_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome; EXT=$BASE/extension_analysis_20260822; DADAIMG=$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif; mkdir -p "$EXT/results/funguild_v3"
module load apps/apptainer/1.3.4
apptainer exec --cleanenv --bind "$BASE:$BASE" "$DADAIMG" python "$EXT/scripts/04_funguild_offline.py" "$BASE/results/its_unite_v2/its_asv_taxonomy.tsv" "$EXT/software/funguild/funguild_db_response.txt" "$EXT/results/funguild_v3/ITS_ASV_FUNGuild_offline.tsv"
test -s "$EXT/results/funguild_v3/ITS_ASV_FUNGuild_offline.tsv"
echo MORCHELLA_FUNGUILD_COMPLETE

