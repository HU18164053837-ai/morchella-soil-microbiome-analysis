#!/usr/bin/env bash
#SBATCH --job-name=mor_picrust2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_picrust2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_picrust2_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome; EXT=$BASE/extension_analysis_20260822; IMG=$EXT/containers/picrust2_2.5.3.sif; DADAIMG=$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif; mkdir -p "$EXT/containers"
module load apps/apptainer/1.3.4
apptainer exec --cleanenv --bind "$BASE:$BASE" "$DADAIMG" python "$EXT/scripts/02_prepare_picrust_input.py" "$EXT/input/16S_clean_asv_count_table.tsv" "$BASE/results/dada2_16S_v1/asv_sequences.fasta" "$EXT/input/16S_clean_asv_sequences.fasta"
if [ ! -s "$IMG" ]; then apptainer pull "$IMG" docker://quay.io/biocontainers/picrust2:2.5.3--pyhdfd78af_0; fi
apptainer exec --cleanenv --bind "$BASE:$BASE" "$IMG" picrust2_pipeline.py -s "$EXT/input/16S_clean_asv_sequences.fasta" -i "$EXT/input/16S_clean_asv_count_table.tsv" -o "$EXT/results/picrust2_v4" -p 10
test -s "$EXT/results/picrust2_v4/pathways_out/path_abun_unstrat.tsv.gz"
echo MORCHELLA_PICRUST2_COMPLETE

