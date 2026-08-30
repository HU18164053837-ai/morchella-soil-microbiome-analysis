#!/usr/bin/env bash
#SBATCH --job-name=mor_16s_silva1382_v2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/16s_silva1382_v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/16s_silva1382_v2_%j.err

set -euo pipefail
PROJECT=/path/to/morchella_microbiome
CONTAINER="$PROJECT/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif"
SCRIPT="$PROJECT/scripts/25_16s_silva_taxonomy_v2.R"
INPUT="$PROJECT/input/16S_corrected_v2"
DB="$PROJECT/databases/silva_138_2_dada2"
OUT="$PROJECT/results/16S_silva1382_v2"

module load apps/apptainer/1.3.4
test -s "$CONTAINER"
test -s "$SCRIPT"
test -s "$INPUT/seqtab_nochim.rds"
test -s "$INPUT/asv_count_table.tsv"
echo '1764e2a36b4500ccb1c7d5261948a414  '"$DB/silva_nr99_v138.2_toGenus_trainset.fa.gz" | md5sum -c -
echo '62fd939bd1aa9832d2089847f9f0b87d  '"$DB/silva_v138.2_assignSpecies.fa.gz" | md5sum -c -
test ! -e "$OUT"

cd "$PROJECT"
apptainer exec --cleanenv --bind "$PROJECT:/project" "$CONTAINER" \
  Rscript /project/scripts/25_16s_silva_taxonomy_v2.R \
  /project "$SLURM_CPUS_PER_TASK"


