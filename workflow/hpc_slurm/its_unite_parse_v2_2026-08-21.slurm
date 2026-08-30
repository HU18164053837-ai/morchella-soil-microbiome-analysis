#!/usr/bin/env bash
#SBATCH --job-name=mor_its_uparse2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/its_unite_parse_v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/its_unite_parse_v2_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
UNITE="$PROJECT/databases/unite_v10_20250219_fungi2"
QUERY="$PROJECT/results/dada2_ITS_v1/asv_sequences.fasta"
OUT="$PROJECT/results/its_unite_v1"
HITS="$OUT/unite_blast_hits.tsv"
PARSER="$PROJECT/scripts/15_parse_unite_blast.py"
STAGE="$PROJECT/work/its_unite_parse_v2_${SLURM_JOB_ID}"

module load apps/python/3.8.10
test -s "$QUERY"
test -s "$HITS"
test -s "$PARSER"
test -s "$UNITE/unite_database_manifest.tsv"
test -s "$UNITE/SHA256SUMS.txt"
test ! -e "$STAGE"
mkdir -p "$STAGE"

python3 "$PARSER" "$HITS" "$QUERY" "$STAGE"
test "$(awk 'END{print NR-1}' "$STAGE/its_asv_taxonomy.tsv")" -eq 8056
test -s "$STAGE/taxonomy_annotation_summary.tsv"

python3 --version > "$STAGE/pipeline_versions.txt" 2>&1
printf 'BLAST source job\t119637542\n' >> "$STAGE/pipeline_versions.txt"
printf 'BLAST hit rows\t%s\n' "$(wc -l < "$HITS")" >> "$STAGE/pipeline_versions.txt"
cp "$UNITE/unite_database_manifest.tsv" "$STAGE/unite_database_manifest.tsv"
cp "$UNITE/SHA256SUMS.txt" "$STAGE/UNITE_SHA256SUMS.txt"

for name in its_asv_taxonomy.tsv unclassified_asvs.fasta \
  taxonomy_annotation_summary.tsv pipeline_versions.txt \
  unite_database_manifest.tsv UNITE_SHA256SUMS.txt; do
  test ! -e "$OUT/$name"
  mv "$STAGE/$name" "$OUT/$name"
done
rmdir "$STAGE"

echo "ITS_UNITE_PARSE_V2_COMPLETE"

