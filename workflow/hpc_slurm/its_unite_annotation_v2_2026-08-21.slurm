#!/usr/bin/env bash
#SBATCH --job-name=mor_its_unite2
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/its_unite_v2_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/its_unite_v2_%j.err

set -euo pipefail

PROJECT=/path/to/morchella_microbiome
UNITE="$PROJECT/databases/unite_v10_20250219_fungi2"
QUERY="$PROJECT/results/dada2_ITS_v1/asv_sequences.fasta"
HITS="$PROJECT/results/its_unite_v1/unite_blast_hits.tsv"
OUT="$PROJECT/results/its_unite_v2"
PARSER="$PROJECT/scripts/15_parse_unite_blast_v2.py"
STAGE="$PROJECT/work/its_unite_v2_${SLURM_JOB_ID}"

module load apps/python/3.8.10
test -s "$QUERY"
test -s "$HITS"
test -s "$PARSER"
test -s "$UNITE/unite_database_manifest.tsv"
test -s "$UNITE/SHA256SUMS.txt"
test ! -e "$OUT"
test ! -e "$STAGE"
mkdir -p "$STAGE"

python3 "$PARSER" "$HITS" "$QUERY" "$STAGE"
test "$(awk 'END{print NR-1}' "$STAGE/its_asv_taxonomy.tsv")" -eq 8056
test "$(awk -F '\t' '$1=="annotated_kingdom"{print $2}' "$STAGE/taxonomy_annotation_summary.tsv")" -gt 0
test "$(awk -F '\t' '$1=="annotated_genus"{print $2}' "$STAGE/taxonomy_annotation_summary.tsv")" -gt 0

cp "$HITS" "$STAGE/unite_blast_hits.tsv"
cp "$UNITE/unite_database_manifest.tsv" "$STAGE/unite_database_manifest.tsv"
cp "$UNITE/SHA256SUMS.txt" "$STAGE/UNITE_SHA256SUMS.txt"
python3 --version > "$STAGE/pipeline_versions.txt" 2>&1
printf 'BLAST source job\t119637542\n' >> "$STAGE/pipeline_versions.txt"
printf 'Parser fix\taccept pipe or semicolon before UNITE rank prefixes\n' >> "$STAGE/pipeline_versions.txt"
sha256sum "$STAGE"/* > "$STAGE/result_SHA256SUMS.txt"

mv "$STAGE" "$OUT"
echo "ITS_UNITE_V2_ANNOTATION_COMPLETE"

