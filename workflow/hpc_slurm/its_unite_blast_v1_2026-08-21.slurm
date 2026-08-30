#!/usr/bin/env bash
#SBATCH --job-name=mor_its_unite1
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=16G
#SBATCH --time=08:00:00
#SBATCH --output=/path/to/morchella_microbiome/logs/its_unite_blast_v1_%j.out
#SBATCH --error=/path/to/morchella_microbiome/logs/its_unite_blast_v1_%j.err

set -euo pipefail
PROJECT=/path/to/morchella_microbiome
UNITE="$PROJECT/databases/unite_v10_20250219_fungi2"
QUERY="$PROJECT/results/dada2_ITS_v1/asv_sequences.fasta"
WORK="$PROJECT/work/its_unite_v1"
OUT="$PROJECT/results/its_unite_v1"
PARSER="$PROJECT/scripts/15_parse_unite_blast.py"

module load apps/blast/2.10.0
test -s "$QUERY"
test -s "$PARSER"
test -s "$UNITE/unite_database_manifest.tsv"
test ! -e "$WORK"
test ! -e "$OUT"
mkdir -p "$WORK/db" "$OUT"

FASTA_REL=$(awk -F '\t' '$1=="fasta_relative_path"{print $2}' "$UNITE/unite_database_manifest.tsv")
REFERENCE="$UNITE/$FASTA_REL"
test -s "$REFERENCE"
makeblastdb -in "$REFERENCE" -dbtype nucl -out "$WORK/db/unite_v10"

blastn -task blastn -query "$QUERY" -db "$WORK/db/unite_v10" \
  -num_threads "$SLURM_CPUS_PER_TASK" -evalue 1e-20 -perc_identity 80 \
  -qcov_hsp_perc 70 -max_target_seqs 20 -dust no \
  -outfmt '6 qseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs' \
  -out "$OUT/unite_blast_hits.tsv"

python3 "$PARSER" "$OUT/unite_blast_hits.tsv" "$QUERY" "$OUT"
cp "$UNITE/unite_database_manifest.tsv" "$OUT/unite_database_manifest.tsv"
cp "$UNITE/SHA256SUMS.txt" "$OUT/UNITE_SHA256SUMS.txt"
blastn -version > "$OUT/pipeline_versions.txt"
python3 --version >> "$OUT/pipeline_versions.txt" 2>&1
test "$(awk 'END{print NR-1}' "$OUT/its_asv_taxonomy.tsv")" -eq 8056
echo "ITS_UNITE_ANNOTATION_COMPLETE"

