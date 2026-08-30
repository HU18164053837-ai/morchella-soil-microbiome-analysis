#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/morchella_microbiome
TARGET="$PROJECT/databases/unite_v10_20250219_fungi2"
TMP="$TARGET.tmp_download"
URL="https://s3.hpc.ut.ee/plutof-public/original/1c8647cd-31f2-44ff-8c77-dbd470e6154c.tgz"
ARCHIVE="$TMP/sh_general_release_s_19.02.2025.tgz"
EXPECTED_SIZE=38985832
FASTA="$TMP/extracted/sh_general_release_dynamic_s_19.02.2025.fasta"
DEV_FASTA="$TMP/extracted/sh_general_release_dynamic_s_19.02.2025_dev.fasta"

test ! -e "$TARGET"
test -s "$ARCHIVE"
test "$(stat -c %s "$ARCHIVE")" -eq "$EXPECTED_SIZE"

if tar -tzf "$ARCHIVE" | awk '/^\// || /(^|\/)\.\.(\/|$)/ { bad=1 } END { exit bad }'; then
  :
else
  echo "Unsafe archive member detected" >&2
  exit 1
fi

# The official archive contains both the standard dynamic FASTA and a _dev
# variant. Use the standard file for production taxonomy annotation.
if test ! -s "$FASTA" || test ! -s "$DEV_FASTA"; then
  test ! -e "$TMP/extracted"
  mkdir -p "$TMP/extracted"
  tar -xzf "$ARCHIVE" -C "$TMP/extracted"
fi
test -s "$FASTA"
test -s "$DEV_FASTA"
RECORDS=$(grep -c '^>' "$FASTA")
DEV_RECORDS=$(grep -c '^>' "$DEV_FASTA")
test "$RECORDS" -eq 168030
test "$DEV_RECORDS" -eq 168030

sha256sum "$ARCHIVE" "$FASTA" "$DEV_FASTA" > "$TMP/SHA256SUMS.txt"
{
  printf 'field\tvalue\n'
  printf 'doi\t10.15156/BIO/3301230\n'
  printf 'release\tUNITE general FASTA release for Fungi 2 v10.0 2025-02-19\n'
  printf 'license\tCC BY\n'
  printf 'source_url\t%s\n' "$URL"
  printf 'archive_bytes\t%s\n' "$EXPECTED_SIZE"
  printf 'fasta_records\t%s\n' "$RECORDS"
  printf 'dev_fasta_records\t%s\n' "$DEV_RECORDS"
  printf 'selected_variant\tstandard dynamic (non-_dev)\n'
  printf 'downloaded_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'fasta_relative_path\t%s\n' "${FASTA#$TMP/}"
} > "$TMP/unite_database_manifest.tsv"

mv "$TMP" "$TARGET"
echo "UNITE_DOWNLOAD_COMPLETE target=$TARGET records=$RECORDS variant=standard_dynamic"

