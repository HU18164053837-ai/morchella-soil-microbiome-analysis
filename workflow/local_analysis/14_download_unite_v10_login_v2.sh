#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/morchella_microbiome
TARGET="$PROJECT/databases/unite_v10_20250219_fungi2"
TMP="$TARGET.tmp_download"
URL="https://s3.hpc.ut.ee/plutof-public/original/1c8647cd-31f2-44ff-8c77-dbd470e6154c.tgz"
ARCHIVE="$TMP/sh_general_release_s_19.02.2025.tgz"
EXPECTED_SIZE=38985832

test ! -e "$TARGET"
mkdir -p "$TMP"
if test -s "$ARCHIVE"; then
  test "$(stat -c %s "$ARCHIVE")" -eq "$EXPECTED_SIZE"
else
  curl -L --fail --retry 20 --retry-delay 5 --connect-timeout 20 \
    -C - -o "$ARCHIVE.part" "$URL"
  test "$(stat -c %s "$ARCHIVE.part")" -eq "$EXPECTED_SIZE"
  mv "$ARCHIVE.part" "$ARCHIVE"
fi

if tar -tzf "$ARCHIVE" | awk '/^\// || /(^|\/)\.\.(\/|$)/ { bad=1 } END { exit bad }'; then
  :
else
  echo "Unsafe archive member detected" >&2
  exit 1
fi

mkdir -p "$TMP/extracted"
tar -xzf "$ARCHIVE" -C "$TMP/extracted"
mapfile -t FASTA_FILES < <(find "$TMP/extracted" -type f \
  \( -iname '*.fasta' -o -iname '*.fa' -o -iname '*.fna' \) | sort)
test "${#FASTA_FILES[@]}" -eq 1
FASTA="${FASTA_FILES[0]}"
RECORDS=$(grep -c '^>' "$FASTA")
test "$RECORDS" -ge 150000
test "$RECORDS" -le 200000

sha256sum "$ARCHIVE" "$FASTA" > "$TMP/SHA256SUMS.txt"
{
  printf 'field\tvalue\n'
  printf 'doi\t10.15156/BIO/3301230\n'
  printf 'release\tUNITE general FASTA release for Fungi 2 v10.0 2025-02-19\n'
  printf 'license\tCC BY\n'
  printf 'source_url\t%s\n' "$URL"
  printf 'archive_bytes\t%s\n' "$EXPECTED_SIZE"
  printf 'fasta_records\t%s\n' "$RECORDS"
  printf 'downloaded_utc\t%s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
  printf 'fasta_relative_path\t%s\n' "${FASTA#$TMP/}"
} > "$TMP/unite_database_manifest.tsv"

mv "$TMP" "$TARGET"
echo "UNITE_DOWNLOAD_COMPLETE target=$TARGET records=$RECORDS"

