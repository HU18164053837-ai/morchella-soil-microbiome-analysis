#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/morchella_microbiome/public_validation
LOG="$ROOT/logs/representative_download_$(date +%Y%m%d_%H%M%S).log"
MANIFEST="$ROOT/metadata/representative_download_manifest.tsv"
mkdir -p "$ROOT/logs" "$ROOT/metadata"

printf 'project\trun\turl\tfile\texpected_md5\tobserved_md5\tsize_bytes\tstatus\ttimestamp\n' > "$MANIFEST"

download_run() {
  local project="$1" run="$2"
  local dest="$ROOT/projects/$project/raw/representative"
  mkdir -p "$dest"
  local api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
  local record
  record=$(curl --fail --silent --show-error --location --retry 4 --connect-timeout 20 "$api" | tail -n 1)
  if [[ -z "$record" || "$record" == run_accession* ]]; then
    echo "ERROR no ENA FASTQ record for $run" | tee -a "$LOG" >&2
    return 1
  fi
  local accession urls md5s sizes
  IFS=$'\t' read -r accession urls md5s sizes <<< "$record"
  IFS=';' read -r -a url_array <<< "$urls"
  IFS=';' read -r -a md5_array <<< "$md5s"
  IFS=';' read -r -a size_array <<< "$sizes"
  for i in "${!url_array[@]}"; do
    local url="https://${url_array[$i]}"
    local file="$dest/$(basename "${url_array[$i]}")"
    local expected="${md5_array[$i]}"
    local expected_size="${size_array[$i]}"
    if [[ -s "$file" ]]; then
      local have
      have=$(md5sum "$file" | awk '{print $1}')
      if [[ "$have" == "$expected" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\treused_verified\t%s\n' "$project" "$run" "$url" "$file" "$expected" "$have" "$(stat -c %s "$file")" "$(date -Is)" >> "$MANIFEST"
        continue
      fi
      echo "ERROR existing file checksum mismatch; refusing overwrite: $file" | tee -a "$LOG" >&2
      return 1
    fi
    echo "Downloading $project $run $(basename "$file") expected_bytes=$expected_size" | tee -a "$LOG"
    wget --continue --tries=6 --timeout=30 --output-document="$file.part" "$url" 2>> "$LOG"
    local observed
    observed=$(md5sum "$file.part" | awk '{print $1}')
    if [[ "$observed" != "$expected" ]]; then
      printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tchecksum_failed\t%s\n' "$project" "$run" "$url" "$file.part" "$expected" "$observed" "$(stat -c %s "$file.part")" "$(date -Is)" >> "$MANIFEST"
      echo "ERROR checksum mismatch for $file.part; retained without rename" | tee -a "$LOG" >&2
      return 1
    fi
    mv "$file.part" "$file"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tdownloaded_verified\t%s\n' "$project" "$run" "$url" "$file" "$expected" "$observed" "$(stat -c %s "$file")" "$(date -Is)" >> "$MANIFEST"
  done
}

download_run PRJNA935967 SRR23502889
download_run PRJNA993383 SRR25224740
download_run PRJNA822526 SRR18577946
download_run PRJNA822526 SRR18578611

echo "Representative downloads complete: $MANIFEST" | tee -a "$LOG"


