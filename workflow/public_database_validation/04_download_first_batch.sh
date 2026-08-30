#!/usr/bin/env bash
set -euo pipefail

ROOT=/path/to/morchella_microbiome/public_validation
META="$ROOT/metadata/sample_metadata_validated.tsv"
STATE="$ROOT/metadata/first_batch_download_state"
MANIFEST="$ROOT/metadata/first_batch_download_manifest.tsv"
LOG="$ROOT/logs/first_batch_download_$(date +%Y%m%d_%H%M%S).log"
mkdir -p "$STATE" "$ROOT/logs"

download_one() {
  local project="$1" run="$2"
  local dest="$ROOT/projects/$project/raw/fastq"
  local rep="$ROOT/projects/$project/raw/representative"
  mkdir -p "$dest"
  local api="https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${run}&result=read_run&fields=run_accession,fastq_ftp,fastq_md5,fastq_bytes&format=tsv"
  local record accession urls md5s sizes
  record=$(curl --fail --silent --show-error --location --retry 5 --connect-timeout 20 "$api" | tail -n 1)
  IFS=$'\t' read -r accession urls md5s sizes <<< "$record"
  if [[ "$accession" != "$run" || -z "$urls" ]]; then
    printf '%s\t%s\tNA\tNA\tNA\tNA\t0\tmetadata_failed\t%s\n' "$project" "$run" "$(date -Is)" > "$STATE/$run.tsv"
    return 1
  fi
  IFS=';' read -r -a ua <<< "$urls"; IFS=';' read -r -a ma <<< "$md5s"; IFS=';' read -r -a sa <<< "$sizes"
  : > "$STATE/$run.tsv.tmp"
  for i in "${!ua[@]}"; do
    local url="https://${ua[$i]}" base file expected observed status
    base=$(basename "${ua[$i]}"); file="$dest/$base"; expected="${ma[$i]}"
    if [[ -s "$file" ]]; then
      observed=$(md5sum "$file" | awk '{print $1}')
      if [[ "$observed" != "$expected" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\texisting_checksum_failed\t%s\n' "$project" "$run" "$url" "$file" "$expected" "$observed" "$(stat -c %s "$file")" "$(date -Is)" >> "$STATE/$run.tsv.tmp"
        return 1
      fi
      status=reused_verified
    elif [[ -s "$rep/$base" ]] && [[ "$(md5sum "$rep/$base" | awk '{print $1}')" == "$expected" ]]; then
      ln "$rep/$base" "$file"
      observed="$expected"; status=hardlinked_representative_verified
    else
      local tmp="$file.v2.part"
      if ! curl --fail --location --retry 8 --connect-timeout 30 --continue-at - --output "$tmp" "$url" >/dev/null 2>&1; then
        printf '%s\t%s\t%s\t%s\t%s\tNA\t%s\tdownload_failed_retryable\t%s\n' "$project" "$run" "$url" "$tmp" "$expected" "$(stat -c %s "$tmp" 2>/dev/null || echo 0)" "$(date -Is)" >> "$STATE/$run.tsv.tmp"
        return 1
      fi
      observed=$(md5sum "$tmp" | awk '{print $1}')
      if [[ "$observed" != "$expected" ]]; then
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tchecksum_failed\t%s\n' "$project" "$run" "$url" "$tmp" "$expected" "$observed" "$(stat -c %s "$tmp")" "$(date -Is)" >> "$STATE/$run.tsv.tmp"
        return 1
      fi
      mv "$tmp" "$file"; status=downloaded_verified
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$project" "$run" "$url" "$file" "$expected" "$observed" "$(stat -c %s "$file")" "$status" "$(date -Is)" >> "$STATE/$run.tsv.tmp"
  done
  mv "$STATE/$run.tsv.tmp" "$STATE/$run.tsv"
}
export -f download_one
export ROOT STATE

awk -F '\t' 'NR>1 && ($1=="PRJNA935967" || $1=="PRJNA993383") {print $1"\t"$2}' "$META" | sort -u > "$STATE/queue_all.tsv"
: > "$STATE/queue.tsv"
while IFS=$'\t' read -r project run; do
  if [[ -s "$STATE/$run.tsv" ]] && awk -F '\t' '$8 !~ /verified$/ {bad=1} END {exit bad}' "$STATE/$run.tsv"; then
    continue
  fi
  printf '%s\t%s\n' "$project" "$run" >> "$STATE/queue.tsv"
done < "$STATE/queue_all.tsv"
wave_size=${WAVE_SIZE:-16}
head -n "$wave_size" "$STATE/queue.tsv" > "$STATE/queue_wave.tsv"
concurrency=${CONCURRENCY:-4}
echo "Starting/resuming first-batch download: total=60 remaining=$(wc -l < "$STATE/queue.tsv") wave=$(wc -l < "$STATE/queue_wave.tsv") concurrent=$concurrency" | tee -a "$LOG"
if [[ -s "$STATE/queue_wave.tsv" ]]; then
  xargs -P "$concurrency" -n 2 bash -c 'download_one "$0" "$1"' < "$STATE/queue_wave.tsv" || true
fi

printf 'project\trun\turl\tfile\texpected_md5\tobserved_md5\tsize_bytes\tstatus\ttimestamp\n' > "$MANIFEST"
for f in "$STATE"/SRR*.tsv; do [[ -s "$f" ]] && cat "$f"; done | sort -k1,1 -k2,2 >> "$MANIFEST"

expected=$(wc -l < "$STATE/queue_all.tsv")
completed=$(awk -F '\t' 'NR>1 && $8 ~ /verified$/ {a[$2]=1} END {print length(a)}' "$MANIFEST")
failed=$(awk -F '\t' 'NR>1 && $8 !~ /verified$/ {a[$2]=1} END {print length(a)}' "$MANIFEST")
bytes=$(awk -F '\t' 'NR>1 && $8 ~ /verified$/ {s+=$7} END {printf "%.0f",s}' "$MANIFEST")
printf 'expected_runs\tcompleted_runs\tfailed_runs\tverified_bytes\tfinished_at\n%s\t%s\t%s\t%s\t%s\n' "$expected" "$completed" "$failed" "$bytes" "$(date -Is)" > "$ROOT/metadata/first_batch_download_summary.tsv"
cat "$ROOT/metadata/first_batch_download_summary.tsv" | tee -a "$LOG"
if [[ "$completed" -eq "$expected" && "$failed" -eq 0 ]]; then
  echo COMPLETE
else
  echo INCOMPLETE_RESUME_REQUIRED
fi

