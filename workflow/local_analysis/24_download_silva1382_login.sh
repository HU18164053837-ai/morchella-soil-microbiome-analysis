#!/usr/bin/env bash
set -euo pipefail

PROJECT=/path/to/morchella_microbiome
DBDIR="$PROJECT/databases/silva_138_2_dada2"
mkdir -p "$DBDIR"
cd "$DBDIR"

download_one() {
  local name="$1" url="$2" md5="$3"
  if [[ -s "$name" ]] && echo "$md5  $name" | md5sum -c -; then
    return 0
  fi
  curl --fail --location --retry 8 --retry-delay 10 --continue-at - \
    --output "$name" "$url"
  echo "$md5  $name" | md5sum -c -
}

download_one \
  silva_nr99_v138.2_toGenus_trainset.fa.gz \
  'https://zenodo.org/records/14169026/files/silva_nr99_v138.2_toGenus_trainset.fa.gz?download=1' \
  1764e2a36b4500ccb1c7d5261948a414

download_one \
  silva_v138.2_assignSpecies.fa.gz \
  'https://zenodo.org/records/14169026/files/silva_v138.2_assignSpecies.fa.gz?download=1' \
  62fd939bd1aa9832d2089847f9f0b87d

md5sum *.fa.gz > SILVA_138_2.md5
date -Is > DOWNLOAD_COMPLETED.txt


