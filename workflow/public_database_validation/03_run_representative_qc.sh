#!/usr/bin/env bash
set -euo pipefail
ROOT=/path/to/morchella_microbiome/public_validation
for f in "$ROOT"/projects/*/raw/representative/*.fastq.gz; do
  gzip -t "$f"
done
if [[ -x /usr/bin/python3 ]] && /usr/bin/python3 --version >/dev/null 2>&1; then
  PYTHON=/usr/bin/python3
else
  source /etc/profile.d/modules.sh
  module load python/3.8.10
  PYTHON=$(command -v python3)
fi
"$PYTHON" "$ROOT/scripts/02_inspect_representative_fastq.py"
sha256sum "$ROOT"/projects/*/raw/representative/*.fastq.gz > "$ROOT/qc/representative_fastq.sha256"
echo "QC complete: $ROOT/qc/representative_structure_qc.tsv"

