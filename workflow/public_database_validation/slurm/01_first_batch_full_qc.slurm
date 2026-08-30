#!/usr/bin/env bash
#SBATCH --job-name=morchella_pv_qc
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/path/to/morchella_microbiome/public_validation/logs/%x_%A_%a.out
#SBATCH --error=/path/to/morchella_microbiome/public_validation/logs/%x_%A_%a.err

set -euo pipefail
ROOT=/path/to/morchella_microbiome/public_validation
PROJECTS=(PRJNA935967 PRJNA993383)
PROJECT=${PROJECTS[$SLURM_ARRAY_TASK_ID]}

source /etc/profile.d/modules.sh
module purge
module load python/3.8.10
cd "$ROOT/projects/$PROJECT"

find raw/fastq -maxdepth 1 -type f -name '*.fastq.gz' -print0 | sort -z | xargs -0 -n1 gzip -t
python3 "$ROOT/scripts/05_full_fastq_qc.py" --project "$PROJECT"
sha256sum raw/fastq/*.fastq.gz > qc/full_fastq.sha256
echo "QC_COMPLETE project=$PROJECT files=$(find raw/fastq -maxdepth 1 -type f -name '*.fastq.gz' | wc -l)"

