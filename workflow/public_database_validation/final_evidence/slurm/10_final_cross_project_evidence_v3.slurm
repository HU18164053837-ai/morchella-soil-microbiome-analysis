#!/usr/bin/env bash
#SBATCH --job-name=mor_final_evid_v3
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=3G
#SBATCH --time=01:00:00
#SBATCH --output=/path/to/morchella_microbiome/public_validation/logs/mor_final_evid_v3_%j.out
#SBATCH --error=/path/to/morchella_microbiome/public_validation/logs/mor_final_evid_v3_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome
PV=$BASE/public_validation
IMG=$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif
SCRIPT=$PV/scripts/15_final_cross_project_evidence_v3.R
FIRST=$PV/results/first_batch_evidence_v2/first_batch_effect_sizes.tsv
SECOND=$PV/projects/PRJNA822526/tables/external_validation_v1/core_genera_external_validation.tsv
OUT=$PV/results/final_cross_project_evidence_v3
test ! -e "$OUT"; test -s "$SCRIPT"; test -s "$FIRST"; test -s "$SECOND"
module load apps/apptainer/1.3.4
apptainer exec --cleanenv --bind "$BASE:$BASE" "$IMG" Rscript "$SCRIPT" "$FIRST" "$SECOND" "$OUT"
test "$(awk 'END{print NR-1}' "$OUT/final_cross_project_evidence_matrix.tsv")" -eq 8
test -s "$OUT/Figure_final_cross_project_evidence.pdf"
test -s "$OUT/Figure_final_cross_project_evidence.svg"
test -s "$OUT/Figure_final_cross_project_evidence.tiff"
test -s "$OUT/羊肚菌公共数据库外部验证_最终证据报告_v3.md"
test -s "$OUT/sessionInfo.txt"
echo MORCHELLA_FINAL_CROSS_PROJECT_EVIDENCE_V2_COMPLETE



