#!/usr/bin/env bash
#SBATCH --job-name=mor_pic_coda
#SBATCH --partition=kshctest02
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --output=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_pic_coda_%j.out
#SBATCH --error=/path/to/morchella_microbiome/extension_analysis_20260822/logs/mor_pic_coda_%j.err
set -euo pipefail
BASE=/path/to/morchella_microbiome
EXT=$BASE/extension_analysis_20260822
MAMBA=$EXT/software/bin/micromamba
MAMBA_ROOT_PREFIX=$EXT/software/mamba_root_offline
ENV=$EXT/software/envs/picrust2_coda_offline_v2
OFFLINE_TAR=$EXT/software/picrust2_offline_packages.tar
OFFLINE_DIR=$EXT/software/picrust2_offline_bundle/packages
EXPLICIT=$EXT/software/picrust2_offline_explicit_v2.txt
PATCH_DIR=$EXT/software/picrust2_coda_patch_v1
OUT=$EXT/results/picrust2_coda_v1
mkdir -p "$MAMBA_ROOT_PREFIX" "$EXT/software/envs" "$OUT"
export MAMBA_ROOT_PREFIX

if [ "$(find "$OFFLINE_DIR" -maxdepth 1 -type f 2>/dev/null | wc -l)" -ne 292 ]; then
  mkdir -p "$EXT/software/picrust2_offline_bundle"
  tar -xf "$OFFLINE_TAR" -C "$EXT/software/picrust2_offline_bundle"
fi

{
  echo '@EXPLICIT'
  find "$OFFLINE_DIR" -maxdepth 1 -type f \( -name '*.conda' -o -name '*.tar.bz2' \) \
    ! -name 'r-cvxr-*' ! -name 'r-lme4-*' ! -name 'r-matrix-*' \
    | sort | sed 's#^#file://#'
  find "$PATCH_DIR" -maxdepth 1 -type f \( -name '*.conda' -o -name '*.tar.bz2' \) \
    | sort | sed 's#^#file://#'
} > "$EXPLICIT"
test "$(grep -c '^file://' "$EXPLICIT")" -eq 296

if [ ! -x "$ENV/bin/Rscript" ]; then
  "$MAMBA" create -y -p "$ENV" --offline --file "$EXPLICIT"
fi

if ! "$MAMBA" run -p "$ENV" Rscript -e "quit(status=ifelse(requireNamespace('GenomeInfoDbData',quietly=TRUE),0,1))"; then
  module load apps/apptainer/1.3.4
  apptainer exec --cleanenv --bind "$BASE:$BASE" \
    "$BASE/containers/bioconductor-dada2_1.34.0--r44he5774e6_0.sif" \
    cp -a /usr/local/lib/R/library/GenomeInfoDbData "$ENV/lib/R/library/"
fi

"$MAMBA" run -p "$ENV" Rscript -e "stopifnot(as.character(packageVersion('Matrix'))=='1.7.4',as.character(packageVersion('lme4'))=='1.1.38',as.character(packageVersion('CVXR'))=='1.0.15',requireNamespace('ALDEx2',quietly=TRUE),requireNamespace('data.table',quietly=TRUE)); cat('LOCKED_OFFLINE_R_ENV_V2_OK\\n')"

"$MAMBA" run -p "$ENV" Rscript \
  "$EXT/scripts/06_picrust2_coda.R" \
  "$EXT/results/picrust2_v4/pathways_out/path_abun_unstrat.tsv.gz" \
  "$EXT/results/picrust2_v4/KO_metagenome_out/weighted_nsti.tsv.gz" \
  "$OUT"

test -s "$OUT/MetaCyc_consensus_and_NSTI_audit.tsv"
echo MORCHELLA_PICRUST2_CODA_JOB_COMPLETE

