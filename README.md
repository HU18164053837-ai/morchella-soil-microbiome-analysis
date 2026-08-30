# Morel greenhouse soil microbiome analysis

Reproducible code archive for the bacterial 16S rRNA gene and fungal ITS amplicon analyses supporting a study of greenhouse soils associated with morel cultivation, continuous-cropping history, and a spatially restricted white-mold case comparison.

## Study design represented in the code

- Paired cultivation comparison: soil collected before (`M`) and after (`A`) morel cultivation at matched spatial positions in three greenhouses with 2-, 3-, and 4-year continuous-cropping histories.
- White-mold case comparison: healthy (`H`) and diseased (`D`) points sampled concurrently within one greenhouse, with three spatial subsamples per point. These subsamples are not independent greenhouse replicates.
- Each physical soil sample has separate full-length bacterial 16S and fungal ITS amplicon libraries.

## Repository layout

- `workflow/local_analysis/`: raw-read inventory, primer processing, DADA2, taxonomy, diversity, compositional analyses, soil integration, sensitivity analyses, and reviewer-requested statistics.
- `workflow/hpc_slurm/`: Slurm submission scripts used during pipeline development and formal runs.
- `workflow/functional_and_statistical_extensions/`: statistical extensions, FUNGuild abundance analysis, PICRUSt2 preparation and compositional analysis.
- `workflow/public_database_validation/`: public-read download/QC, DADA2 processing, and contrast-level cross-project evidence synthesis.
- `figures/`: scripts used to construct and academically restyle manuscript figures.
- `manuscript_tools/`: utilities used to assemble or inspect manuscript artifacts; these are not part of the biological analysis pipeline.
- `config/sample_metadata_columns.tsv`: metadata schema only. Sample-specific local paths are deliberately excluded.
- `code_manifest.tsv`: SHA-256 manifest of the public code snapshot.

## Important reproducibility notes

This repository is a sanitized code snapshot. Machine-specific paths were replaced with:

- Windows project root: `D:\path\to\morchella_microbiome`
- HPC project root: `/path/to/morchella_microbiome`
- Local auxiliary workspace: `C:\path\to\workspace`

Set these paths for your environment before running a script. Slurm partitions, container locations, memory, CPUs, database paths, and module names must also be adapted to the target cluster.

The repository intentionally does not contain raw sequences, private metadata, large derived objects, database files, third-party FUNGuild source code, manuscripts, or unpublished result tables. Raw reads will be linked through NCBI accessions after public release.

## Main computational stages

1. Raw-file inventory, checksum/QC, primer orientation and trimming.
2. Marker-specific PacBio CCS denoising with DADA2.
3. ITS taxonomy with UNITE and bacterial taxonomy with SILVA.
4. Alpha/beta diversity, paired design-aware analyses, compositional testing, and soil-variable integration.
5. Candidate-genus screening and robustness/sensitivity analyses.
6. FUNGuild trophic-mode analysis and PICRUSt2 functional prediction.
7. Public-dataset contrast synthesis using explicitly documented post hoc qualification rules.
8. Manuscript figure generation.

## Software

Core software includes R, DADA2, phyloseq, vegan, ggplot2, Python 3, Cutadapt, BLAST+, UNITE, SILVA, FUNGuild and PICRUSt2. Exact package snapshots from individual analyses are recorded in the associated project archive; cluster environments should be rebuilt using the versions described in the manuscript Methods and public release notes.

## Data availability

NCBI BioProject: `TO_BE_ADDED`

GitHub release/DOI: `TO_BE_ADDED`

## Citation

Please cite the associated manuscript when it becomes available. Candidate directions and external-evidence grades are study-specific and should not be interpreted as validated biomarkers or causal effects.
