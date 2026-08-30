#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 06_dada2_16s_band32_diagnostic.R <project_root> <threads>")
}
project <- normalizePath(args[[1]], mustWork = TRUE)
threads <- as.integer(args[[2]])
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer")

benchmark <- file.path(project, "results", "benchmark_a2_2_full_v2", "16S")
filtered <- file.path(benchmark, "A2-2_16S_filtered.fastq.gz")
error_path <- file.path(benchmark, "error_model_benchmark.rds")
baseline_path <- file.path(benchmark, "seqtab_benchmark.rds")
outdir <- file.path(project, "results", "diagnostic_16s_band32_v1")

for (path in c(filtered, error_path, baseline_path)) test <- stopifnot(file.exists(path))
if (dir.exists(outdir)) stop("Refusing to overwrite existing diagnostic directory: ", outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

err <- readRDS(error_path)
baseline <- readRDS(baseline_path)
derep <- derepFastq(filtered, verbose = TRUE)

started <- Sys.time()
dada_band32 <- dada(
  derep,
  err = err,
  pool = FALSE,
  BAND_SIZE = 32,
  multithread = threads,
  verbose = TRUE
)
finished <- Sys.time()
seqtab_band32 <- makeSequenceTable(setNames(list(dada_band32), "A2-2"))

saveRDS(dada_band32, file.path(outdir, "dada_A2-2_16S_band32.rds"))
saveRDS(seqtab_band32, file.path(outdir, "seqtab_A2-2_16S_band32.rds"))

filtered_reads <- sum(derep$uniques)
summary <- data.frame(
  configuration = c("baseline_default_band", "pacbio_band32"),
  filtered_reads = filtered_reads,
  reads_denoised = c(sum(baseline), sum(seqtab_band32)),
  denoise_retention = c(sum(baseline), sum(seqtab_band32)) / filtered_reads,
  asvs = c(ncol(baseline), ncol(seqtab_band32)),
  elapsed_seconds = c(NA_real_, as.numeric(difftime(finished, started, units = "secs"))),
  threads = c(1L, threads),
  stringsAsFactors = FALSE
)
write.table(summary, file.path(outdir, "band32_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  paste0("error_quality_min\t", min(as.integer(colnames(err)))),
  paste0("error_quality_max\t", max(as.integer(colnames(err)))),
  paste0("filtered_unique_sequences\t", length(derep$uniques)),
  paste0("filtered_reads\t", filtered_reads),
  paste0("threads\t", threads),
  paste0("BAND_SIZE\t32")
), file.path(outdir, "diagnostic_parameters.tsv"))

writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Completed A2-2 16S BAND_SIZE=32 diagnostic")

