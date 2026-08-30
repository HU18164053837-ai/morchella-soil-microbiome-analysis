#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 08_dada2_16s_balanced_A_resume_v2.R <project_root> <threads>")
}
project <- normalizePath(args[[1]], mustWork = TRUE)
threads <- as.integer(args[[2]])
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer")

sample_ids <- c("A2-1", "A3-1", "A4-1")
source_dir <- file.path(project, "results", "diagnostic_16s_balanced_A_v1")
outdir <- file.path(project, "results", "diagnostic_16s_balanced_A_v2")
filtered <- file.path(source_dir, "filtered", paste0(sample_ids, "_16S_filtered.fastq.gz"))
names(filtered) <- sample_ids
error_path <- file.path(source_dir, "error_model_A_balanced.rds")
independent_dada_path <- file.path(source_dir, "dada_A_independent.rds")
independent_seqtab_path <- file.path(source_dir, "seqtab_A_independent.rds")
filter_stats_path <- file.path(source_dir, "filter_stats.tsv")

required <- c(filtered, error_path, independent_dada_path,
              independent_seqtab_path, filter_stats_path)
stopifnot(all(file.exists(required)))
if (dir.exists(outdir)) stop("Refusing to overwrite existing v2 directory: ", outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

filter_stats <- read.delim(filter_stats_path, check.names = FALSE)
filter_stats$sample_id <- sub("_16S_trimmed.fastq.gz$", "", filter_stats$sample_id)
if (!setequal(filter_stats$sample_id, sample_ids)) {
  stop("Unexpected sample IDs in filter_stats.tsv")
}
filtered_counts <- setNames(filter_stats$reads.out, filter_stats$sample_id)

err <- readRDS(error_path)
independent <- readRDS(independent_dada_path)
if (length(independent) != length(sample_ids)) stop("Unexpected independent result length")
names(independent) <- sample_ids
independent_seqtab <- readRDS(independent_seqtab_path)

summarize_dada <- function(ans, label) {
  counts <- vapply(ans, function(x) sum(getUniques(x)), numeric(1))
  data.frame(
    sample_id = sample_ids,
    mode = label,
    reads_filtered = as.numeric(filtered_counts[sample_ids]),
    reads_denoised = as.numeric(counts[sample_ids]),
    denoise_retention = as.numeric(counts[sample_ids]) /
      as.numeric(filtered_counts[sample_ids]),
    stringsAsFactors = FALSE
  )
}

independent_summary <- summarize_dada(independent, "independent")
saveRDS(independent_seqtab, file.path(outdir, "seqtab_A_independent_reused.rds"))

pseudo_started <- Sys.time()
pseudo <- dada(
  filtered,
  err = err,
  pool = "pseudo",
  BAND_SIZE = 32,
  multithread = threads,
  verbose = TRUE
)
pseudo_finished <- Sys.time()
names(pseudo) <- sample_ids
pseudo_seqtab <- makeSequenceTable(pseudo)
saveRDS(pseudo, file.path(outdir, "dada_A_pseudo.rds"))
saveRDS(pseudo_seqtab, file.path(outdir, "seqtab_A_pseudo.rds"))

comparison <- rbind(independent_summary, summarize_dada(pseudo, "pseudo"))
write.table(comparison, file.path(outdir, "inference_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(filter_stats, file.path(outdir, "filter_stats_reused.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(err, file.path(outdir, "error_model_A_balanced_reused.rds"))

writeLines(c(
  "marker\t16S",
  "batch\tA",
  "source_job\t119575319",
  "source_result\tdiagnostic_16s_balanced_A_v1",
  "error_model\treused",
  "independent_inference\treused",
  "pseudo_pooling\tnew",
  "BAND_SIZE\t32",
  paste0("threads\t", threads),
  paste0("pseudo_seconds\t", as.numeric(difftime(pseudo_finished, pseudo_started, units = "secs")))
), file.path(outdir, "diagnostic_parameters.tsv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Completed balanced A-batch 16S v2 resume diagnostic")

