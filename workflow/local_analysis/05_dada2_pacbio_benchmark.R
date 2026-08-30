#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) {
  stop("Usage: Rscript 05_dada2_pacbio_benchmark.R <project_root>")
}
project <- normalizePath(args[[1]], mustWork = TRUE)
outdir <- file.path(project, "results", "benchmark_a2_2")
if (dir.exists(outdir)) stop("Refusing to overwrite existing benchmark directory: ", outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

threads <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1"))
if (is.na(threads) || threads < 1L) threads <- 1L

run_marker <- function(marker, min_len, max_len) {
  input <- file.path(
    project, "input", "trimmed", marker,
    paste0("A2-2_", marker, "_trimmed.fastq.gz")
  )
  if (!file.exists(input)) stop("Missing benchmark input: ", input)
  marker_out <- file.path(outdir, marker)
  dir.create(marker_out, recursive = TRUE, showWarnings = FALSE)
  filtered <- file.path(marker_out, paste0("A2-2_", marker, "_filtered.fastq.gz"))

  started <- Sys.time()
  filt <- filterAndTrim(
    input, filtered,
    maxN = 0, maxEE = 2, truncQ = 2,
    minLen = min_len, maxLen = max_len,
    rm.phix = FALSE, compress = TRUE,
    multithread = threads, verbose = TRUE
  )
  filtered_at <- Sys.time()
  err <- learnErrors(
    filtered,
    nbases = 2e7,
    errorEstimationFunction = PacBioErrfun,
    multithread = threads,
    randomize = TRUE,
    verbose = TRUE
  )
  learned_at <- Sys.time()
  derep <- derepFastq(filtered, verbose = TRUE)
  dada_out <- dada(derep, err = err, pool = FALSE,
                   multithread = threads, verbose = TRUE)
  seqtab <- makeSequenceTable(setNames(list(dada_out), "A2-2"))
  finished <- Sys.time()

  saveRDS(err, file.path(marker_out, "error_model_benchmark.rds"))
  saveRDS(seqtab, file.path(marker_out, "seqtab_benchmark.rds"))
  pdf(file.path(marker_out, "error_model_benchmark.pdf"), width = 8, height = 6)
  print(plotErrors(err, nominalQ = TRUE))
  dev.off()

  data.frame(
    marker = marker,
    sample_id = "A2-2",
    reads_in = unname(filt[1, "reads.in"]),
    reads_filtered = unname(filt[1, "reads.out"]),
    asvs = ncol(seqtab),
    reads_denoised = sum(seqtab),
    filter_seconds = as.numeric(difftime(filtered_at, started, units = "secs")),
    learn_errors_seconds = as.numeric(difftime(learned_at, filtered_at, units = "secs")),
    denoise_seconds = as.numeric(difftime(finished, learned_at, units = "secs")),
    total_seconds = as.numeric(difftime(finished, started, units = "secs")),
    threads = threads,
    stringsAsFactors = FALSE
  )
}

set.seed(20260819)
summary <- rbind(
  run_marker("16S", 1200, 1800),
  run_marker("ITS", 400, 1500)
)
write.table(summary, file.path(outdir, "benchmark_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
writeLines(c(
  paste0("R\t", R.version.string),
  paste0("dada2\t", as.character(packageVersion("dada2")))
), file.path(outdir, "pipeline_versions.txt"))


