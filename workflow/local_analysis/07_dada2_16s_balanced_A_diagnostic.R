#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 07_dada2_16s_balanced_A_diagnostic.R <project_root> <threads>")
}
project <- normalizePath(args[[1]], mustWork = TRUE)
threads <- as.integer(args[[2]])
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer")

sample_ids <- c("A2-1", "A3-1", "A4-1")
reads_per_sample <- 10000L
seed <- 20260820L
outdir <- file.path(project, "results", "diagnostic_16s_balanced_A_v1")
filtered_dir <- file.path(outdir, "filtered")
training_dir <- file.path(outdir, "error_training")

if (dir.exists(outdir)) stop("Refusing to overwrite existing diagnostic directory: ", outdir)
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(training_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- file.path(project, "input", "trimmed", "16S",
                    paste0(sample_ids, "_16S_trimmed.fastq.gz"))
filtered <- file.path(filtered_dir, paste0(sample_ids, "_16S_filtered.fastq.gz"))
training <- file.path(training_dir, paste0(sample_ids, "_16S_errortrain_10000.fastq.gz"))
names(inputs) <- names(filtered) <- names(training) <- sample_ids
stopifnot(all(file.exists(inputs)))

filter_stats <- filterAndTrim(
  inputs, filtered,
  maxN = 0, maxEE = 2, truncQ = 2,
  minLen = 1200, maxLen = 1800,
  rm.phix = FALSE, compress = TRUE,
  multithread = threads, verbose = TRUE
)
write.table(
  data.frame(sample_id = rownames(filter_stats), filter_stats, check.names = FALSE),
  file.path(outdir, "filter_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

set.seed(seed)
for (id in sample_ids) {
  fq <- readFastq(filtered[[id]])
  if (length(fq) < reads_per_sample) stop("Too few filtered reads for ", id)
  idx <- sort(sample.int(length(fq), reads_per_sample, replace = FALSE))
  writeFastq(fq[idx], training[[id]], compress = TRUE)
  rm(fq)
  gc()
}

learn_started <- Sys.time()
err <- learnErrors(
  unname(training),
  nbases = 1e8,
  errorEstimationFunction = PacBioErrfun,
  BAND_SIZE = 32,
  multithread = threads,
  randomize = TRUE,
  verbose = TRUE
)
learn_finished <- Sys.time()
saveRDS(err, file.path(outdir, "error_model_A_balanced.rds"))
pdf(file.path(outdir, "error_model_A_balanced.pdf"), width = 8, height = 6)
print(plotErrors(err, nominalQ = TRUE))
dev.off()

run_dada <- function(pool_mode, label) {
  started <- Sys.time()
  ans <- dada(
    unname(filtered), err = err, pool = pool_mode,
    BAND_SIZE = 32, multithread = threads, verbose = TRUE
  )
  names(ans) <- sample_ids
  seqtab <- makeSequenceTable(ans)
  saveRDS(ans, file.path(outdir, paste0("dada_A_", label, ".rds")))
  saveRDS(seqtab, file.path(outdir, paste0("seqtab_A_", label, ".rds")))
  counts <- vapply(ans, function(x) sum(getUniques(x)), numeric(1))
  list(
    summary = data.frame(
      sample_id = sample_ids,
      mode = label,
      reads_filtered = as.numeric(filter_stats[sample_ids, "reads.out"]),
      reads_denoised = as.numeric(counts[sample_ids]),
      denoise_retention = as.numeric(counts[sample_ids]) /
        as.numeric(filter_stats[sample_ids, "reads.out"]),
      stringsAsFactors = FALSE
    ),
    elapsed = as.numeric(difftime(Sys.time(), started, units = "secs"))
  )
}

independent <- run_dada(FALSE, "independent")
pseudo <- run_dada("pseudo", "pseudo")
comparison <- rbind(independent$summary, pseudo$summary)
write.table(comparison, file.path(outdir, "inference_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  paste0("marker\t16S"),
  paste0("batch\tA"),
  paste0("training_samples\t", paste(sample_ids, collapse = ",")),
  paste0("training_reads_per_sample\t", reads_per_sample),
  paste0("training_reads_total\t", reads_per_sample * length(sample_ids)),
  paste0("random_seed\t", seed),
  paste0("BAND_SIZE\t32"),
  paste0("threads\t", threads),
  paste0("learn_errors_seconds\t", as.numeric(difftime(learn_finished, learn_started, units = "secs"))),
  paste0("independent_seconds\t", independent$elapsed),
  paste0("pseudo_seconds\t", pseudo$elapsed)
), file.path(outdir, "diagnostic_parameters.tsv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Completed balanced A-batch 16S diagnostic")

