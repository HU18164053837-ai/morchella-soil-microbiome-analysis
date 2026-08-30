#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 09_dada2_16s_A3_A4_quality_diagnostic.R <project_root> <threads>")
}

project <- normalizePath(args[[1]], mustWork = TRUE)
threads <- as.integer(args[[2]])
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer")

sample_ids <- c("A3-1", "A4-1")
trimmed_dir <- file.path(project, "input", "trimmed", "16S")
filtered_dir <- file.path(project, "results", "diagnostic_16s_balanced_A_v1", "filtered")
comparison_path <- file.path(project, "results", "diagnostic_16s_balanced_A_v2", "inference_comparison.tsv")
outdir <- file.path(project, "results", "diagnostic_16s_A3_A4_quality_v1")

if (dir.exists(outdir)) stop("Refusing to overwrite existing output: ", outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

fastq_metrics <- function(path, sample_id, stage) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  rows <- list()
  block_id <- 0L
  repeat {
    x <- readLines(con, n = 40000L, warn = FALSE)
    if (!length(x)) break
    if (length(x) %% 4L != 0L) stop("Malformed FASTQ: ", path)
    idx <- seq.int(1L, length(x), by = 4L)
    seqs <- x[idx + 1L]
    quals <- x[idx + 3L]
    qints <- lapply(quals, function(q) utf8ToInt(q) - 33L)
    block_id <- block_id + 1L
    rows[[block_id]] <- data.frame(
      sample_id = sample_id,
      stage = stage,
      length = nchar(seqs, type = "bytes"),
      n_count = lengths(regmatches(seqs, gregexpr("N", seqs, fixed = TRUE))),
      mean_q = vapply(qints, mean, numeric(1)),
      min_q = vapply(qints, min, numeric(1)),
      expected_errors = vapply(qints, function(q) sum(10^(-q / 10)), numeric(1)),
      sequence = seqs,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

process_sample <- function(sample_id) {
  trimmed <- file.path(trimmed_dir, paste0(sample_id, "_16S_trimmed.fastq.gz"))
  filtered <- file.path(filtered_dir, paste0(sample_id, "_16S_filtered.fastq.gz"))
  if (!file.exists(trimmed) || !file.exists(filtered)) stop("Missing FASTQ for ", sample_id)
  rbind(
    fastq_metrics(trimmed, sample_id, "trimmed"),
    fastq_metrics(filtered, sample_id, "filtered")
  )
}

metrics_list <- parallel::mclapply(sample_ids, process_sample,
                                   mc.cores = min(threads, length(sample_ids)))
metrics <- do.call(rbind, metrics_list)

summarize_group <- function(d) {
  abundance <- sort(table(d$sequence), decreasing = TRUE)
  data.frame(
    sample_id = d$sample_id[[1]],
    stage = d$stage[[1]],
    reads = nrow(d),
    unique_sequences = length(abundance),
    singleton_fraction = mean(abundance == 1L),
    top_sequence_fraction = as.numeric(abundance[[1]]) / nrow(d),
    length_median = median(d$length),
    length_q05 = unname(quantile(d$length, 0.05)),
    length_q95 = unname(quantile(d$length, 0.95)),
    mean_q_median = median(d$mean_q),
    expected_errors_median = median(d$expected_errors),
    expected_errors_q95 = unname(quantile(d$expected_errors, 0.95)),
    pct_with_N = mean(d$n_count > 0) * 100,
    pct_EE_gt_2 = mean(d$expected_errors > 2) * 100,
    pct_length_outside_1200_1800 = mean(d$length < 1200 | d$length > 1800) * 100,
    stringsAsFactors = FALSE
  )
}

groups <- split(metrics, interaction(metrics$sample_id, metrics$stage, drop = TRUE))
summary_table <- do.call(rbind, lapply(groups, summarize_group))
summary_table <- summary_table[order(summary_table$sample_id, summary_table$stage), ]
write.table(summary_table, file.path(outdir, "quality_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(metrics[, setdiff(names(metrics), "sequence")],
            file.path(outdir, "per_read_quality_metrics.tsv.gz"),
            sep = "\t", quote = FALSE, row.names = FALSE)

comparison <- read.delim(comparison_path, check.names = FALSE)
comparison <- comparison[comparison$sample_id %in% sample_ids, ]
write.table(comparison, file.path(outdir, "dada_inference_comparison_reused.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(outdir, "quality_diagnostic_plots.pdf"), width = 10, height = 8)
par(mfrow = c(2, 2), mar = c(4, 4, 3, 1))
for (sample_id in sample_ids) {
  d <- metrics[metrics$sample_id == sample_id & metrics$stage == "trimmed", ]
  hist(d$expected_errors, breaks = 80, xlim = c(0, unname(quantile(d$expected_errors, 0.99))),
       main = paste(sample_id, "trimmed expected errors"), xlab = "Expected errors")
  abline(v = 2, col = "red", lty = 2)
  hist(d$length, breaks = 80, main = paste(sample_id, "trimmed lengths"), xlab = "Length (bp)")
  abline(v = c(1200, 1800), col = "red", lty = 2)
}
dev.off()

writeLines(c(
  "marker\t16S",
  "samples\tA3-1,A4-1",
  "analysis\tread-level quality and abundance diagnostic",
  "filter_thresholds\tmaxN=0,maxEE=2,truncQ=2,minLen=1200,maxLen=1800",
  paste0("threads\t", threads)
), file.path(outdir, "diagnostic_parameters.tsv"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
message("Completed A3-1/A4-1 16S quality diagnostic")

