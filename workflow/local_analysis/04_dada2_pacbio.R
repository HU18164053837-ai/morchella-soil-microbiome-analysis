#!/usr/bin/env Rscript

.libPaths(c(
  "C:/path/to/workspace/Rlib",
  .libPaths()
))
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[[1]] %in% c("16S", "ITS")) {
  stop("Usage: Rscript 04_dada2_pacbio.R <16S|ITS>")
}
marker <- args[[1]]

root <- "C:/path/to/workspace"
input_dir <- file.path(root, "dada2_runtime", "input", marker)
runtime <- file.path(root, "dada2_runtime", marker)
filtered_dir <- file.path(runtime, "filtered")
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)

if (marker == "16S") {
  min_len <- 1200
  max_len <- 1800
} else {
  min_len <- 400
  max_len <- 1500
}

files <- sort(list.files(input_dir, pattern = "_trimmed.fastq.gz$", full.names = TRUE))
if (length(files) != 24) stop("Expected 24 trimmed FASTQ files for ", marker)
sample_ids <- sub(paste0("_", marker, "_trimmed.fastq.gz$"), "", basename(files))
filtered <- file.path(filtered_dir, paste0(sample_ids, "_", marker, "_filtered.fastq.gz"))
names(files) <- sample_ids
names(filtered) <- sample_ids

filter_stats_path <- file.path(runtime, "filter_stats.tsv")
if (file.exists(filter_stats_path)) {
  missing_filtered <- !file.exists(filtered)
  if (any(missing_filtered)) {
    message("Regenerating ", sum(missing_filtered), " missing filtered FASTQ file(s) for ", marker)
    filterAndTrim(
      files[missing_filtered], filtered[missing_filtered], maxN = 0, maxEE = 2,
      truncQ = 2, minLen = min_len, maxLen = max_len, rm.phix = FALSE,
      compress = TRUE, multithread = FALSE, verbose = TRUE
    )
  } else {
    message("Reusing completed filtered FASTQ files for ", marker)
  }
  saved_filter <- read.delim(filter_stats_path, check.names = FALSE)
  filter_stats <- as.matrix(saved_filter[, c("reads.in", "reads.out")])
  rownames(filter_stats) <- saved_filter$sample_id
} else {
  message("Filtering ", marker, " reads with minLen=", min_len, ", maxLen=", max_len)
  filter_stats <- filterAndTrim(
    files, filtered,
    maxN = 0,
    maxEE = 2,
    truncQ = 2,
    minLen = min_len,
    maxLen = max_len,
    rm.phix = FALSE,
    compress = TRUE,
    multithread = FALSE,
    verbose = TRUE
  )
  write.table(
    data.frame(sample_id = rownames(filter_stats), filter_stats, check.names = FALSE),
    filter_stats_path, sep = "\t", quote = FALSE, row.names = FALSE
  )
}

batch_for_sample <- function(x) {
  if (grepl("^M", x)) return("M")
  if (grepl("^A", x)) return("A")
  return("DH")
}
batches <- vapply(sample_ids, batch_for_sample, character(1))

seqtabs <- list()
denoised_counts <- setNames(numeric(length(sample_ids)), sample_ids)
error_models <- list()

for (batch in c("M", "A", "DH")) {
  ids <- sample_ids[batches == batch]
  batch_files <- filtered[ids]
  training_files <- sort(list.files(
    file.path(runtime, "error_training", batch), pattern = "fastq.gz$", full.names = TRUE
  ))
  if (length(training_files) < 3) {
    stop("Expected at least three balanced error-training subsets for ", marker, " batch ", batch)
  }
  error_path <- file.path(runtime, paste0("error_model_", batch, ".rds"))
  if (file.exists(error_path)) {
    message("Reusing PacBio error model for ", marker, " batch ", batch)
    err <- readRDS(error_path)
  } else {
    message("Learning PacBio error model for ", marker, " batch ", batch, " (", length(ids), " samples)")
    err <- learnErrors(
      training_files,
      nbases = 1e8,
      errorEstimationFunction = PacBioErrfun,
      multithread = FALSE,
      randomize = TRUE,
      verbose = TRUE
    )
  }
  error_models[[batch]] <- err
  saveRDS(err, error_path)
  pdf(file.path(runtime, paste0("error_model_", batch, ".pdf")), width = 8, height = 6)
  print(plotErrors(err, nominalQ = TRUE))
  dev.off()

  derep <- derepFastq(batch_files, verbose = TRUE)
  names(derep) <- ids
  dada_out <- dada(derep, err = err, pool = FALSE, multithread = TRUE, verbose = TRUE)
  seqtab <- makeSequenceTable(dada_out)
  seqtabs[[batch]] <- seqtab
  denoised_counts[rownames(seqtab)] <- rowSums(seqtab)
  saveRDS(dada_out, file.path(runtime, paste0("dada_", batch, ".rds")))
  saveRDS(seqtab, file.path(runtime, paste0("seqtab_", batch, ".rds")))
}

seqtab_all <- Reduce(function(x, y) mergeSequenceTables(x, y), seqtabs)
seqtab_nochim <- removeBimeraDenovo(
  seqtab_all, method = "consensus", multithread = TRUE, verbose = TRUE
)
saveRDS(seqtab_all, file.path(runtime, "seqtab_all.rds"))
saveRDS(seqtab_nochim, file.path(runtime, "seqtab_nochim.rds"))

input_counts <- setNames(filter_stats[, "reads.in"], rownames(filter_stats))
filtered_counts <- setNames(filter_stats[, "reads.out"], rownames(filter_stats))
nonchim_counts <- setNames(rep(0, length(sample_ids)), sample_ids)
nonchim_counts[rownames(seqtab_nochim)] <- rowSums(seqtab_nochim)
track <- data.frame(
  sample_id = sample_ids,
  marker = marker,
  batch = batches,
  input = as.numeric(input_counts[sample_ids]),
  filtered = as.numeric(filtered_counts[sample_ids]),
  denoised = as.numeric(denoised_counts[sample_ids]),
  nonchim = as.numeric(nonchim_counts[sample_ids])
)
track$filter_retention <- track$filtered / track$input
track$final_retention <- track$nonchim / track$input
write.table(track, file.path(runtime, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

asv_sequences <- colnames(seqtab_nochim)
asv_ids <- paste0(marker, "_ASV", seq_along(asv_sequences))
count_table <- t(seqtab_nochim)
rownames(count_table) <- asv_ids
write.table(
  data.frame(ASV_ID = rownames(count_table), count_table, check.names = FALSE),
  file.path(runtime, "asv_count_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
con <- file(file.path(runtime, "asv_sequences.fasta"), "w")
for (i in seq_along(asv_sequences)) {
  writeLines(c(paste0(">", asv_ids[[i]]), asv_sequences[[i]]), con)
}
close(con)

lengths <- nchar(asv_sequences)
summary_lines <- c(
  paste0("marker\t", marker),
  paste0("samples\t", nrow(seqtab_nochim)),
  paste0("ASVs_before_chimera\t", ncol(seqtab_all)),
  paste0("ASVs_after_chimera\t", ncol(seqtab_nochim)),
  paste0("reads_input\t", sum(track$input)),
  paste0("reads_filtered\t", sum(track$filtered)),
  paste0("reads_denoised\t", sum(track$denoised)),
  paste0("reads_nonchim\t", sum(track$nonchim)),
  paste0("final_retention\t", signif(sum(track$nonchim) / sum(track$input), 6)),
  paste0("ASV_length_min\t", min(lengths)),
  paste0("ASV_length_median\t", median(lengths)),
  paste0("ASV_length_max\t", max(lengths))
)
writeLines(summary_lines, file.path(runtime, "dada2_summary.tsv"))
message("Completed DADA2 PacBio pipeline for ", marker)

