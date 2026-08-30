#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dada2)
  library(ShortRead)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: Rscript 10_dada2_16s_formal_v1.R <project_root> <threads>")
project <- normalizePath(args[[1]], mustWork = TRUE)
threads <- as.integer(args[[2]])
if (is.na(threads) || threads < 1L) stop("threads must be a positive integer")

marker <- "16S"
seed <- 20260820L
reads_per_sample <- 10000L
input_dir <- file.path(project, "input", "trimmed", marker)
workdir <- file.path(project, "work", "dada2_16S_v1")
filtered_dir <- file.path(workdir, "filtered")
training_root <- file.path(workdir, "error_training")
outdir <- file.path(project, "results", "dada2_16S_v1")
checkpoint_dir <- file.path(project, "checkpoints", "dada2_16S_v1")

for (d in c(workdir, outdir, checkpoint_dir)) {
  if (dir.exists(d)) stop("Refusing to overwrite existing directory: ", d)
}
dir.create(filtered_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

files <- sort(list.files(input_dir, pattern = "_16S_trimmed.fastq.gz$", full.names = TRUE))
if (length(files) != 24L) stop("Expected 24 trimmed 16S FASTQ files, found ", length(files))
sample_ids <- sub("_16S_trimmed.fastq.gz$", "", basename(files), fixed = TRUE)
names(files) <- sample_ids
filtered <- file.path(filtered_dir, paste0(sample_ids, "_16S_filtered.fastq.gz"))
names(filtered) <- sample_ids

batch_for_sample <- function(x) {
  if (grepl("^M", x)) "M" else if (grepl("^A", x)) "A" else "DH"
}
batches <- vapply(sample_ids, batch_for_sample, character(1))
expected_batch_sizes <- c(M = 9L, A = 9L, DH = 6L)
observed_batch_sizes <- table(factor(batches, levels = names(expected_batch_sizes)))
if (!all(as.integer(observed_batch_sizes) == expected_batch_sizes)) {
  stop("Unexpected batch sizes: ", paste(names(observed_batch_sizes), observed_batch_sizes, collapse = ", "))
}

message("Filtering 24 full-length 16S FASTQ files")
filter_matrix <- filterAndTrim(
  unname(files), unname(filtered),
  maxN = 0, maxEE = 2, truncQ = 2,
  minLen = 1200, maxLen = 1800,
  rm.phix = FALSE, compress = TRUE,
  multithread = threads, verbose = TRUE
)
filter_stats <- data.frame(
  sample_id = sample_ids,
  reads.in = filter_matrix[, "reads.in"],
  reads.out = filter_matrix[, "reads.out"],
  stringsAsFactors = FALSE
)
write.table(filter_stats, file.path(outdir, "filter_stats.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(seed)
training_files <- setNames(vector("list", length(sample_ids)), sample_ids)
training_manifest <- list()
for (id in sample_ids) {
  batch <- batches[[id]]
  batch_dir <- file.path(training_root, batch)
  dir.create(batch_dir, recursive = TRUE, showWarnings = FALSE)
  fq <- readFastq(filtered[[id]])
  n_available <- length(fq)
  n_take <- min(reads_per_sample, n_available)
  if (n_take < 1000L) stop("Too few filtered reads for balanced training: ", id)
  idx <- sort(sample.int(n_available, n_take, replace = FALSE))
  target <- file.path(batch_dir, paste0(id, "_16S_errortrain_", n_take, ".fastq.gz"))
  writeFastq(fq[idx], target, compress = TRUE)
  training_files[[id]] <- target
  training_manifest[[id]] <- data.frame(
    sample_id = id, batch = batch, filtered_reads = n_available,
    training_reads = n_take, training_file = target,
    stringsAsFactors = FALSE
  )
  rm(fq)
  gc()
}
training_manifest <- do.call(rbind, training_manifest)
write.table(training_manifest, file.path(outdir, "error_training_manifest.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

seqtabs <- list()
denoised_counts <- setNames(numeric(length(sample_ids)), sample_ids)
batch_runtime <- list()

for (batch in c("M", "A", "DH")) {
  ids <- sample_ids[batches == batch]
  batch_training <- unname(unlist(training_files[ids]))
  batch_filtered <- unname(filtered[ids])

  learn_started <- Sys.time()
  err <- learnErrors(
    batch_training,
    nbases = 1e8,
    errorEstimationFunction = PacBioErrfun,
    BAND_SIZE = 32,
    multithread = threads,
    randomize = TRUE,
    verbose = TRUE
  )
  learn_seconds <- as.numeric(difftime(Sys.time(), learn_started, units = "secs"))
  saveRDS(err, file.path(checkpoint_dir, paste0("error_model_", batch, ".rds")))
  saveRDS(err, file.path(outdir, paste0("error_model_", batch, ".rds")))
  pdf(file.path(outdir, paste0("error_model_", batch, ".pdf")), width = 8, height = 6)
  print(plotErrors(err, nominalQ = TRUE))
  dev.off()

  infer_started <- Sys.time()
  dada_out <- dada(
    batch_filtered,
    err = err,
    pool = "pseudo",
    BAND_SIZE = 32,
    multithread = threads,
    verbose = TRUE
  )
  names(dada_out) <- ids
  infer_seconds <- as.numeric(difftime(Sys.time(), infer_started, units = "secs"))
  seqtab <- makeSequenceTable(dada_out)
  seqtabs[[batch]] <- seqtab
  denoised_counts[ids] <- vapply(dada_out, function(x) sum(getUniques(x)), numeric(1))
  saveRDS(dada_out, file.path(checkpoint_dir, paste0("dada_", batch, ".rds")))
  saveRDS(seqtab, file.path(checkpoint_dir, paste0("seqtab_", batch, ".rds")))
  batch_runtime[[batch]] <- data.frame(
    batch = batch, samples = length(ids),
    learn_errors_seconds = learn_seconds,
    pseudo_inference_seconds = infer_seconds,
    stringsAsFactors = FALSE
  )
  rm(err, dada_out)
  gc()
}

seqtab_all <- Reduce(function(x, y) mergeSequenceTables(x, y), seqtabs)
seqtab_nochim <- removeBimeraDenovo(
  seqtab_all, method = "consensus", multithread = threads, verbose = TRUE
)
saveRDS(seqtab_all, file.path(outdir, "seqtab_all.rds"))
saveRDS(seqtab_nochim, file.path(outdir, "seqtab_nochim.rds"))
saveRDS(seqtab_all, file.path(checkpoint_dir, "seqtab_all.rds"))
saveRDS(seqtab_nochim, file.path(checkpoint_dir, "seqtab_nochim.rds"))

input_counts <- setNames(filter_stats$reads.in, filter_stats$sample_id)
filtered_counts <- setNames(filter_stats$reads.out, filter_stats$sample_id)
nonchim_counts <- setNames(rep(0, length(sample_ids)), sample_ids)
nonchim_counts[rownames(seqtab_nochim)] <- rowSums(seqtab_nochim)
track <- data.frame(
  sample_id = sample_ids,
  marker = marker,
  batch = batches,
  input = as.numeric(input_counts[sample_ids]),
  filtered = as.numeric(filtered_counts[sample_ids]),
  denoised = as.numeric(denoised_counts[sample_ids]),
  nonchim = as.numeric(nonchim_counts[sample_ids]),
  stringsAsFactors = FALSE
)
track$filter_retention <- track$filtered / track$input
track$denoise_retention <- track$denoised / track$filtered
track$final_retention <- track$nonchim / track$input
write.table(track, file.path(outdir, "read_tracking.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

depth_summary <- data.frame(
  sample_id = track$sample_id,
  batch = track$batch,
  nonchim_reads = track$nonchim,
  final_retention = track$final_retention,
  low_depth_flag = track$nonchim < 5000,
  low_retention_flag = track$final_retention < 0.10,
  stringsAsFactors = FALSE
)
write.table(depth_summary, file.path(outdir, "sample_depth_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

asv_sequences <- colnames(seqtab_nochim)
asv_ids <- paste0("16S_ASV", seq_along(asv_sequences))
count_table <- t(seqtab_nochim)
rownames(count_table) <- asv_ids
write.table(data.frame(ASV_ID = rownames(count_table), count_table, check.names = FALSE),
            file.path(outdir, "asv_count_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
con <- file(file.path(outdir, "asv_sequences.fasta"), "w")
for (i in seq_along(asv_sequences)) {
  writeLines(c(paste0(">", asv_ids[[i]]), asv_sequences[[i]]), con)
}
close(con)

batch_runtime <- do.call(rbind, batch_runtime)
write.table(batch_runtime, file.path(outdir, "batch_runtime.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

lengths <- nchar(asv_sequences)
summary_lines <- c(
  "metric\tvalue",
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
writeLines(summary_lines, file.path(outdir, "dada2_summary.tsv"))

writeLines(c(
  paste0("R\t", R.version.string),
  paste0("dada2\t", as.character(packageVersion("dada2"))),
  paste0("ShortRead\t", as.character(packageVersion("ShortRead"))),
  "container\tbioconductor-dada2_1.34.0--r44he5774e6_0.sif"
), file.path(outdir, "pipeline_versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

flagged <- depth_summary$sample_id[depth_summary$low_depth_flag | depth_summary$low_retention_flag]
report <- c(
  "# 羊肚菌土壤微生物组16S DADA2正式运行报告",
  "",
  "## 方法",
  "",
  "采用PacBioErrfun、BAND_SIZE=32和pseudo-pooling。M、A、DH三个批次分别按每个样本最多10000条过滤后序列进行平衡误差训练；样本推断后合并序列表，并以consensus方法去除嵌合体。",
  "",
  "过滤参数：maxN=0、maxEE=2、truncQ=2、minLen=1200、maxLen=1800、rm.phix=FALSE。",
  "",
  "## 结果摘要",
  "",
  paste0("输入reads：", sum(track$input)),
  paste0("过滤后reads：", sum(track$filtered)),
  paste0("去噪后reads：", sum(track$denoised)),
  paste0("非嵌合reads：", sum(track$nonchim)),
  paste0("非嵌合ASV：", ncol(seqtab_nochim)),
  "",
  "## 风险样本",
  "",
  if (length(flagged)) paste(flagged, collapse = "、") else "无",
  "",
  "低深度或低保留率样本不在本步骤自动删除，需结合后续物种注释和样本深度分布决定是否纳入统计分析。"
)
writeLines(report, file.path(outdir, "run_report_zh.md"), useBytes = TRUE)
message("Completed formal PacBio 16S DADA2 pipeline")

