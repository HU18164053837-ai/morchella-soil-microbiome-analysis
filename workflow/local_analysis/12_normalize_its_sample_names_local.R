#!/usr/bin/env Rscript

project <- "."
source_dir <- file.path(
  project, "analysis", "hpc_results", "ITS_dada2_v1_remote_original"
)
target_dir <- file.path(
  project, "analysis", "hpc_results", "ITS_dada2_v1_corrected_v1"
)

if (!dir.exists(source_dir)) stop("Missing source directory: ", source_dir)
if (dir.exists(target_dir)) stop("Refusing to overwrite target directory: ", target_dir)

source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE,
                           all.files = TRUE, no.. = TRUE)
source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
source_rel <- substring(source_files, nchar(source_dir) + 2L)
source_md5_before <- unname(tools::md5sum(source_files))
names(source_md5_before) <- source_rel

dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
for (rel in source_rel) {
  src <- file.path(source_dir, rel)
  dst <- file.path(target_dir, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dst, overwrite = FALSE, copy.date = TRUE)) {
    stop("Failed to copy: ", rel)
  }
}

strip_suffix <- function(x) sub("_ITS_trimmed\\.fastq\\.gz$", "", x)
results_dir <- file.path(target_dir, "results")

mapping_source <- read.delim(file.path(results_dir, "read_tracking.tsv"),
                             check.names = FALSE, stringsAsFactors = FALSE)$sample_id
mapping <- data.frame(
  original_sample_id = mapping_source,
  normalized_sample_id = strip_suffix(mapping_source),
  stringsAsFactors = FALSE
)
if (nrow(mapping) != 24L || anyDuplicated(mapping$normalized_sample_id)) {
  stop("Expected 24 unique normalized sample IDs")
}

metadata <- read.delim(file.path(project, "analysis", "metadata", "sample_metadata.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
expected_ids <- sort(unique(metadata$biological_sample_id[metadata$marker == "ITS"]))
if (!identical(sort(mapping$normalized_sample_id), expected_ids)) {
  stop("Normalized ITS IDs do not match sample_metadata.tsv")
}

sample_id_tables <- c(
  "filter_stats.tsv",
  "read_tracking.tsv",
  "sample_depth_summary.tsv",
  "error_training_manifest.tsv"
)
for (filename in sample_id_tables) {
  path <- file.path(results_dir, filename)
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"sample_id" %in% names(tab)) stop("Missing sample_id in ", filename)
  tab$sample_id <- strip_suffix(tab$sample_id)
  if (anyDuplicated(tab$sample_id)) stop("Duplicate normalized IDs in ", filename)
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

count_path <- file.path(results_dir, "asv_count_table.tsv")
count_table <- read.delim(count_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_columns <- setdiff(names(count_table), "ASV_ID")
normalized_columns <- strip_suffix(sample_columns)
if (length(sample_columns) != 24L || anyDuplicated(normalized_columns)) {
  stop("Unexpected sample columns in asv_count_table.tsv")
}
names(count_table)[match(sample_columns, names(count_table))] <- normalized_columns
write.table(count_table, count_path, sep = "\t", quote = FALSE, row.names = FALSE)

for (filename in c("seqtab_all.rds", "seqtab_nochim.rds")) {
  path <- file.path(results_dir, filename)
  seqtab <- readRDS(path)
  new_ids <- strip_suffix(rownames(seqtab))
  if (nrow(seqtab) != 24L || anyDuplicated(new_ids)) {
    stop("Unexpected sample rows in ", filename)
  }
  rownames(seqtab) <- new_ids
  saveRDS(seqtab, path)
}

seqtab_nochim <- readRDS(file.path(results_dir, "seqtab_nochim.rds"))
corrected_counts <- read.delim(count_path, check.names = FALSE,
                               stringsAsFactors = FALSE)
count_matrix <- as.matrix(corrected_counts[, -1L, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
if (!identical(colnames(count_matrix), rownames(seqtab_nochim))) {
  stop("Sample order differs between count table and seqtab_nochim.rds")
}
if (!isTRUE(all.equal(unname(count_matrix), unname(t(seqtab_nochim)),
                      check.attributes = FALSE))) {
  stop("Counts differ between count table and seqtab_nochim.rds")
}

tracking <- read.delim(file.path(results_dir, "read_tracking.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
tracked_nonchim <- setNames(tracking$nonchim, tracking$sample_id)
if (!setequal(rownames(seqtab_nochim), tracking$sample_id)) {
  stop("Sample IDs differ between tracking and seqtab_nochim.rds")
}
if (!all(rowSums(seqtab_nochim) == tracked_nonchim[rownames(seqtab_nochim)])) {
  stop("Nonchimera read totals differ between tracking and seqtab_nochim.rds")
}

write.table(mapping, file.path(target_dir, "sample_name_mapping.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

source_md5_after <- unname(tools::md5sum(source_files))
names(source_md5_after) <- source_rel
if (!identical(source_md5_before, source_md5_after)) {
  stop("Original downloaded directory changed during correction")
}

report <- c(
  "# ITS样本名规范化报告",
  "",
  "- 原始下载目录保持不变：`ITS_dada2_v1_remote_original`。",
  "- 修正版目录：`ITS_dada2_v1_corrected_v1`。",
  "- 规则：仅删除样本名末尾的 `_ITS_trimmed.fastq.gz`。",
  paste0("- 已规范化样本数：", nrow(mapping), "。"),
  paste0("- ASV数：", nrow(count_matrix), "。"),
  paste0("- 非嵌合reads总数：", sum(count_matrix), "。"),
  "- 已验证计数表与 `seqtab_nochim.rds` 数值、样本顺序完全一致。",
  "- 已验证 `read_tracking.tsv` 的nonchim计数与序列表逐样本一致。",
  "- 已验证24个规范化样本名与 `sample_metadata.tsv` 完全匹配。",
  "- 已验证修正前后原始下载目录全部文件MD5不变。",
  "- ASV序列、丰度数值、误差模型、日志及Slurm脚本均未改变。"
)
writeLines(report, file.path(target_dir, "CORRECTION_REPORT.md"), useBytes = TRUE)

cat("NORMALIZATION_OK\n")
cat("samples=", nrow(mapping), "\n", sep = "")
cat("asvs=", nrow(count_matrix), "\n", sep = "")
cat("nonchim_reads=", sum(count_matrix), "\n", sep = "")
cat("source_files_unchanged=", length(source_files), "\n", sep = "")

