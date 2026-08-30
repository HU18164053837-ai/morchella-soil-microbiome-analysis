#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(ggplot2))

project <- "."
source_dir <- file.path(project, "analysis", "hpc_results", "16S_dada2_v1_remote_original")
target_dir <- file.path(project, "analysis", "hpc_results", "16S_dada2_v1_corrected_v2")
diag_dir <- file.path(target_dir, "diagnostics_v1")

if (!dir.exists(source_dir)) stop("Missing source directory: ", source_dir)
if (dir.exists(target_dir)) stop("Refusing to overwrite target directory: ", target_dir)

source_files <- list.files(source_dir, recursive = TRUE, full.names = TRUE, all.files = TRUE, no.. = TRUE)
source_files <- source_files[file.info(source_files)$isdir %in% FALSE]
source_rel <- substring(source_files, nchar(source_dir) + 2L)
source_md5_before <- unname(tools::md5sum(source_files))
names(source_md5_before) <- source_rel

dir.create(target_dir, recursive = TRUE, showWarnings = FALSE)
for (rel in source_rel) {
  src <- file.path(source_dir, rel)
  dst <- file.path(target_dir, rel)
  dir.create(dirname(dst), recursive = TRUE, showWarnings = FALSE)
  if (!file.copy(src, dst, overwrite = FALSE, copy.date = TRUE)) stop("Failed to copy: ", rel)
}

strip_suffix <- function(x) sub("_16S_trimmed\\.fastq\\.gz$", "", x)
results_dir <- file.path(target_dir, "results")
checkpoint_dir <- file.path(target_dir, "checkpoints")

mapping_source <- read.delim(file.path(results_dir, "read_tracking.tsv"), check.names = FALSE,
                             stringsAsFactors = FALSE)$sample_id
mapping <- data.frame(original_sample_id = mapping_source,
                      normalized_sample_id = strip_suffix(mapping_source),
                      stringsAsFactors = FALSE)
if (nrow(mapping) != 24L || anyDuplicated(mapping$normalized_sample_id)) stop("Expected 24 unique normalized IDs")

metadata <- read.delim(file.path(project, "analysis", "metadata", "sample_metadata.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
expected_ids <- sort(unique(metadata$biological_sample_id[metadata$marker == "16S"]))
if (!identical(sort(mapping$normalized_sample_id), expected_ids)) stop("Normalized IDs do not match metadata")

for (filename in c("filter_stats.tsv", "read_tracking.tsv", "sample_depth_summary.tsv", "error_training_manifest.tsv")) {
  path <- file.path(results_dir, filename)
  tab <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  tab$sample_id <- strip_suffix(tab$sample_id)
  if (anyDuplicated(tab$sample_id)) stop("Duplicate IDs in ", filename)
  write.table(tab, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

count_path <- file.path(results_dir, "asv_count_table.tsv")
count_table <- read.delim(count_path, check.names = FALSE, stringsAsFactors = FALSE)
sample_columns <- setdiff(names(count_table), "ASV_ID")
normalized_columns <- strip_suffix(sample_columns)
if (length(sample_columns) != 24L || anyDuplicated(normalized_columns)) stop("Unexpected count-table samples")
names(count_table)[match(sample_columns, names(count_table))] <- normalized_columns
write.table(count_table, count_path, sep = "\t", quote = FALSE, row.names = FALSE)

normalize_seqtab <- function(path) {
  x <- readRDS(path)
  new_ids <- strip_suffix(rownames(x))
  if (nrow(x) < 1L || anyDuplicated(new_ids)) stop("Invalid seqtab IDs: ", path)
  rownames(x) <- new_ids
  saveRDS(x, path)
}
for (path in c(file.path(results_dir, "seqtab_all.rds"), file.path(results_dir, "seqtab_nochim.rds"),
               list.files(checkpoint_dir, pattern = "^seqtab_.*\\.rds$", full.names = TRUE))) normalize_seqtab(path)

for (path in list.files(checkpoint_dir, pattern = "^dada_.*\\.rds$", full.names = TRUE)) {
  x <- readRDS(path)
  names(x) <- strip_suffix(names(x))
  if (anyDuplicated(names(x))) stop("Duplicate dada IDs: ", path)
  saveRDS(x, path)
}

seqtab <- readRDS(file.path(results_dir, "seqtab_nochim.rds"))
counts <- as.matrix(count_table[, -1L, drop = FALSE])
storage.mode(counts) <- "numeric"
if (!identical(colnames(counts), rownames(seqtab))) stop("Count-table and seqtab sample order differ")
if (!isTRUE(all.equal(unname(counts), unname(t(seqtab)), check.attributes = FALSE))) stop("Counts differ")

tracking <- read.delim(file.path(results_dir, "read_tracking.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
tracked <- setNames(tracking$nonchim, tracking$sample_id)
if (!all(rowSums(seqtab) == tracked[rownames(seqtab)])) stop("Tracking totals differ")

write.table(mapping, file.path(target_dir, "sample_name_mapping.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)

asv_total <- colSums(seqtab)
asv_prev <- colSums(seqtab > 0)
asv_len <- nchar(colnames(seqtab))
asv_diag <- data.frame(
  ASV_ID = count_table$ASV_ID,
  sequence_length = asv_len,
  total_abundance = asv_total,
  prevalence_samples = asv_prev,
  prevalence_fraction = asv_prev / nrow(seqtab),
  max_sample_abundance = apply(seqtab, 2, max),
  mean_abundance_when_present = asv_total / asv_prev,
  global_relative_abundance = asv_total / sum(asv_total),
  singleton_read = asv_total == 1,
  low_abundance_le2 = asv_total <= 2,
  low_abundance_le5 = asv_total <= 5,
  low_abundance_le10 = asv_total <= 10,
  single_sample = asv_prev == 1,
  stringsAsFactors = FALSE
)
write.table(asv_diag, file.path(diag_dir, "asv_abundance_prevalence_length.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

thresholds <- c(1, 2, 5, 10, 20, 50, 100)
low_summary <- do.call(rbind, lapply(thresholds, function(z) data.frame(
  maximum_total_abundance = z,
  ASVs = sum(asv_total <= z),
  ASV_fraction = mean(asv_total <= z),
  reads = sum(asv_total[asv_total <= z]),
  read_fraction = sum(asv_total[asv_total <= z]) / sum(asv_total)
)))
write.table(low_summary, file.path(diag_dir, "low_abundance_threshold_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

prevalence_summary <- do.call(rbind, lapply(c(1, 2, 3, 6, 12, 24), function(z) data.frame(
  maximum_prevalence_samples = z,
  ASVs = sum(asv_prev <= z),
  ASV_fraction = mean(asv_prev <= z),
  reads = sum(asv_total[asv_prev <= z]),
  read_fraction = sum(asv_total[asv_prev <= z]) / sum(asv_total)
)))
write.table(prevalence_summary, file.path(diag_dir, "prevalence_threshold_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sample_diag <- data.frame(
  sample_id = rownames(seqtab),
  batch = tracking$batch[match(rownames(seqtab), tracking$sample_id)],
  input_reads = tracking$input[match(rownames(seqtab), tracking$sample_id)],
  filtered_reads = tracking$filtered[match(rownames(seqtab), tracking$sample_id)],
  nonchim_reads = rowSums(seqtab),
  observed_ASVs = rowSums(seqtab > 0),
  final_retention = tracking$final_retention[match(rownames(seqtab), tracking$sample_id)],
  reads_in_global_singleton_ASVs = rowSums(seqtab[, asv_total == 1, drop = FALSE]),
  reads_in_global_le5_ASVs = rowSums(seqtab[, asv_total <= 5, drop = FALSE]),
  low_depth_lt5000 = rowSums(seqtab) < 5000,
  low_retention_lt010 = tracking$final_retention[match(rownames(seqtab), tracking$sample_id)] < 0.10,
  stringsAsFactors = FALSE
)
write.table(sample_diag, file.path(diag_dir, "sample_depth_asv_diagnostics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

length_summary <- data.frame(
  metric = c("min", "q01", "q05", "median", "q95", "q99", "max", "outside_1300_1600"),
  value = c(min(asv_len), quantile(asv_len, .01), quantile(asv_len, .05), median(asv_len),
            quantile(asv_len, .95), quantile(asv_len, .99), max(asv_len), sum(asv_len < 1300 | asv_len > 1600))
)
write.table(length_summary, file.path(diag_dir, "asv_length_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

theme_qc <- theme_classic(base_size = 9, base_family = "sans") +
  theme(axis.line = element_line(linewidth = .35), axis.ticks = element_line(linewidth = .35),
        plot.title = element_text(face = "bold"), legend.title = element_blank())

save_png <- function(path, plot, width = 2274, height = 1447) {
  png(path, width = width, height = height, res = 350, type = "cairo", bg = "white")
  print(plot)
  dev.off()
}

p1 <- ggplot(asv_diag, aes(total_abundance, prevalence_samples)) +
  geom_point(alpha = .28, size = .7, colour = "#2878B5") +
  scale_x_log10(labels = scales::comma) +
  labs(title = "16S ASV abundance-prevalence diagnostic", x = "Total abundance (log10 scale)",
       y = "Samples detected (of 24)") + theme_qc
save_png(file.path(diag_dir, "QC_asv_abundance_prevalence.png"), p1)

p2 <- ggplot(asv_diag, aes(sequence_length)) +
  geom_histogram(binwidth = 5, boundary = 0, fill = "#57A6A1", colour = "white", linewidth = .15) +
  geom_vline(xintercept = c(1300, 1600), linetype = 2, colour = "#D9822B", linewidth = .45) +
  labs(title = "Full-length 16S ASV length distribution", x = "ASV length (bp)", y = "ASVs") + theme_qc
save_png(file.path(diag_dir, "QC_asv_length_distribution.png"), p2)

sample_diag$sample_id <- factor(sample_diag$sample_id, levels = sample_diag$sample_id[order(sample_diag$batch, sample_diag$nonchim_reads)])
p3 <- ggplot(sample_diag, aes(sample_id, nonchim_reads, fill = batch)) +
  geom_col(width = .78) + geom_hline(yintercept = 5000, linetype = 2, colour = "#C84A3A") +
  scale_y_continuous(labels = scales::comma) +
  scale_fill_manual(values = c(M = "#5B8DB8", A = "#E3A24A", DH = "#6BAF92")) +
  labs(title = "Non-chimeric 16S depth by sample", subtitle = "Dashed line: 5,000-read diagnostic threshold",
       x = NULL, y = "Non-chimeric reads") + theme_qc +
  theme(axis.text.x = element_text(angle = 55, hjust = 1), legend.position = "top")
save_png(file.path(diag_dir, "QC_sample_nonchim_depth.png"), p3, width = 2480, height = 1654)

p4 <- ggplot(low_summary, aes(factor(maximum_total_abundance), ASV_fraction * 100)) +
  geom_col(fill = "#8C8C8C", width = .72) +
  geom_text(aes(label = sprintf("%.1f%%", ASV_fraction * 100)), vjust = -.3, size = 3) +
  scale_y_continuous(expand = expansion(mult = c(0, .12))) +
  labs(title = "Fraction of ASVs below abundance thresholds", x = "Maximum total abundance", y = "ASVs (%)") + theme_qc
save_png(file.path(diag_dir, "QC_low_abundance_ASV_fraction.png"), p4)

batch_summary <- aggregate(cbind(input_reads, filtered_reads, nonchim_reads) ~ batch, sample_diag, sum)
batch_summary$filter_retention <- batch_summary$filtered_reads / batch_summary$input_reads
batch_summary$final_retention <- batch_summary$nonchim_reads / batch_summary$input_reads
write.table(batch_summary, file.path(diag_dir, "batch_depth_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

report <- c(
  "# 16S sample-name normalization and ASV diagnostic report",
  "", "## Data protection and normalization", "",
  "- The original downloaded directory was kept unchanged.",
  "- Only the terminal _16S_trimmed.fastq.gz suffix was removed from sample IDs.",
  paste0("- Samples: ", nrow(seqtab), "; ASVs: ", ncol(seqtab), "; non-chimeric reads: ", sum(seqtab), "."),
  "- IDs match metadata; the count table, RDS objects and tracking totals agree.",
  "", "## Diagnostic summary", "",
  paste0("- Samples below 5,000 non-chimeric reads: ", sum(sample_diag$low_depth_lt5000), "/24."),
  paste0("- Samples below 10% final retention: ", sum(sample_diag$low_retention_lt010), "/24."),
  paste0("- Global singleton ASVs: ", sum(asv_total == 1), "."),
  paste0("- ASVs with total abundance <=5: ", sum(asv_total <= 5), "."),
  paste0("- ASVs detected in one sample: ", sum(asv_prev == 1), "."),
  paste0("- Median ASV length: ", median(asv_len), " bp; outside 1300-1600 bp: ", sum(asv_len < 1300 | asv_len > 1600), "."),
  "", "No sample or ASV was removed in this diagnostic step."
)
writeLines(report, file.path(diag_dir, "16S_DIAGNOSTIC_REPORT.md"))

source_md5_after <- unname(tools::md5sum(source_files))
names(source_md5_after) <- source_rel
if (!identical(source_md5_before, source_md5_after)) stop("Original downloaded directory changed")

writeLines(c(
  "# 16S sample-name correction report", "",
  "- The original downloaded directory was not modified.",
  "- Sample IDs only were normalized; ASV sequences and abundances were unchanged.",
  "- Result RDS, checkpoint RDS, count-table columns and sample_id TSV fields were synchronized.",
  "- All 24 IDs match metadata and all original downloaded files retained their MD5 values."
), file.path(target_dir, "CORRECTION_REPORT.md"))

cat("NORMALIZATION_AND_DIAGNOSTICS_OK\n")
cat("samples=", nrow(seqtab), " asvs=", ncol(seqtab), " reads=", sum(seqtab), "\n", sep = "")
cat("low_depth=", sum(sample_diag$low_depth_lt5000), " low_retention=", sum(sample_diag$low_retention_lt010), "\n", sep = "")

