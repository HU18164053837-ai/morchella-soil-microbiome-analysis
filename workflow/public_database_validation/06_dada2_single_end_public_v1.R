#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(dada2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7L) stop("Usage: Rscript 06_dada2_single_end_public_v1.R project marker input_dir output_dir minLen maxLen threads")
project <- args[1]; marker <- args[2]; input_dir <- args[3]; output_dir <- args[4]
min_len <- as.integer(args[5]); max_len <- as.integer(args[6]); threads <- as.integer(args[7])
if (dir.exists(output_dir)) stop("Refusing to overwrite output directory: ", output_dir)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "filtered"), showWarnings = FALSE)

fn <- sort(list.files(input_dir, pattern = "\\.fastq\\.gz$", full.names = TRUE))
if (!length(fn)) stop("No FASTQ files in ", input_dir)
sample_id <- sub("\\.fastq\\.gz$", "", basename(fn))
if (anyDuplicated(sample_id)) stop("Duplicate sample IDs")
filt <- file.path(output_dir, "filtered", paste0(sample_id, ".filt.fastq.gz"))

filter_stats <- filterAndTrim(fn, filt, truncLen = 0, maxN = 0, maxEE = 2,
                              truncQ = 2, minLen = min_len, maxLen = max_len,
                              rm.phix = TRUE, compress = TRUE,
                              multithread = threads, verbose = TRUE)
filter_df <- data.frame(sample = rownames(filter_stats), input = filter_stats[, 1],
                        filtered = filter_stats[, 2], check.names = FALSE)
filter_df$retained_fraction <- filter_df$filtered / filter_df$input
write.table(filter_df, file.path(output_dir, "filter_stats.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

set.seed(20260822)
err <- learnErrors(filt, nbases = 1e8, randomize = TRUE,
                   multithread = threads, verbose = TRUE)
saveRDS(err, file.path(output_dir, "error_model.rds"))
pdf(file.path(output_dir, "error_model_diagnostic.pdf"), width = 7, height = 5, useDingbats = FALSE)
print(plotErrors(err, nominalQ = TRUE)); dev.off()

derep <- derepFastq(filt, verbose = TRUE)
names(derep) <- sample_id
dada_out <- dada(derep, err = err, pool = FALSE, multithread = threads, verbose = TRUE)
seqtab_all <- makeSequenceTable(dada_out)
seqtab_nochim <- removeBimeraDenovo(seqtab_all, method = "consensus",
                                    multithread = threads, verbose = TRUE)
saveRDS(seqtab_all, file.path(output_dir, "seqtab_all.rds"))
saveRDS(seqtab_nochim, file.path(output_dir, "seqtab_nochim.rds"))

get_unique <- function(x) sum(getUniques(x))
tracking <- data.frame(sample = sample_id,
                       input = filter_df$input[match(sample_id, filter_df$sample)],
                       filtered = filter_df$filtered[match(sample_id, filter_df$sample)],
                       denoised = vapply(dada_out, get_unique, numeric(1)),
                       nonchimera = rowSums(seqtab_nochim), check.names = FALSE)
write.table(tracking, file.path(output_dir, "read_tracking.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

asv_seq <- colnames(seqtab_nochim)
asv_id <- sprintf("ASV%06d", seq_along(asv_seq))
count_table <- as.data.frame(seqtab_nochim, check.names = FALSE)
colnames(count_table) <- asv_id
count_table <- cbind(sample = rownames(count_table), count_table)
write.table(count_table, file.path(output_dir, "asv_count_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
con <- file(file.path(output_dir, "asv_sequences.fasta"), "w")
for (i in seq_along(asv_seq)) writeLines(c(paste0(">", asv_id[i]), asv_seq[i]), con)
close(con)

depth <- data.frame(sample = rownames(seqtab_nochim), reads = rowSums(seqtab_nochim))
write.table(depth, file.path(output_dir, "sample_depth_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(paste("project", project), paste("marker", marker),
             paste("DADA2", as.character(packageVersion("dada2"))),
             paste("R", R.version.string), paste("minLen", min_len),
             paste("maxLen", max_len), "truncLen 0", "maxEE 2", "maxN 0", "truncQ 2", "rm.phix TRUE"),
           file.path(output_dir, "pipeline_versions.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
writeLines(c("DADA2 single-end public-validation run complete",
             paste("Project:", project), paste("Marker:", marker),
             paste("Samples:", nrow(seqtab_nochim)), paste("ASVs before chimera removal:", ncol(seqtab_all)),
             paste("ASVs after chimera removal:", ncol(seqtab_nochim)),
             sprintf("Overall filtered retention: %.4f", sum(filter_df$filtered) / sum(filter_df$input)),
             sprintf("Overall nonchimera retention: %.4f", sum(seqtab_nochim) / sum(filter_df$input))),
           file.path(output_dir, "run_summary.txt"))


