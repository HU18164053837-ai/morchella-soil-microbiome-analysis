#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

set.seed(20260821L)
project <- "."
source_dir <- file.path(project, "analysis", "downstream", "ITS_v2_taxonomy", "tables")
count_file <- file.path(source_dir, "fungal_asv_count_table.tsv")
tax_file <- file.path(source_dir, "fungal_asv_taxonomy.tsv")
metadata_file <- file.path(project, "analysis", "metadata", "sample_metadata.tsv")
result_dir <- file.path(project, "analysis", "downstream", "ITS_v4_genus_screening")
if (dir.exists(result_dir)) stop("Refusing to overwrite ITS_v4_genus_screening")
for (d in file.path(result_dir, c("tables", "figures", "qc"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

count_tab <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.delim(tax_file, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = character())
meta_all <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta_all[meta_all$marker == "ITS", ]
stopifnot(nrow(count_tab) == nrow(tax), !anyDuplicated(count_tab$ASV_ID),
          !anyDuplicated(tax$ASV_ID), setequal(count_tab$ASV_ID, tax$ASV_ID),
          nrow(meta) == 24L)
tax <- tax[match(count_tab$ASV_ID, tax$ASV_ID), ]

asv_counts <- t(as.matrix(count_tab[, -1L, drop = FALSE]))
storage.mode(asv_counts) <- "numeric"
colnames(asv_counts) <- count_tab$ASV_ID
stopifnot(setequal(rownames(asv_counts), meta$biological_sample_id))
meta <- meta[match(rownames(asv_counts), meta$biological_sample_id), ]
meta$sample_id <- meta$biological_sample_id

genus <- tax$Genus
genus[!nzchar(genus)] <- "Unclassified_at_genus"
genus_counts <- t(rowsum(t(asv_counts), group = genus, reorder = TRUE))
stopifnot(all(rowSums(genus_counts) == rowSums(asv_counts)))
genus_rel <- genus_counts / rowSums(genus_counts)
genus_clr <- log2(genus_counts + 0.5)
genus_clr <- genus_clr - rowMeans(genus_clr)

genus_out <- data.frame(sample_id = rownames(genus_counts), genus_counts,
                        check.names = FALSE)
write.table(genus_out, file.path(result_dir, "tables", "fungal_genus_count_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
rel_out <- data.frame(sample_id = rownames(genus_rel), genus_rel, check.names = FALSE)
write.table(rel_out, file.path(result_dir, "tables", "fungal_genus_relative_abundance.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ma_meta <- meta[meta$cohort == "cultivation_paired", ]
before_meta <- ma_meta[ma_meta$condition == "before_cultivation", ]
after_meta <- ma_meta[ma_meta$condition == "after_cultivation", ]
before_meta <- before_meta[match(after_meta$pair_id, before_meta$pair_id), ]
stopifnot(nrow(before_meta) == 9L, nrow(after_meta) == 9L,
          identical(before_meta$pair_id, after_meta$pair_id))
b_ids <- before_meta$sample_id
a_ids <- after_meta$sample_id

ma_rows <- vector("list", ncol(genus_counts))
gh_rows <- list()
for (j in seq_len(ncol(genus_counts))) {
  g <- colnames(genus_counts)[j]
  delta <- genus_clr[a_ids, j] - genus_clr[b_ids, j]
  gh <- aggregate(delta, by = list(greenhouse_id = after_meta$greenhouse_id), FUN = mean)
  names(gh)[2] <- "mean_paired_clr_change"
  gh$genus <- g
  gh_rows[[j]] <- gh[, c("genus", "greenhouse_id", "mean_paired_clr_change")]
  p <- tryCatch(wilcox.test(genus_clr[a_ids, j], genus_clr[b_ids, j],
                            paired = TRUE, exact = FALSE)$p.value,
                error = function(e) NA_real_)
  pos <- sum(delta > 0); neg <- sum(delta < 0)
  gh_signs <- sign(gh$mean_paired_clr_change)
  gh_consistent <- all(gh_signs > 0) || all(gh_signs < 0)
  ma_rows[[j]] <- data.frame(
    genus = g, total_reads_24 = sum(genus_counts[, j]),
    prevalence_MA = sum(genus_counts[c(b_ids, a_ids), j] > 0),
    mean_relative_before = mean(genus_rel[b_ids, j]),
    mean_relative_after = mean(genus_rel[a_ids, j]),
    log2_mean_relative_ratio = log2((mean(genus_rel[a_ids, j]) + 1e-6) /
                                      (mean(genus_rel[b_ids, j]) + 1e-6)),
    mean_paired_clr_change = mean(delta), median_paired_clr_change = median(delta),
    positive_pairs = pos, negative_pairs = neg, zero_pairs = sum(delta == 0),
    pair_direction_consistency = max(pos, neg) / length(delta),
    greenhouse_direction_consistent = gh_consistent,
    wilcoxon_p_exploratory = p, stringsAsFactors = FALSE)
}
ma_screen <- do.call(rbind, ma_rows)
ma_screen$BH_q_exploratory <- p.adjust(ma_screen$wilcoxon_p_exploratory, method = "BH")
eligible <- ma_screen$genus != "Unclassified_at_genus" &
  ma_screen$prevalence_MA >= 6L & ma_screen$total_reads_24 >= 100L
ma_screen$priority <- "not_prioritized"
moderate <- eligible & abs(ma_screen$median_paired_clr_change) >= 0.5 &
  abs(ma_screen$log2_mean_relative_ratio) >= 0.5 &
  sign(ma_screen$median_paired_clr_change) == sign(ma_screen$log2_mean_relative_ratio) &
  ma_screen$pair_direction_consistency >= 6/9 & ma_screen$greenhouse_direction_consistent
high <- eligible & abs(ma_screen$median_paired_clr_change) >= 1 &
  abs(ma_screen$log2_mean_relative_ratio) >= 1 &
  sign(ma_screen$median_paired_clr_change) == sign(ma_screen$log2_mean_relative_ratio) &
  ma_screen$pair_direction_consistency >= 7/9 & ma_screen$greenhouse_direction_consistent
ma_screen$priority[moderate] <- "moderate"
ma_screen$priority[high] <- "high"
ma_screen <- ma_screen[order(factor(ma_screen$priority,
                                    levels = c("high", "moderate", "not_prioritized")),
                             -abs(ma_screen$median_paired_clr_change)), ]
write.table(ma_screen, file.path(result_dir, "tables", "MA_genus_screening.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

gh_effects <- do.call(rbind, gh_rows)
write.table(gh_effects, file.path(result_dir, "tables", "MA_genus_greenhouse_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ma_screen[ma_screen$priority != "not_prioritized", ],
            file.path(result_dir, "tables", "MA_prioritized_genera.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

h_ids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "healthy"]
d_ids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "diseased"]
stopifnot(length(h_ids) == 3L, length(d_ids) == 3L)
hd_rows <- vector("list", ncol(genus_counts))
for (j in seq_len(ncol(genus_counts))) {
  h <- genus_rel[h_ids, j]; d <- genus_rel[d_ids, j]
  h_mean <- mean(h); d_mean <- mean(d)
  effect <- log2((d_mean + 1e-6) / (h_mean + 1e-6))
  hd_rows[[j]] <- data.frame(
    genus = colnames(genus_counts)[j], total_reads_24 = sum(genus_counts[, j]),
    detected_healthy = sum(genus_counts[h_ids, j] > 0),
    detected_diseased = sum(genus_counts[d_ids, j] > 0),
    mean_relative_healthy = h_mean, mean_relative_diseased = d_mean,
    log2_mean_relative_ratio_D_vs_H = effect,
    diseased_above_healthy_median = sum(d > median(h)),
    diseased_below_healthy_median = sum(d < median(h)), stringsAsFactors = FALSE)
}
hd_screen <- do.call(rbind, hd_rows)
hd_screen$priority <- "not_prioritized"
hd_eligible <- hd_screen$genus != "Unclassified_at_genus" &
  (hd_screen$detected_healthy + hd_screen$detected_diseased) >= 2L &
  pmax(hd_screen$mean_relative_healthy, hd_screen$mean_relative_diseased) >= 0.001 &
  abs(hd_screen$log2_mean_relative_ratio_D_vs_H) >= 1
hd_screen$priority[hd_eligible] <- "descriptive_candidate"
hd_robust <- hd_eligible &
  pmax(hd_screen$mean_relative_healthy, hd_screen$mean_relative_diseased) >= 0.005 &
  abs(hd_screen$log2_mean_relative_ratio_D_vs_H) >= 2 &
  ((hd_screen$log2_mean_relative_ratio_D_vs_H > 0 & hd_screen$diseased_above_healthy_median == 3L) |
   (hd_screen$log2_mean_relative_ratio_D_vs_H < 0 & hd_screen$diseased_below_healthy_median == 3L))
hd_screen$priority[hd_robust] <- "robust_descriptive_candidate"
hd_screen <- hd_screen[order(factor(hd_screen$priority,
                                    levels = c("robust_descriptive_candidate",
                                               "descriptive_candidate", "not_prioritized")),
                             -abs(hd_screen$log2_mean_relative_ratio_D_vs_H)), ]
write.table(hd_screen, file.path(result_dir, "tables", "HD_genus_screening.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(hd_screen[hd_screen$priority != "not_prioritized", ],
            file.path(result_dir, "tables", "HD_descriptive_candidate_genera.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ma_plot <- ma_screen[ma_screen$genus != "Unclassified_at_genus", ]
label_ma <- head(ma_plot[order(-abs(ma_plot$median_paired_clr_change)), ], 12L)$genus
p_ma <- ggplot(ma_plot, aes(median_paired_clr_change, pair_direction_consistency,
                            color = priority, size = log10(total_reads_24 + 1))) +
  geom_hline(yintercept = 6/9, linetype = 2, color = "grey60") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey60") +
  geom_point(alpha = 0.75) +
  geom_text(data = ma_plot[ma_plot$genus %in% label_ma, ], aes(label = genus),
            size = 2.6, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  theme_bw(base_size = 10) +
  labs(x = "Median paired CLR change (After - Before)",
       y = "Pair direction consistency", color = "Priority", size = "log10 reads")
ggsave(file.path(result_dir, "figures", "MA_genus_effect_screen.pdf"),
       p_ma, width = 9, height = 6)

hd_plot <- hd_screen[hd_screen$genus != "Unclassified_at_genus" &
                       pmax(hd_screen$mean_relative_healthy,
                            hd_screen$mean_relative_diseased) > 0, ]
label_hd <- head(hd_plot[order(-abs(hd_plot$log2_mean_relative_ratio_D_vs_H)), ], 12L)$genus
p_hd <- ggplot(hd_plot, aes(log2_mean_relative_ratio_D_vs_H,
                            pmax(mean_relative_healthy, mean_relative_diseased),
                            color = priority, size = log10(total_reads_24 + 1))) +
  geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey60") +
  geom_point(alpha = 0.75) + scale_y_log10() +
  geom_text(data = hd_plot[hd_plot$genus %in% label_hd, ], aes(label = genus),
            size = 2.6, vjust = -0.7, check_overlap = TRUE, show.legend = FALSE) +
  theme_bw(base_size = 10) +
  labs(x = "log2 mean relative abundance ratio (Diseased / Healthy)",
       y = "Maximum mean relative abundance", color = "Priority", size = "log10 reads")
ggsave(file.path(result_dir, "figures", "HD_genus_effect_screen.pdf"),
       p_hd, width = 9, height = 6)

metrics <- data.frame(
  metric = c("fungal_genera_total", "MA_high_priority", "MA_moderate_priority",
             "HD_robust_descriptive_candidates", "HD_other_descriptive_candidates"),
  value = c(ncol(genus_counts), sum(ma_screen$priority == "high"),
            sum(ma_screen$priority == "moderate"),
            sum(hd_screen$priority == "robust_descriptive_candidate"),
            sum(hd_screen$priority == "descriptive_candidate")),
  stringsAsFactors = FALSE)
write.table(metrics, file.path(result_dir, "qc", "screening_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
cat("ITS_GENUS_SCREENING_V1_OK\n")

