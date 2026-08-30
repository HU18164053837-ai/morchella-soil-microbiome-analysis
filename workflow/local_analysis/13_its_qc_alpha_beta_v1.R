#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(ape)
  library(ggplot2)
})

set.seed(20260821L)
project <- "."
source_dir <- file.path(project, "analysis", "hpc_results", "ITS_dada2_v1_corrected_v1")
result_dir <- file.path(project, "analysis", "downstream", "ITS_v1")
if (!dir.exists(source_dir)) stop("Missing corrected ITS source directory")
if (dir.exists(result_dir)) stop("Refusing to overwrite existing downstream directory")
dirs <- file.path(result_dir, c("qc", "diversity", "figures", "taxonomy"))
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)

count_tab <- read.delim(file.path(source_dir, "results", "asv_count_table.tsv"),
                        check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(count_tab$ASV_ID)) stop("Duplicated ASV IDs")
counts <- t(as.matrix(count_tab[, -1L, drop = FALSE]))
storage.mode(counts) <- "numeric"
colnames(counts) <- count_tab$ASV_ID

metadata_all <- read.delim(file.path(project, "analysis", "metadata", "sample_metadata.tsv"),
                           check.names = FALSE, stringsAsFactors = FALSE)
metadata <- metadata_all[metadata_all$marker == "ITS", ]
if (nrow(metadata) != 24L || anyDuplicated(metadata$biological_sample_id)) {
  stop("Expected 24 unique ITS metadata rows")
}
if (!setequal(rownames(counts), metadata$biological_sample_id)) {
  stop("Count-table samples do not match metadata")
}
metadata <- metadata[match(rownames(counts), metadata$biological_sample_id), ]
metadata$sample_id <- metadata$biological_sample_id

tracking <- read.delim(file.path(source_dir, "results", "read_tracking.tsv"),
                       check.names = FALSE, stringsAsFactors = FALSE)
tracking <- tracking[match(rownames(counts), tracking$sample_id), ]
if (anyNA(tracking$sample_id) || !all(rowSums(counts) == tracking$nonchim)) {
  stop("Tracking and count totals disagree")
}

qc <- merge(metadata, tracking, by = "sample_id", sort = FALSE,
            suffixes = c("_metadata", "_tracking"))
qc <- qc[match(rownames(counts), qc$sample_id), ]
qc$library_depth_flag <- qc$nonchim < 5000
qc$retention_flag <- qc$final_retention < 0.10
qc$metadata_missing_fields <- apply(qc, 1, function(x) sum(is.na(x) | trimws(x) == ""))
qc$analysis_status <- ifelse(qc$library_depth_flag | qc$retention_flag,
                             "review", "retain")
write.table(qc, file.path(result_dir, "qc", "sample_qc_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

alpha <- data.frame(
  sample_id = rownames(counts),
  observed_asv = specnumber(counts),
  shannon = diversity(counts, index = "shannon"),
  simpson = diversity(counts, index = "simpson"),
  invsimpson = diversity(counts, index = "invsimpson"),
  stringsAsFactors = FALSE
)
alpha$pielou <- alpha$shannon / log(alpha$observed_asv)
min_depth <- min(rowSums(counts))
rare_obs <- rare_shannon <- matrix(NA_real_, nrow(counts), 100L)
for (i in seq_len(100L)) {
  rare <- rrarefy(counts, sample = min_depth)
  rare_obs[, i] <- specnumber(rare)
  rare_shannon[, i] <- diversity(rare, index = "shannon")
}
alpha$rarefied_depth <- min_depth
alpha$rarefied_observed_mean <- rowMeans(rare_obs)
alpha$rarefied_observed_sd <- apply(rare_obs, 1, sd)
alpha$rarefied_shannon_mean <- rowMeans(rare_shannon)
alpha$rarefied_shannon_sd <- apply(rare_shannon, 1, sd)
alpha <- merge(alpha, metadata[, c("sample_id", "cohort", "condition",
                                   "continuous_cropping_years", "greenhouse_id",
                                   "spatial_point", "pair_id", "independent_unit")],
               by = "sample_id", sort = FALSE)
alpha <- alpha[match(rownames(counts), alpha$sample_id), ]
write.table(alpha, file.path(result_dir, "diversity", "alpha_diversity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ma_alpha <- alpha[grepl("^[MA]", alpha$sample_id), ]
ma_before <- ma_alpha[ma_alpha$condition == "before_cultivation", ]
ma_after <- ma_alpha[ma_alpha$condition != "before_cultivation", ]
ma_before <- ma_before[match(ma_after$pair_id, ma_before$pair_id), ]
if (anyNA(ma_before$pair_id) || !identical(ma_before$pair_id, ma_after$pair_id)) {
  stop("M/A pairing is incomplete")
}
paired_alpha <- data.frame(
  pair_id = ma_after$pair_id,
  greenhouse_id = ma_after$greenhouse_id,
  continuous_cropping_years = ma_after$continuous_cropping_years,
  before_sample = ma_before$sample_id,
  after_sample = ma_after$sample_id,
  observed_change = ma_after$observed_asv - ma_before$observed_asv,
  shannon_change = ma_after$shannon - ma_before$shannon,
  simpson_change = ma_after$simpson - ma_before$simpson,
  stringsAsFactors = FALSE
)
write.table(paired_alpha, file.path(result_dir, "diversity", "alpha_MA_paired_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

relative <- counts / rowSums(counts)
bray <- vegdist(relative, method = "bray")
jaccard <- vegdist(counts > 0, method = "jaccard", binary = TRUE)

pcoa_table <- function(distance, metric) {
  fit <- cmdscale(distance, k = 2, eig = TRUE, add = TRUE)
  axes <- fit$points[, 1:2, drop = FALSE]
  labels <- rownames(relative)
  if (length(labels) != nrow(axes)) stop("PCoA labels and coordinates differ")
  explained <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  data.frame(sample_id = labels, metric = metric,
             axis1 = axes[, 1], axis2 = axes[, 2],
             axis1_percent = explained[1], axis2_percent = explained[2],
             stringsAsFactors = FALSE)
}
pcoa <- rbind(pcoa_table(bray, "Bray-Curtis"), pcoa_table(jaccard, "Jaccard"))
pcoa <- merge(pcoa, metadata[, c("sample_id", "cohort", "condition",
                                 "continuous_cropping_years", "greenhouse_id",
                                 "spatial_point", "pair_id", "independent_unit")],
              by = "sample_id", sort = FALSE)
write.table(pcoa, file.path(result_dir, "diversity", "beta_pcoa_coordinates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ma_ids <- metadata$sample_id[grepl("^[MA]", metadata$sample_id)]
ma_meta <- metadata[match(ma_ids, metadata$sample_id), ]
ma_meta$condition <- factor(ma_meta$condition,
                            levels = c("before_cultivation", "after_cultivation"))
bray_ma <- as.dist(as.matrix(bray)[ma_ids, ma_ids])
perm_control <- how(nperm = 999, blocks = factor(ma_meta$pair_id))
perm <- adonis2(bray_ma ~ condition, data = ma_meta,
                permutations = perm_control, by = "margin")
permanova <- data.frame(
  comparison = "M_vs_A_exploratory_pair_blocked",
  distance = "Bray-Curtis_relative_abundance",
  term = "condition",
  df = perm$Df[1],
  R2 = perm$R2[1],
  F = perm$F[1],
  p_value_exploratory = perm$`Pr(>F)`[1],
  permutation_constraint = "within_pair_id",
  interpretation = "exploratory_only_nested_subsamples_not_independent_greenhouses",
  stringsAsFactors = FALSE
)
write.table(permanova, file.path(result_dir, "diversity", "permanova_exploratory.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pdf(file.path(result_dir, "figures", "rarefaction_curves.pdf"), width = 9, height = 7)
rarecurve(counts, step = 1000, sample = min_depth, label = FALSE,
          col = as.integer(factor(metadata$condition)),
          xlab = "Reads", ylab = "Observed ASVs")
legend("bottomright", legend = levels(factor(metadata$condition)),
       col = seq_along(levels(factor(metadata$condition))), lty = 1, cex = 0.7)
dev.off()

p_alpha <- ggplot(alpha, aes(x = condition, y = shannon, color = greenhouse_id)) +
  geom_point(size = 2.5) +
  facet_wrap(~cohort, scales = "free_x") +
  theme_bw(base_size = 11) +
  labs(x = NULL, y = "Shannon diversity", color = "Greenhouse") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(result_dir, "figures", "alpha_shannon_by_design.pdf"),
       p_alpha, width = 9, height = 5)

plot_pcoa <- pcoa[pcoa$metric == "Bray-Curtis", ]
p_beta <- ggplot(plot_pcoa, aes(axis1, axis2, color = condition, shape = greenhouse_id)) +
  geom_point(size = 3) + theme_bw(base_size = 11) +
  labs(x = sprintf("PCoA1 (%.1f%%)", unique(plot_pcoa$axis1_percent)),
       y = sprintf("PCoA2 (%.1f%%)", unique(plot_pcoa$axis2_percent)),
       color = "Condition", shape = "Greenhouse")
ggsave(file.path(result_dir, "figures", "beta_bray_pcoa.pdf"),
       p_beta, width = 8, height = 6)

qc_report <- c(
  "# ITS样本质控报告", "",
  paste0("- 样本数：", nrow(qc)), "- ASV数：8056",
  paste0("- 非嵌合reads范围：", min(qc$nonchim), "–", max(qc$nonchim)),
  paste0("- 最终保留率范围：", sprintf("%.2f%%", 100 * min(qc$final_retention)),
         "–", sprintf("%.2f%%", 100 * max(qc$final_retention))),
  paste0("- 低深度或低保留率风险样本：", sum(qc$analysis_status != "retain")),
  paste0("- 统一稀释深度：", min_depth, " reads；重复100次用于敏感性评估。"),
  "", "## 实验设计约束", "",
  "- M与A按pair_id配对；棚内3个空间亚样本不视为3个独立大棚。",
  "- 连作年限与大棚完全混杂，不进行独立年限因果推断。",
  "- H/D来自同一大棚两个点位，相关比较仅作描述性探索，不报告常规独立重复检验。",
  "- 当前不自动删除任何样本。"
)
writeLines(qc_report, file.path(result_dir, "qc", "QC_REPORT_ZH.md"), useBytes = TRUE)

methods <- c(
  "# ITS下游分析方法记录", "",
  "输入为DADA2去嵌合体ASV计数表，共24个样本、8056个ASV。",
  "Alpha多样性包括Observed ASV、Shannon、Simpson、inverse Simpson和Pielou；并在最小样本深度下进行100次随机稀释敏感性分析。",
  "Beta多样性包括相对丰度Bray-Curtis距离和presence/absence Jaccard距离，使用Lingoes校正PCoA。",
  "M/A的PERMANOVA仅作为探索性结果，置换限制在pair_id内；由于亚样本嵌套于大棚，不将其解释为独立大棚重复。",
  "H/D不执行常规PERMANOVA显著性检验。连作年限与大棚混杂，不作独立年限效应检验。",
  "软件：R 4.4.1、vegan 2.7-3、permute、ape、ggplot2；随机种子20260821。"
)
writeLines(methods, file.path(result_dir, "METHODS_ZH.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
cat("ITS_QC_ALPHA_BETA_OK\n")

