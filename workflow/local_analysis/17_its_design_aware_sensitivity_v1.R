#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(ggplot2)
})

set.seed(20260821L)
project <- "."
source_dir <- file.path(project, "analysis", "downstream", "ITS_v2_taxonomy")
count_file <- file.path(project, "analysis", "hpc_results", "ITS_dada2_v1_corrected_v1",
                        "results", "asv_count_table.tsv")
fungal_count_file <- file.path(source_dir, "tables", "fungal_asv_count_table.tsv")
metadata_file <- file.path(project, "analysis", "metadata", "sample_metadata.tsv")
result_dir <- file.path(project, "analysis", "downstream", "ITS_v3_design_aware")
if (dir.exists(result_dir)) stop("Refusing to overwrite ITS_v3_design_aware")
for (d in file.path(result_dir, c("alpha", "beta", "figures", "qc"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

read_counts <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(x$ASV_ID))
  m <- t(as.matrix(x[, -1L, drop = FALSE]))
  storage.mode(m) <- "numeric"
  colnames(m) <- x$ASV_ID
  m
}
all_counts <- read_counts(count_file)
fungal_counts <- read_counts(fungal_count_file)
stopifnot(identical(rownames(all_counts), rownames(fungal_counts)))

meta_all <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta_all[meta_all$marker == "ITS", ]
stopifnot(nrow(meta) == 24L, !anyDuplicated(meta$biological_sample_id),
          setequal(rownames(all_counts), meta$biological_sample_id))
meta <- meta[match(rownames(all_counts), meta$biological_sample_id), ]
meta$sample_id <- meta$biological_sample_id

alpha_one <- function(x, subset_name) {
  data.frame(sample_id = rownames(x), subset = subset_name,
             depth = rowSums(x), observed_asv = specnumber(x),
             shannon = diversity(x, "shannon"),
             simpson = diversity(x, "simpson"), stringsAsFactors = FALSE)
}
alpha <- rbind(alpha_one(all_counts, "all_ASVs"),
               alpha_one(fungal_counts, "UNITE_kingdom_Fungi"))
alpha <- merge(alpha, meta[, c("sample_id", "cohort", "condition",
                               "continuous_cropping_years", "greenhouse_id",
                               "spatial_point", "pair_id", "independent_unit")],
               by = "sample_id", sort = FALSE)

ma <- alpha[alpha$cohort == "cultivation_paired", ]
paired_rows <- list()
for (ss in unique(ma$subset)) {
  z <- ma[ma$subset == ss, ]
  before <- z[z$condition == "before_cultivation", ]
  after <- z[z$condition == "after_cultivation", ]
  before <- before[match(after$pair_id, before$pair_id), ]
  stopifnot(nrow(before) == 9L, nrow(after) == 9L,
            identical(before$pair_id, after$pair_id))
  paired_rows[[ss]] <- data.frame(
    subset = ss, pair_id = after$pair_id, greenhouse_id = after$greenhouse_id,
    continuous_cropping_years = after$continuous_cropping_years,
    before_sample = before$sample_id, after_sample = after$sample_id,
    observed_before = before$observed_asv, observed_after = after$observed_asv,
    observed_change = after$observed_asv - before$observed_asv,
    shannon_before = before$shannon, shannon_after = after$shannon,
    shannon_change = after$shannon - before$shannon,
    simpson_before = before$simpson, simpson_after = after$simpson,
    simpson_change = after$simpson - before$simpson,
    stringsAsFactors = FALSE)
}
paired <- do.call(rbind, paired_rows)
write.table(paired, file.path(result_dir, "alpha", "MA_paired_changes_all_vs_fungal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

greenhouse_summary <- aggregate(cbind(observed_change, shannon_change, simpson_change) ~
                                   subset + greenhouse_id + continuous_cropping_years,
                                 data = paired, FUN = mean)
greenhouse_summary$n_spatial_pairs <- 3L
greenhouse_summary$interpretation <- "descriptive_greenhouse_mean_of_3_nested_spatial_pairs"
write.table(greenhouse_summary, file.path(result_dir, "alpha", "MA_greenhouse_mean_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

paired_summary <- do.call(rbind, lapply(split(paired, paired$subset), function(z) {
  data.frame(subset = unique(z$subset), n_spatial_pairs = nrow(z), n_greenhouses = 3L,
             observed_mean_change = mean(z$observed_change),
             observed_median_change = median(z$observed_change),
             observed_positive_pairs = sum(z$observed_change > 0),
             shannon_mean_change = mean(z$shannon_change),
             shannon_median_change = median(z$shannon_change),
             shannon_positive_pairs = sum(z$shannon_change > 0),
             simpson_mean_change = mean(z$simpson_change),
             stringsAsFactors = FALSE)
}))
write.table(paired_summary, file.path(result_dir, "alpha", "MA_paired_change_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hd <- alpha[alpha$cohort == "white_mold_case_example", ]
hd_summary <- do.call(rbind, lapply(split(hd, hd$subset), function(z) {
  h <- z[z$condition == "healthy", ]
  d <- z[z$condition == "diseased", ]
  data.frame(subset = unique(z$subset), healthy_n_subsamples = nrow(h),
             diseased_n_subsamples = nrow(d),
             observed_healthy_mean = mean(h$observed_asv),
             observed_diseased_mean = mean(d$observed_asv),
             observed_difference_D_minus_H = mean(d$observed_asv) - mean(h$observed_asv),
             shannon_healthy_mean = mean(h$shannon),
             shannon_diseased_mean = mean(d$shannon),
             shannon_difference_D_minus_H = mean(d$shannon) - mean(h$shannon),
             interpretation = "descriptive_two_points_in_one_greenhouse",
             stringsAsFactors = FALSE)
}))
write.table(hd_summary, file.path(result_dir, "alpha", "HD_descriptive_alpha.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bray_one <- function(x) vegdist(x / rowSums(x), method = "bray")
bray_all <- bray_one(all_counts)
bray_fungal <- bray_one(fungal_counts)
distance_sensitivity <- data.frame(
  comparison = "all_ASVs_vs_UNITE_kingdom_Fungi",
  n_samples = nrow(all_counts),
  pearson_distance_correlation = cor(as.vector(bray_all), as.vector(bray_fungal), method = "pearson"),
  spearman_distance_correlation = cor(as.vector(bray_all), as.vector(bray_fungal), method = "spearman"),
  stringsAsFactors = FALSE)
write.table(distance_sensitivity, file.path(result_dir, "beta", "bray_distance_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pcoa_one <- function(d, subset_name) {
  fit <- cmdscale(d, k = 2, eig = TRUE, add = TRUE)
  explained <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  data.frame(sample_id = rownames(fit$points), subset = subset_name,
             axis1 = fit$points[, 1], axis2 = fit$points[, 2],
             axis1_percent = explained[1], axis2_percent = explained[2],
             stringsAsFactors = FALSE)
}
pcoa <- rbind(pcoa_one(bray_all, "all_ASVs"),
              pcoa_one(bray_fungal, "UNITE_kingdom_Fungi"))
pcoa <- merge(pcoa, meta[, c("sample_id", "cohort", "condition", "greenhouse_id",
                             "pair_id", "independent_unit")], by = "sample_id", sort = FALSE)
write.table(pcoa, file.path(result_dir, "beta", "bray_pcoa_all_vs_fungal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

permanova_rows <- list()
for (ss in c("all_ASVs", "UNITE_kingdom_Fungi")) {
  d <- if (ss == "all_ASVs") bray_all else bray_fungal
  ids <- meta$sample_id[meta$cohort == "cultivation_paired"]
  zmeta <- meta[match(ids, meta$sample_id), ]
  zmeta$condition <- factor(zmeta$condition,
                            levels = c("before_cultivation", "after_cultivation"))
  zd <- as.dist(as.matrix(d)[ids, ids])
  control <- how(nperm = 999, blocks = factor(zmeta$pair_id))
  fit <- adonis2(zd ~ condition, data = zmeta, permutations = control, by = "margin")
  permanova_rows[[ss]] <- data.frame(
    subset = ss, comparison = "M_vs_A_pair_blocked_exploratory",
    R2 = fit$R2[1], F = fit$F[1], p_value_exploratory = fit$`Pr(>F)`[1],
    permutation_constraint = "within_pair_id",
    limitation = "9_spatial_pairs_nested_in_3_greenhouses_year_confounded_with_greenhouse",
    stringsAsFactors = FALSE)
}
permanova <- do.call(rbind, permanova_rows)
write.table(permanova, file.path(result_dir, "beta", "MA_permanova_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hd_distance <- function(d, subset_name) {
  ids_h <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "healthy"]
  ids_d <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "diseased"]
  m <- as.matrix(d)
  within_h <- m[ids_h, ids_h][lower.tri(m[ids_h, ids_h])]
  within_d <- m[ids_d, ids_d][lower.tri(m[ids_d, ids_d])]
  between <- as.vector(m[ids_h, ids_d])
  data.frame(subset = subset_name, mean_within_healthy = mean(within_h),
             mean_within_diseased = mean(within_d), mean_between_points = mean(between),
             between_minus_pooled_within = mean(between) - mean(c(within_h, within_d)),
             interpretation = "descriptive_two_points_in_one_greenhouse", stringsAsFactors = FALSE)
}
hd_beta <- rbind(hd_distance(bray_all, "all_ASVs"),
                 hd_distance(bray_fungal, "UNITE_kingdom_Fungi"))
write.table(hd_beta, file.path(result_dir, "beta", "HD_descriptive_bray.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p_pair <- ggplot(paired, aes(x = factor(c(0, 1)), group = pair_id, color = greenhouse_id))
plot_pair_metric <- function(metric, ylab, filename) {
  before <- paired[, c("subset", "pair_id", "greenhouse_id", paste0(metric, "_before"))]
  after <- paired[, c("subset", "pair_id", "greenhouse_id", paste0(metric, "_after"))]
  names(before)[4] <- "value"; names(after)[4] <- "value"
  before$time <- "Before"; after$time <- "After"
  z <- rbind(before, after)
  z$time <- factor(z$time, levels = c("Before", "After"))
  p <- ggplot(z, aes(time, value, group = pair_id, color = greenhouse_id)) +
    geom_line(alpha = 0.65) + geom_point(size = 2) + facet_wrap(~subset, scales = "free_y") +
    theme_bw(base_size = 10) + labs(x = NULL, y = ylab, color = "Greenhouse")
  ggsave(file.path(result_dir, "figures", filename), p, width = 9, height = 5)
}
plot_pair_metric("observed", "Observed ASVs", "MA_paired_observed_all_vs_fungal.pdf")
plot_pair_metric("shannon", "Shannon diversity", "MA_paired_shannon_all_vs_fungal.pdf")

p_beta <- ggplot(pcoa, aes(axis1, axis2, color = condition, shape = greenhouse_id)) +
  geom_point(size = 2.8) + facet_wrap(~subset, scales = "free") +
  theme_bw(base_size = 10) + labs(x = "PCoA1", y = "PCoA2", color = "Condition", shape = "Greenhouse")
ggsave(file.path(result_dir, "figures", "bray_pcoa_all_vs_fungal.pdf"),
       p_beta, width = 10, height = 5)

qc <- c("ITS design-aware sensitivity analysis",
        paste0("Samples=", nrow(all_counts)),
        paste0("All_ASVs=", ncol(all_counts)),
        paste0("Fungal_ASVs=", ncol(fungal_counts)),
        "M/A: 9 paired spatial points nested in 3 greenhouses.",
        "H/D: two points in one greenhouse, 3 spatial subsamples per point.",
        "All p-values are exploratory; H/D has no independent-replicate test.")
writeLines(qc, file.path(result_dir, "qc", "design_constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
cat("ITS_DESIGN_AWARE_SENSITIVITY_V1_OK\n")

