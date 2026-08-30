#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(vegan)
  library(permute)
  library(ggplot2)
  library(scales)
})

set.seed(20260822L)
project <- "."
out16 <- file.path(project, "analysis", "downstream", "16S_v1_silva_diversity")
outjoint <- file.path(project, "analysis", "downstream", "joint_16S_ITS_v1")
if (dir.exists(out16) || dir.exists(outjoint)) stop("Refusing to overwrite versioned output directories")
for (d in file.path(out16, c("objects", "qc", "alpha", "beta", "taxonomy", "figures")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (d in file.path(outjoint, c("tables", "figures", "qc")))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

count_file <- file.path(project, "analysis", "hpc_results", "16S_dada2_v1_corrected_v2",
                        "results", "asv_count_table.tsv")
tax_file <- file.path(project, "analysis", "hpc_results", "16S_silva1382_v2",
                      "asv_taxonomy_silva1382.tsv")
flag_file <- file.path(project, "analysis", "hpc_results", "16S_silva1382_v2",
                       "asv_taxonomic_cleaning_flags.tsv")
meta_file <- file.path(project, "analysis", "metadata", "sample_metadata.tsv")
its_frozen_root <- file.path(project, "analysis", "frozen", "ITS_STAGE1_FREEZE_2026-08-21_v1",
                             "evidence", "analysis")
its_count_file <- file.path(its_frozen_root, "downstream", "ITS_v2_taxonomy", "tables",
                            "fungal_asv_count_table.tsv")
its_tax_file <- file.path(its_frozen_root, "downstream", "ITS_v2_taxonomy", "tables",
                          "fungal_asv_taxonomy.tsv")
its_priority_file <- file.path(its_frozen_root, "downstream", "ITS_v4_genus_screening", "tables",
                               "MA_prioritized_genera.tsv")
required <- c(count_file, tax_file, flag_file, meta_file, its_count_file, its_tax_file,
              its_priority_file)
if (!all(file.exists(required))) stop("Missing required input: ", paste(required[!file.exists(required)], collapse = "; "))

read_asv_counts <- function(path) {
  x <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"ASV_ID" %in% names(x) || anyDuplicated(x$ASV_ID)) stop("Invalid ASV table: ", path)
  m <- t(as.matrix(x[, setdiff(names(x), "ASV_ID"), drop = FALSE]))
  storage.mode(m) <- "numeric"
  colnames(m) <- x$ASV_ID
  if (any(!is.finite(m)) || any(m < 0) || any(m != round(m))) stop("Non-integer counts: ", path)
  m
}

counts_all <- read_asv_counts(count_file)
tax <- read.delim(tax_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
flags <- read.delim(flag_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!setequal(colnames(counts_all), tax$ASV_ID) || !setequal(colnames(counts_all), flags$ASV_ID))
  stop("16S ASV IDs differ between counts, taxonomy and cleaning flags")
tax <- tax[match(colnames(counts_all), tax$ASV_ID), ]
flags <- flags[match(colnames(counts_all), flags$ASV_ID), ]
keep <- flags$proposed_keep_for_prokaryote_analysis
if (anyNA(keep) || !is.logical(keep)) stop("Invalid proposed_keep flag")
counts <- counts_all[, keep, drop = FALSE]
tax_clean <- tax[keep, , drop = FALSE]
flags_clean <- flags[keep, , drop = FALSE]
if (!identical(colnames(counts), tax_clean$ASV_ID)) stop("Clean table order mismatch")

meta_all <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta_all[meta_all$marker == "16S", ]
if (nrow(meta) != 24L || anyDuplicated(meta$biological_sample_id) ||
    !setequal(rownames(counts), meta$biological_sample_id)) stop("16S metadata mismatch")
meta <- meta[match(rownames(counts), meta$biological_sample_id), ]
meta$sample_id <- meta$biological_sample_id
meta$condition_short <- ifelse(meta$condition == "before_cultivation", "M_before",
                               ifelse(meta$condition == "after_cultivation", "A_after",
                                      ifelse(meta$condition == "healthy", "H_healthy", "D_diseased")))

depth_before <- rowSums(counts_all)
depth_after <- rowSums(counts)
sample_qc <- data.frame(
  sample_id = rownames(counts), cohort = meta$cohort, condition = meta$condition,
  greenhouse_id = meta$greenhouse_id, pair_id = meta$pair_id,
  reads_before_organelle_filter = depth_before,
  reads_after_organelle_filter = depth_after,
  reads_removed = depth_before - depth_after,
  removed_fraction = (depth_before - depth_after) / depth_before,
  observed_asvs_after_filter = specnumber(counts),
  stringsAsFactors = FALSE)
write.table(sample_qc, file.path(out16, "qc", "sample_qc_after_organelle_filter.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

filter_summary <- data.frame(
  metric = c("input_asvs", "retained_asvs", "removed_asvs_union", "chloroplast_flagged_asvs",
             "mitochondria_flagged_asvs", "input_reads", "retained_reads", "removed_reads_union",
             "retained_read_fraction"),
  value = c(ncol(counts_all), ncol(counts), sum(!keep), sum(flags$is_chloroplast),
            sum(flags$is_mitochondria), sum(counts_all), sum(counts),
            sum(counts_all) - sum(counts), sum(counts) / sum(counts_all)),
  stringsAsFactors = FALSE)
write.table(filter_summary, file.path(out16, "qc", "organelle_filter_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(ASV_ID = colnames(counts), t(counts), check.names = FALSE),
            file.path(out16, "objects", "16S_clean_asv_count_table.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(tax_clean, file.path(out16, "objects", "16S_clean_taxonomy_silva1382.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(out16, "objects", "16S_sample_metadata.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

depths <- c(1000L, 1500L, 1800L, 2500L, 3000L, 5000L)
depth_eval <- do.call(rbind, lapply(depths, function(d) {
  ok <- depth_after >= d
  data.frame(depth = d, retained_samples = sum(ok), dropped_samples = sum(!ok),
             dropped_sample_ids = paste(rownames(counts)[!ok], collapse = ";"),
             min_observed_depth = min(depth_after), stringsAsFactors = FALSE)
}))
write.table(depth_eval, file.path(out16, "qc", "rarefaction_depth_retention.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

primary_depth <- 1500L
set.seed(20260822L)
rare <- rrarefy(counts, sample = primary_depth)
if (!all(rowSums(rare) == primary_depth)) stop("Primary rarefaction did not retain all samples")
write.table(data.frame(sample_id = rownames(rare), rare, check.names = FALSE),
            file.path(out16, "objects", "16S_rarefied_counts_depth1500.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
rel <- counts / rowSums(counts)

alpha_metrics <- function(x, depth_label) data.frame(
  sample_id = rownames(x), rarefaction_depth = depth_label, observed_asv = specnumber(x),
  shannon = diversity(x, "shannon"), shannon_effective = exp(diversity(x, "shannon")),
  invsimpson = diversity(x, "invsimpson"),
  pielou_evenness = ifelse(specnumber(x) > 1, diversity(x, "shannon") / log(specnumber(x)), NA_real_),
  stringsAsFactors = FALSE)
alpha <- alpha_metrics(rare, primary_depth)
alpha <- merge(alpha, meta[, c("sample_id", "cohort", "condition", "condition_short",
                               "continuous_cropping_years", "greenhouse_id", "spatial_point",
                               "pair_id", "independent_unit")], by = "sample_id", sort = FALSE)
alpha <- alpha[match(rownames(rare), alpha$sample_id), ]
write.table(alpha, file.path(out16, "alpha", "alpha_diversity_depth1500.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

alpha_sensitivity <- list()
for (d in c(1500L, 1800L, 2500L)) {
  ids <- rownames(counts)[depth_after >= d]
  set.seed(20260822L + d)
  rd <- rrarefy(counts[ids, , drop = FALSE], sample = d)
  z <- alpha_metrics(rd, d)
  z$cohort <- meta$cohort[match(z$sample_id, meta$sample_id)]
  z$condition <- meta$condition[match(z$sample_id, meta$sample_id)]
  alpha_sensitivity[[as.character(d)]] <- z
}
write.table(do.call(rbind, alpha_sensitivity),
            file.path(out16, "alpha", "alpha_rarefaction_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

ma <- alpha[alpha$cohort == "cultivation_paired", ]
before <- ma[ma$condition == "before_cultivation", ]
after <- ma[ma$condition == "after_cultivation", ]
before <- before[match(after$pair_id, before$pair_id), ]
if (nrow(after) != 9L || !identical(before$pair_id, after$pair_id)) stop("M/A pairing mismatch")
paired <- data.frame(
  pair_id = after$pair_id, greenhouse_id = after$greenhouse_id,
  continuous_cropping_years = after$continuous_cropping_years,
  before_sample = before$sample_id, after_sample = after$sample_id,
  observed_before = before$observed_asv, observed_after = after$observed_asv,
  observed_change = after$observed_asv - before$observed_asv,
  shannon_before = before$shannon, shannon_after = after$shannon,
  shannon_change = after$shannon - before$shannon,
  shannon_effective_before = before$shannon_effective,
  shannon_effective_after = after$shannon_effective,
  shannon_effective_change = after$shannon_effective - before$shannon_effective,
  invsimpson_change = after$invsimpson - before$invsimpson,
  stringsAsFactors = FALSE)
write.table(paired, file.path(out16, "alpha", "MA_paired_alpha_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
gh_alpha <- aggregate(cbind(observed_change, shannon_change, shannon_effective_change,
                            invsimpson_change) ~ greenhouse_id + continuous_cropping_years,
                      data = paired, FUN = mean)
gh_alpha$n_nested_spatial_pairs <- 3L
gh_alpha$interpretation <- "descriptive_greenhouse_mean_year_confounded_with_greenhouse"
write.table(gh_alpha, file.path(out16, "alpha", "MA_greenhouse_mean_alpha_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
alpha_summary <- data.frame(
  metric = c("observed_asv", "shannon", "shannon_effective", "invsimpson"),
  mean_change_A_minus_M = c(mean(paired$observed_change), mean(paired$shannon_change),
                            mean(paired$shannon_effective_change), mean(paired$invsimpson_change)),
  median_change_A_minus_M = c(median(paired$observed_change), median(paired$shannon_change),
                              median(paired$shannon_effective_change), median(paired$invsimpson_change)),
  positive_pairs = c(sum(paired$observed_change > 0), sum(paired$shannon_change > 0),
                     sum(paired$shannon_effective_change > 0), sum(paired$invsimpson_change > 0)),
  negative_pairs = c(sum(paired$observed_change < 0), sum(paired$shannon_change < 0),
                     sum(paired$shannon_effective_change < 0), sum(paired$invsimpson_change < 0)),
  stringsAsFactors = FALSE)
write.table(alpha_summary, file.path(out16, "alpha", "MA_alpha_change_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hd <- alpha[alpha$cohort == "white_mold_case_example", ]
hd_alpha <- do.call(rbind, lapply(c("observed_asv", "shannon", "shannon_effective", "invsimpson"), function(v) {
  h <- hd[hd$condition == "healthy", v]; d <- hd[hd$condition == "diseased", v]
  data.frame(metric = v, healthy_point_mean = mean(h), diseased_point_mean = mean(d),
             D_minus_H = mean(d) - mean(h), healthy_subsamples = length(h),
             diseased_subsamples = length(d),
             interpretation = "descriptive_two_points_in_one_greenhouse_no_independent_replicate",
             stringsAsFactors = FALSE)
}))
write.table(hd_alpha, file.path(out16, "alpha", "HD_descriptive_alpha.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bray <- vegdist(rare, method = "bray")
jaccard <- vegdist(rare, method = "jaccard", binary = TRUE)
write.table(as.matrix(bray), file.path(out16, "beta", "bray_distance_matrix_depth1500.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)
write.table(as.matrix(jaccard), file.path(out16, "beta", "jaccard_distance_matrix_depth1500.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

pcoa_df <- function(d, metric) {
  fit <- cmdscale(d, eig = TRUE, add = TRUE, k = 2)
  pct <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
  data.frame(sample_id = rownames(fit$points), metric = metric,
             axis1 = fit$points[, 1], axis2 = fit$points[, 2],
             axis1_percent = pct[1], axis2_percent = pct[2], stringsAsFactors = FALSE)
}
pcoa <- rbind(pcoa_df(bray, "Bray_Curtis"), pcoa_df(jaccard, "Jaccard_binary"))
pcoa <- merge(pcoa, meta[, c("sample_id", "cohort", "condition", "condition_short",
                             "greenhouse_id", "pair_id")], by = "sample_id", sort = FALSE)
write.table(pcoa, file.path(out16, "beta", "pcoa_coordinates.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

beta_tests <- list()
for (metric in c("Bray_Curtis", "Jaccard_binary")) {
  d <- if (metric == "Bray_Curtis") bray else jaccard
  ids <- meta$sample_id[meta$cohort == "cultivation_paired"]
  zm <- meta[match(ids, meta$sample_id), ]
  zm$condition <- factor(zm$condition, levels = c("before_cultivation", "after_cultivation"))
  zd <- as.dist(as.matrix(d)[ids, ids])
  ctrl <- how(nperm = 9999, blocks = factor(zm$pair_id))
  fit <- adonis2(zd ~ condition, data = zm, permutations = ctrl)
  disp <- betadisper(zd, zm$condition)
  dfit <- permutest(disp, permutations = 9999)
  beta_tests[[metric]] <- data.frame(
    metric = metric, comparison = "M_vs_A_pair_blocked_exploratory", n = length(ids),
    R2 = fit$R2[1], pseudo_F = fit$F[1], permanova_p = fit$`Pr(>F)`[1],
    dispersion_F = dfit$tab[1, "F"], dispersion_p = dfit$tab[1, "Pr(>F)"],
    permutation_constraint = "within_pair_id_9999",
    limitation = "9_spatial_pairs_nested_in_3_greenhouses;year_fully_confounded_with_greenhouse",
    stringsAsFactors = FALSE)
}
write.table(do.call(rbind, beta_tests), file.path(out16, "beta", "MA_permanova_and_dispersion.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

hd_beta_one <- function(d, metric) {
  h <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "healthy"]
  z <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "diseased"]
  m <- as.matrix(d); wh <- m[h, h][lower.tri(m[h, h])]; wd <- m[z, z][lower.tri(m[z, z])]
  bt <- as.vector(m[h, z])
  data.frame(metric = metric, mean_within_healthy = mean(wh), mean_within_diseased = mean(wd),
             mean_between_points = mean(bt), between_minus_pooled_within = mean(bt) - mean(c(wh, wd)),
             interpretation = "descriptive_two_points_in_one_greenhouse_no_independent_test",
             stringsAsFactors = FALSE)
}
write.table(rbind(hd_beta_one(bray, "Bray_Curtis"), hd_beta_one(jaccard, "Jaccard_binary")),
            file.path(out16, "beta", "HD_descriptive_beta.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

clean_label <- function(x, fallback) {
  x[is.na(x) | !nzchar(x)] <- fallback
  x
}
aggregate_rank <- function(rank, fallback) {
  lab <- clean_label(tax_clean[[rank]], fallback)
  t(rowsum(t(counts), group = lab, reorder = TRUE))
}
phylum_counts <- aggregate_rank("Phylum", "Unclassified_at_phylum")
genus_counts <- aggregate_rank("Genus", "Unclassified_at_genus")
phylum_rel <- phylum_counts / rowSums(phylum_counts)
genus_rel <- genus_counts / rowSums(genus_counts)
write.table(data.frame(sample_id = rownames(phylum_rel), phylum_rel, check.names = FALSE),
            file.path(out16, "taxonomy", "phylum_relative_abundance.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(data.frame(sample_id = rownames(genus_rel), genus_rel, check.names = FALSE),
            file.path(out16, "taxonomy", "genus_relative_abundance.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

rank_summary <- function(m, rank) {
  z <- sort(colSums(m), decreasing = TRUE)
  data.frame(rank = rank, taxon = names(z), total_reads = unname(z),
             global_relative_abundance = unname(z / sum(z)),
             prevalence_samples = colSums(m > 0)[names(z)], stringsAsFactors = FALSE)
}
write.table(rbind(rank_summary(phylum_counts, "Phylum"), rank_summary(genus_counts, "Genus")),
            file.path(out16, "taxonomy", "taxon_global_summary.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

genus_clr <- log2(genus_counts + 0.5)
genus_clr <- genus_clr - rowMeans(genus_clr)
b_ids <- before$sample_id; a_ids <- after$sample_id
ma_genus <- do.call(rbind, lapply(seq_len(ncol(genus_counts)), function(j) {
  delta <- genus_clr[a_ids, j] - genus_clr[b_ids, j]
  gh <- aggregate(delta, list(greenhouse_id = after$greenhouse_id), mean)$x
  p <- tryCatch(wilcox.test(genus_clr[a_ids, j], genus_clr[b_ids, j], paired = TRUE,
                            exact = FALSE)$p.value, error = function(e) NA_real_)
  data.frame(genus = colnames(genus_counts)[j], total_reads = sum(genus_counts[, j]),
             prevalence_MA = sum(genus_counts[c(b_ids, a_ids), j] > 0),
             mean_relative_M = mean(genus_rel[b_ids, j]), mean_relative_A = mean(genus_rel[a_ids, j]),
             log2_mean_relative_ratio_A_vs_M = log2((mean(genus_rel[a_ids, j]) + 1e-6) /
                                                      (mean(genus_rel[b_ids, j]) + 1e-6)),
             median_paired_clr_change = median(delta), mean_paired_clr_change = mean(delta),
             positive_pairs = sum(delta > 0), negative_pairs = sum(delta < 0),
             pair_direction_consistency = max(sum(delta > 0), sum(delta < 0)) / 9,
             greenhouse_direction_consistent = all(gh > 0) || all(gh < 0),
             wilcoxon_p_exploratory = p, stringsAsFactors = FALSE)
}))
ma_genus$BH_q_exploratory <- p.adjust(ma_genus$wilcoxon_p_exploratory, "BH")
eligible <- ma_genus$genus != "Unclassified_at_genus" & ma_genus$prevalence_MA >= 6 &
  ma_genus$total_reads >= 100
ma_genus$priority <- "not_prioritized"
moderate <- eligible & abs(ma_genus$median_paired_clr_change) >= 0.5 &
  abs(ma_genus$log2_mean_relative_ratio_A_vs_M) >= 0.5 &
  sign(ma_genus$median_paired_clr_change) == sign(ma_genus$log2_mean_relative_ratio_A_vs_M) &
  ma_genus$pair_direction_consistency >= 6/9 & ma_genus$greenhouse_direction_consistent
high <- eligible & abs(ma_genus$median_paired_clr_change) >= 1 &
  abs(ma_genus$log2_mean_relative_ratio_A_vs_M) >= 1 &
  sign(ma_genus$median_paired_clr_change) == sign(ma_genus$log2_mean_relative_ratio_A_vs_M) &
  ma_genus$pair_direction_consistency >= 7/9 & ma_genus$greenhouse_direction_consistent
ma_genus$priority[moderate] <- "moderate"
ma_genus$priority[high] <- "high"
ma_genus <- ma_genus[order(factor(ma_genus$priority, c("high", "moderate", "not_prioritized")),
                           -abs(ma_genus$median_paired_clr_change)), ]
write.table(ma_genus, file.path(out16, "taxonomy", "MA_genus_screening.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(ma_genus[ma_genus$priority != "not_prioritized", ],
            file.path(out16, "taxonomy", "MA_prioritized_genera.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

h_ids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "healthy"]
d_ids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "diseased"]
hd_genus <- do.call(rbind, lapply(seq_len(ncol(genus_counts)), function(j) {
  hm <- mean(genus_rel[h_ids, j]); dm <- mean(genus_rel[d_ids, j])
  data.frame(genus = colnames(genus_counts)[j], total_reads = sum(genus_counts[, j]),
             detected_H = sum(genus_counts[h_ids, j] > 0), detected_D = sum(genus_counts[d_ids, j] > 0),
             mean_relative_H = hm, mean_relative_D = dm,
             log2_mean_relative_ratio_D_vs_H = log2((dm + 1e-6) / (hm + 1e-6)),
             interpretation = "descriptive_two_points_one_greenhouse", stringsAsFactors = FALSE)
}))
hd_genus <- hd_genus[order(-abs(hd_genus$log2_mean_relative_ratio_D_vs_H)), ]
write.table(hd_genus, file.path(out16, "taxonomy", "HD_descriptive_genus.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

analysis_object <- list(counts_clean = counts, counts_rarefied_depth1500 = rare,
                        taxonomy = tax_clean, flags = flags_clean, metadata = meta,
                        relative_abundance = rel, primary_rarefaction_depth = primary_depth,
                        seed = 20260822L)
saveRDS(analysis_object, file.path(out16, "objects", "16S_analysis_object_v1.rds"))

theme_pub <- theme_bw(base_size = 10) +
  theme(panel.grid.minor = element_blank(), legend.position = "right",
        axis.text.x = element_text(angle = 45, hjust = 1))
cond_cols <- c(M_before = "#3B6FB6", A_after = "#E07A2D", H_healthy = "#2A9D8F", D_diseased = "#C43C39")

p_depth <- ggplot(sample_qc, aes(factor(sample_id, levels = sample_id[order(reads_after_organelle_filter)]),
                                 reads_after_organelle_filter, fill = condition)) +
  geom_col() + geom_hline(yintercept = primary_depth, linetype = 2) +
  scale_y_continuous(labels = comma) + labs(x = NULL, y = "Reads after organelle filtering",
                                            title = "16S sample depth after organelle filtering") + theme_pub
ggsave(file.path(out16, "figures", "16S_sample_depth_after_filter.pdf"), p_depth, width = 9, height = 5)

pair_plot <- function(metric, ylabel, file) {
  b <- paired[, c("pair_id", "greenhouse_id", paste0(metric, "_before"))]
  a <- paired[, c("pair_id", "greenhouse_id", paste0(metric, "_after"))]
  names(b)[3] <- "value"; names(a)[3] <- "value"; b$time <- "M_before"; a$time <- "A_after"
  z <- rbind(b, a); z$time <- factor(z$time, c("M_before", "A_after"))
  p <- ggplot(z, aes(time, value, group = pair_id, color = greenhouse_id)) +
    geom_line(alpha = 0.7) + geom_point(size = 2.2) + theme_pub +
    labs(x = NULL, y = ylabel, color = "Greenhouse")
  ggsave(file.path(out16, "figures", file), p, width = 6.5, height = 5)
}
pair_plot("observed", "Observed ASVs at 1,500 reads", "MA_paired_observed_asv.pdf")
pair_plot("shannon", "Shannon diversity (ln)", "MA_paired_shannon.pdf")

p_b <- ggplot(pcoa[pcoa$metric == "Bray_Curtis", ],
              aes(axis1, axis2, color = condition_short, shape = greenhouse_id)) +
  geom_point(size = 3) + scale_color_manual(values = cond_cols) + theme_pub +
  labs(x = "Bray-Curtis PCoA1", y = "Bray-Curtis PCoA2", color = "Condition", shape = "Greenhouse")
ggsave(file.path(out16, "figures", "16S_bray_pcoa.pdf"), p_b, width = 8, height = 6)

long_top <- function(m, top_n, rank) {
  keep_taxa <- names(head(sort(colMeans(m), decreasing = TRUE), top_n))
  mm <- m[, keep_taxa, drop = FALSE]
  other <- pmax(0, 1 - rowSums(mm))
  mm <- cbind(mm, Other = other)
  z <- as.data.frame(as.table(mm), stringsAsFactors = FALSE)
  names(z) <- c("sample_id", "taxon", "relative_abundance")
  z$rank <- rank
  z$condition_short <- meta$condition_short[match(z$sample_id, meta$sample_id)]
  z
}
phylum_long <- long_top(phylum_rel, 10L, "Phylum")
genus_long <- long_top(genus_rel, 15L, "Genus")
write.table(phylum_long, file.path(out16, "taxonomy", "phylum_top10_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(genus_long, file.path(out16, "taxonomy", "genus_top15_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
plot_comp <- function(z, file, title) {
  p <- ggplot(z, aes(sample_id, relative_abundance, fill = taxon)) + geom_col(width = 0.9) +
    facet_grid(~condition_short, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = percent_format()) + theme_pub +
    theme(axis.text.x = element_text(angle = 70, hjust = 1, size = 7)) +
    labs(x = NULL, y = "Relative abundance", fill = NULL, title = title)
  ggsave(file.path(out16, "figures", file), p, width = 13, height = 6)
}
plot_comp(phylum_long, "16S_phylum_composition_top10.pdf", "16S phylum composition")
plot_comp(genus_long, "16S_genus_composition_top15.pdf", "16S genus composition")

rare_grid <- unique(round(seq(100, min(10000, max(depth_after)), length.out = 20)))
expected_richness <- function(v, n) sum(1 - exp(lchoose(sum(v) - v, n) - lchoose(sum(v), n)), na.rm = TRUE)
rare_curves <- do.call(rbind, lapply(rownames(counts), function(id) {
  ds <- rare_grid[rare_grid <= depth_after[id]]
  data.frame(sample_id = id, depth = ds,
             expected_asvs = vapply(ds, function(n) expected_richness(counts[id, ], n), numeric(1)),
             condition_short = meta$condition_short[match(id, meta$sample_id)], stringsAsFactors = FALSE)
}))
write.table(rare_curves, file.path(out16, "qc", "alpha_rarefaction_curves.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
p_r <- ggplot(rare_curves, aes(depth, expected_asvs, group = sample_id, color = condition_short)) +
  geom_line(alpha = 0.65) + geom_vline(xintercept = primary_depth, linetype = 2) +
  scale_color_manual(values = cond_cols) + theme_pub +
  labs(x = "Rarefaction depth", y = "Expected observed ASVs", color = "Condition")
ggsave(file.path(out16, "figures", "16S_alpha_rarefaction_curves.pdf"), p_r, width = 8, height = 6)

# Cross-marker integration: each marker is normalized independently; raw counts are never concatenated.
its_counts <- read_asv_counts(its_count_file)
its_tax <- read.delim(its_tax_file, check.names = FALSE, stringsAsFactors = FALSE, na.strings = character())
if (!setequal(rownames(its_counts), meta$sample_id)) stop("ITS frozen sample IDs differ from 16S")
its_counts <- its_counts[meta$sample_id, , drop = FALSE]
set.seed(20260823L)
its_rare <- rrarefy(its_counts, sample = 16000L)
its_alpha <- alpha_metrics(its_rare, 16000L)
joint_alpha <- merge(alpha[, c("sample_id", "cohort", "condition", "greenhouse_id", "pair_id",
                               "observed_asv", "shannon", "shannon_effective", "invsimpson")],
                     its_alpha[, c("sample_id", "observed_asv", "shannon", "shannon_effective", "invsimpson")],
                     by = "sample_id", suffixes = c("_16S", "_ITS"), sort = FALSE)
joint_alpha <- joint_alpha[match(meta$sample_id, joint_alpha$sample_id), ]
write.table(joint_alpha, file.path(outjoint, "tables", "joint_alpha_by_sample.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

joint_ma <- joint_alpha[joint_alpha$cohort == "cultivation_paired", ]
jb <- joint_ma[joint_ma$condition == "before_cultivation", ]
ja <- joint_ma[joint_ma$condition == "after_cultivation", ]
jb <- jb[match(ja$pair_id, jb$pair_id), ]
joint_delta <- data.frame(pair_id = ja$pair_id, greenhouse_id = ja$greenhouse_id,
                          observed_delta_16S = ja$observed_asv_16S - jb$observed_asv_16S,
                          observed_delta_ITS = ja$observed_asv_ITS - jb$observed_asv_ITS,
                          shannon_delta_16S = ja$shannon_16S - jb$shannon_16S,
                          shannon_delta_ITS = ja$shannon_ITS - jb$shannon_ITS,
                          stringsAsFactors = FALSE)
write.table(joint_delta, file.path(outjoint, "tables", "MA_cross_marker_alpha_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cor_rows <- do.call(rbind, lapply(c("observed", "shannon"), function(v) {
  x <- joint_delta[[paste0(v, "_delta_16S")]]; y <- joint_delta[[paste0(v, "_delta_ITS")]]
  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))
  data.frame(metric = v, n_pairs = length(x), spearman_rho = unname(ct$estimate),
             p_exploratory = ct$p.value,
             limitation = "9_spatial_pairs_nested_in_3_greenhouses_no_independent_year_effect",
             stringsAsFactors = FALSE)
}))
write.table(cor_rows, file.path(outjoint, "tables", "MA_cross_marker_alpha_correlation.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

bray16 <- vegdist(rare, "bray")
brayits <- vegdist(its_rare, "bray")
mantel_rows <- list()
for (scope in c("all24", "MA18", "HD6")) {
  ids <- if (scope == "all24") meta$sample_id else if (scope == "MA18")
    meta$sample_id[meta$cohort == "cultivation_paired"] else meta$sample_id[meta$cohort == "white_mold_case_example"]
  d16 <- as.dist(as.matrix(bray16)[ids, ids]); dit <- as.dist(as.matrix(brayits)[ids, ids])
  mt <- mantel(d16, dit, method = "spearman", permutations = 9999)
  mantel_rows[[scope]] <- data.frame(scope = scope, n_samples = length(ids),
                                     spearman_mantel_r = unname(mt$statistic),
                                     p_exploratory = mt$signif,
                                     limitation = ifelse(scope == "HD6", "two_points_one_greenhouse",
                                                         "nested_spatial_subsamples_and_greenhouse_confounding"),
                                     stringsAsFactors = FALSE)
}
write.table(do.call(rbind, mantel_rows), file.path(outjoint, "tables", "cross_marker_bray_mantel.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

joint_pcoa <- rbind(transform(pcoa_df(bray16, "16S_Bray"), marker = "16S"),
                    transform(pcoa_df(brayits, "ITS_Bray"), marker = "ITS"))
joint_pcoa <- merge(joint_pcoa, meta[, c("sample_id", "cohort", "condition_short", "greenhouse_id")],
                    by = "sample_id", sort = FALSE)
write.table(joint_pcoa, file.path(outjoint, "tables", "cross_marker_pcoa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
p_joint <- ggplot(joint_pcoa, aes(axis1, axis2, color = condition_short, shape = greenhouse_id)) +
  geom_point(size = 2.7) + facet_wrap(~marker, scales = "free") +
  scale_color_manual(values = cond_cols) + theme_pub +
  labs(x = "PCoA1", y = "PCoA2", color = "Condition", shape = "Greenhouse")
ggsave(file.path(outjoint, "figures", "16S_ITS_bray_pcoa_side_by_side.pdf"), p_joint, width = 11, height = 5)

p_delta <- ggplot(joint_delta, aes(shannon_delta_16S, shannon_delta_ITS, color = greenhouse_id)) +
  geom_hline(yintercept = 0, color = "grey70") + geom_vline(xintercept = 0, color = "grey70") +
  geom_point(size = 3) + geom_text(aes(label = pair_id), vjust = -0.7, size = 3) + theme_pub +
  labs(x = "16S Shannon change (A - M)", y = "ITS Shannon change (A - M)", color = "Greenhouse")
ggsave(file.path(outjoint, "figures", "MA_cross_marker_shannon_changes.pdf"), p_delta, width = 7, height = 6)

its_priority <- read.delim(its_priority_file, check.names = FALSE, stringsAsFactors = FALSE)
write.table(its_priority, file.path(outjoint, "tables", "ITS_frozen_MA_prioritized_genera_snapshot.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(c(
  "# 16S analysis methods and design constraints", "",
  "- Input: normalized DADA2 non-chimeric counts and SILVA 138.2 taxonomy.",
  "- Chloroplast and mitochondria: removed by proposed_keep_for_prokaryote_analysis; source tables unchanged.",
  "- Alpha/beta primary branch: rarefied without replacement to 1,500 reads (seed 20260822), retaining 24/24 samples; sensitivity at 1,800 and 2,500 reads.",
  "- Alpha: Observed ASV, Shannon (natural log), exp(Shannon), Inverse Simpson and Pielou evenness.",
  "- Beta: Bray-Curtis and binary Jaccard; M/A PERMANOVA constrained within pair_id with dispersion tests.",
  "- Composition: relative abundance from unrarefied cleaned counts; genus screening uses CLR(+0.5).",
  "- Nine M/A spatial pairs are nested in three greenhouses; year is fully confounded with greenhouse.",
  "- H/D are two points in one greenhouse, with three subsamples per point; descriptive only.",
  "- No reliable reference placement tree was built, so Faith PD and UniFrac are not reported."
), file.path(out16, "METHODS_ZH.md"))

top_phyla <- head(sort(colSums(phylum_counts), decreasing = TRUE), 8)
top_genera <- head(sort(colSums(genus_counts), decreasing = TRUE), 10)
bt <- do.call(rbind, beta_tests)
report <- c(
  "# Morchella soil 16S stage report", "", "## Cleaning and QC", "",
  sprintf("- Input: %s ASVs and %s reads.", comma(ncol(counts_all)), comma(sum(counts_all))),
  sprintf("- After organelle removal: %s ASVs and %s reads (%.2f%% retained).", comma(ncol(counts)), comma(sum(counts)), 100 * sum(counts) / sum(counts_all)),
  sprintf("- Clean depth range: %s-%s reads; primary depth 1,500 retains 24/24 samples.", comma(min(depth_after)), comma(max(depth_after))),
  "", "## Alpha diversity", "",
  sprintf("- M to A: mean Observed change %.1f; mean Shannon change %.3f; Shannon increased in %d and decreased in %d of 9 pairs.",
          mean(paired$observed_change), mean(paired$shannon_change), sum(paired$shannon_change > 0), sum(paired$shannon_change < 0)),
  sprintf("- H/D descriptive Shannon point difference (D-H): %.3f; not an independent-replicate test.", hd_alpha$D_minus_H[hd_alpha$metric == "shannon"]),
  "", "## Beta diversity", "",
  sprintf("- M/A Bray-Curtis: R2=%.3f, paired PERMANOVA p=%.4f, dispersion p=%.4f.", bt$R2[bt$metric == "Bray_Curtis"], bt$permanova_p[bt$metric == "Bray_Curtis"], bt$dispersion_p[bt$metric == "Bray_Curtis"]),
  sprintf("- M/A Jaccard: R2=%.3f, paired PERMANOVA p=%.4f, dispersion p=%.4f.", bt$R2[bt$metric == "Jaccard_binary"], bt$permanova_p[bt$metric == "Jaccard_binary"], bt$dispersion_p[bt$metric == "Jaccard_binary"]),
  "", "## Composition", "",
  paste0("- Top phyla: ", paste(sprintf("%s (%.1f%%)", names(top_phyla), 100 * top_phyla / sum(phylum_counts)), collapse = "; "), "."),
  paste0("- Top genera: ", paste(sprintf("%s (%.1f%%)", names(top_genera), 100 * top_genera / sum(genus_counts)), collapse = "; "), "."),
  sprintf("- Genus screening: %d high and %d moderate exploratory candidates.", sum(ma_genus$priority == "high"), sum(ma_genus$priority == "moderate")),
  "", "## Limits", "",
  "- Three spatial subsamples are not three independent greenhouse replicates.",
  "- Cropping year cannot be separated from greenhouse.",
  "- H/D is a descriptive within-greenhouse two-point contrast."
)
writeLines(report, file.path(out16, "RUN_REPORT_ZH.md"))

mr <- do.call(rbind, mantel_rows)
joint_report <- c(
  "# 16S-ITS joint interpretation report", "", "## Strategy", "",
  "- 16S and ITS were cleaned, rarefied and analyzed separately; raw marker counts were not concatenated.",
  "- 16S used 1,500 reads and the frozen ITS fungal table used 16,000 reads; both retained 24 samples.",
  "- Joint evidence uses paired M/A alpha changes, Bray-Curtis Mantel correlations and side-by-side PCoA.",
  "", "## Results", "",
  sprintf("- Across 9 M/A spatial pairs, cross-marker Shannon-change Spearman rho=%.3f, exploratory p=%.4f.",
          cor_rows$spearman_rho[cor_rows$metric == "shannon"], cor_rows$p_exploratory[cor_rows$metric == "shannon"]),
  sprintf("- Across all 24 samples, cross-marker Bray-Curtis Mantel r=%.3f, exploratory p=%.4f.",
          mr$spearman_mantel_r[mr$scope == "all24"], mr$p_exploratory[mr$scope == "all24"]),
  sprintf("- M/A Mantel r=%.3f, p=%.4f; H/D Mantel r=%.3f, p=%.4f.",
          mr$spearman_mantel_r[mr$scope == "MA18"], mr$p_exploratory[mr$scope == "MA18"],
          mr$spearman_mantel_r[mr$scope == "HD6"], mr$p_exploratory[mr$scope == "HD6"]),
  sprintf("- Frozen ITS provides %d prioritized fungal genera; 16S provides %d prioritized bacterial/archaeal genera.",
          nrow(its_priority), sum(ma_genus$priority != "not_prioritized")),
  "", "## Interpretation limits", "",
  "- Concordance supports a candidate coordinated community response, not direct interaction or causality.",
  "- Discordance may reflect different trophic or niche responses.",
  "- Greenhouse-year confounding and the H/D design make all cross-marker statistics exploratory.",
  "- Yield/disease prediction is not possible until outcome data are available."
)
writeLines(joint_report, file.path(outjoint, "JOINT_INTERPRETATION_REPORT_ZH.md"))

writeLines(c("16S_primary_depth=1500", "ITS_joint_depth=16000", "seed_16S=20260822",
             "seed_ITS=20260823", "permutations=9999",
             "raw_marker_counts_concatenated=FALSE"), file.path(outjoint, "qc", "analysis_parameters.txt"))
writeLines(capture.output(sessionInfo()), file.path(out16, "sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), file.path(outjoint, "sessionInfo.txt"))
cat("MORCHELLA_16S_AND_JOINT_V1_OK\n")

