#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))

set.seed(20260821L)
project <- "."
base_name <- Sys.getenv("ITS_SOIL_BASE", unset = "ITS_v5_soil_integration")
base_dir <- file.path(project, "analysis", "downstream", base_name)
soil_candidates <- list.files(file.path(base_dir, "input"),
                              pattern = "^soil_metrics_normalized_v[0-9]+\\.tsv$",
                              full.names = TRUE)
if (length(soil_candidates) != 1L) stop("Expected exactly one normalized soil metric table")
soil_file <- soil_candidates[[1L]]
genus_count_file <- file.path(project, "analysis", "downstream", "ITS_v4_genus_screening",
                              "tables", "fungal_genus_count_table.tsv")
screen_file <- file.path(project, "analysis", "downstream", "ITS_v4_genus_screening",
                         "tables", "MA_genus_screening.tsv")
alpha_change_file <- file.path(project, "analysis", "downstream", "ITS_v3_design_aware",
                               "alpha", "MA_paired_changes_all_vs_fungal.tsv")
metadata_file <- file.path(project, "analysis", "metadata", "sample_metadata.tsv")
out_tables <- file.path(base_dir, "tables")
out_figures <- file.path(base_dir, "figures")
out_qc <- file.path(base_dir, "qc")
for (d in c(out_tables, out_figures, out_qc)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
targets <- file.path(out_tables, c("soil_paired_changes.tsv", "genus_soil_spearman.tsv",
                                  "genus_soil_screening_signals.tsv",
                                  "alpha_soil_spearman.tsv",
                                  "soil_greenhouse_mean_changes.tsv"))
if (any(file.exists(targets))) stop("Refusing to overwrite existing soil association outputs")

soil <- read.delim(soil_file, check.names = FALSE, stringsAsFactors = FALSE)
genus_tab <- read.delim(genus_count_file, check.names = FALSE, stringsAsFactors = FALSE)
screen <- read.delim(screen_file, check.names = FALSE, stringsAsFactors = FALSE)
alpha_change <- read.delim(alpha_change_file, check.names = FALSE, stringsAsFactors = FALSE)
meta_all <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta_all[meta_all$marker == "ITS" & meta_all$cohort == "cultivation_paired", ]
stopifnot(nrow(soil) == 18L, nrow(meta) == 18L, setequal(soil$sample_id, meta$biological_sample_id))
meta$sample_id <- meta$biological_sample_id
soil <- merge(soil, meta[, c("sample_id", "pair_id", "greenhouse_id",
                             "continuous_cropping_years", "condition")],
              by = "sample_id", sort = FALSE)

metric_names <- setdiff(names(soil), c("sample_id", "phase", "lab_sample_id", "pair_id",
                                       "greenhouse_id", "continuous_cropping_years", "condition"))
before <- soil[soil$condition == "before_cultivation", ]
after <- soil[soil$condition == "after_cultivation", ]
before <- before[match(after$pair_id, before$pair_id), ]
stopifnot(nrow(before) == 9L, identical(before$pair_id, after$pair_id))
soil_change <- data.frame(pair_id = after$pair_id, greenhouse_id = after$greenhouse_id,
                          continuous_cropping_years = after$continuous_cropping_years,
                          before_sample = before$sample_id, after_sample = after$sample_id,
                          stringsAsFactors = FALSE)
for (m in metric_names) {
  soil_change[[paste0(m, "_before")]] <- before[[m]]
  soil_change[[paste0(m, "_after")]] <- after[[m]]
  soil_change[[paste0(m, "_change")]] <- after[[m]] - before[[m]]
}
write.table(soil_change, file.path(out_tables, "soil_paired_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

change_cols <- paste0(metric_names, "_change")
gh_change <- aggregate(soil_change[, change_cols, drop = FALSE],
                       by = list(greenhouse_id = soil_change$greenhouse_id,
                                 continuous_cropping_years = soil_change$continuous_cropping_years),
                       FUN = mean)
gh_change$n_spatial_pairs <- 3L
gh_change$interpretation <- "descriptive_mean_of_nested_spatial_pairs"
write.table(gh_change, file.path(out_tables, "soil_greenhouse_mean_changes.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

genus_counts <- as.matrix(genus_tab[, -1L, drop = FALSE])
storage.mode(genus_counts) <- "numeric"
rownames(genus_counts) <- genus_tab$sample_id
genus_clr <- log2(genus_counts + 0.5)
genus_clr <- genus_clr - rowMeans(genus_clr)
b_ids <- before$sample_id; a_ids <- after$sample_id
stopifnot(all(c(b_ids, a_ids) %in% rownames(genus_clr)))
priority_genera <- screen$genus[screen$priority != "not_prioritized"]
priority_genera <- intersect(priority_genera, colnames(genus_clr))

corr_rows <- list(); idx <- 1L
for (g in priority_genera) {
  genus_delta <- genus_clr[a_ids, g] - genus_clr[b_ids, g]
  for (m in metric_names) {
    soil_delta <- soil_change[[paste0(m, "_change")]]
    test <- suppressWarnings(cor.test(genus_delta, soil_delta, method = "spearman", exact = FALSE))
    corr_rows[[idx]] <- data.frame(genus = g, genus_priority = screen$priority[match(g, screen$genus)],
                                   soil_metric = m, n_spatial_pairs = 9L,
                                   spearman_rho = unname(test$estimate),
                                   p_value_exploratory = test$p.value,
                                   stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
}
corr <- do.call(rbind, corr_rows)
corr$BH_q_exploratory <- p.adjust(corr$p_value_exploratory, method = "BH")
corr$greenhouse_n <- 3L
corr$limitation <- "spatial_pairs_nested_in_greenhouses_metric_names_D_to_I_pending"
corr <- corr[order(-abs(corr$spearman_rho)), ]
write.table(corr, file.path(out_tables, "genus_soil_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
signals <- corr[abs(corr$spearman_rho) >= 0.7 & corr$p_value_exploratory < 0.05, ]
write.table(signals, file.path(out_tables, "genus_soil_screening_signals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

alpha_rows <- list(); idx <- 1L
for (ss in unique(alpha_change$subset)) {
  z <- alpha_change[alpha_change$subset == ss, ]
  z <- z[match(soil_change$pair_id, z$pair_id), ]
  stopifnot(identical(z$pair_id, soil_change$pair_id))
  for (response in c("observed_change", "shannon_change", "simpson_change")) {
    for (m in metric_names) {
      test <- suppressWarnings(cor.test(z[[response]], soil_change[[paste0(m, "_change")]],
                                        method = "spearman", exact = FALSE))
      alpha_rows[[idx]] <- data.frame(subset = ss, diversity_response = response,
                                      soil_metric = m, n_spatial_pairs = 9L,
                                      spearman_rho = unname(test$estimate),
                                      p_value_exploratory = test$p.value,
                                      stringsAsFactors = FALSE)
      idx <- idx + 1L
    }
  }
}
alpha_corr <- do.call(rbind, alpha_rows)
alpha_corr$BH_q_exploratory <- p.adjust(alpha_corr$p_value_exploratory, method = "BH")
alpha_corr$limitation <- "exploratory_nested_spatial_pairs"
alpha_corr <- alpha_corr[order(-abs(alpha_corr$spearman_rho)), ]
write.table(alpha_corr, file.path(out_tables, "alpha_soil_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

max_by_genus <- aggregate(abs(spearman_rho) ~ genus, data = corr, FUN = max)
names(max_by_genus)[2] <- "max_abs_rho"
top_genera <- head(max_by_genus[order(-max_by_genus$max_abs_rho), "genus"], 20L)
heat <- corr[corr$genus %in% top_genera, ]
heat$genus <- factor(heat$genus, levels = rev(top_genera))
p_heat <- ggplot(heat, aes(soil_metric, genus, fill = spearman_rho)) +
  geom_tile(color = "white") + scale_fill_gradient2(low = "#2166AC", mid = "white",
                                                      high = "#B2182B", midpoint = 0,
                                                      limits = c(-1, 1)) +
  theme_bw(base_size = 9) + labs(x = "Soil metric change", y = NULL,
                                 fill = "Spearman rho") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(out_figures, "genus_soil_correlation_heatmap.pdf"),
       p_heat, width = 10, height = 8)

qc <- data.frame(metric = c("soil_samples", "paired_points", "soil_metrics",
                            "priority_genera_tested", "genus_metric_tests",
                            "unadjusted_screening_signals", "BH_q_below_0.05",
                            "yield_data_present", "disease_severity_numeric_present"),
                 value = c(18L, 9L, length(metric_names), length(priority_genera), nrow(corr),
                           nrow(signals), sum(corr$BH_q_exploratory < 0.05), 0L, 0L),
                 stringsAsFactors = FALSE)
write.table(qc, file.path(out_qc, "soil_integration_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(base_dir, "sessionInfo.txt"))
cat("ITS_SOIL_ASSOCIATION_V1_OK\n")

