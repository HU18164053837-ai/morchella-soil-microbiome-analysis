#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({ library(ggplot2); library(scales) })
set.seed(20260822L)
project <- "."
out16 <- file.path(project, "analysis", "downstream", "16S_v2_compositional_soil")
outjoint <- file.path(project, "analysis", "downstream", "joint_16S_ITS_soil_v2")
if (dir.exists(out16) || dir.exists(outjoint)) stop("Refusing to overwrite versioned outputs")
for (d in file.path(out16, c("tables", "figures", "qc", "objects"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
for (d in file.path(outjoint, c("tables", "figures", "qc"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)

obj_file <- file.path(project, "analysis", "downstream", "16S_v1_silva_diversity", "objects", "16S_analysis_object_v1.rds")
soil_file <- file.path(project, "analysis", "downstream", "ITS_v5_soil_integration_v2", "input", "soil_metrics_normalized_v2.tsv")
its_stage_file <- file.path(project, "analysis", "downstream", "ITS_v4_genus_screening", "tables", "MA_genus_screening.tsv")
its_soil_file <- file.path(project, "analysis", "downstream", "ITS_v5_soil_integration_v2", "tables", "genus_soil_spearman.tsv")
its_hd_file <- file.path(project, "analysis", "downstream", "ITS_v4_genus_screening", "tables", "HD_genus_screening.tsv")
required <- c(obj_file, soil_file, its_stage_file, its_soil_file, its_hd_file)
if (!all(file.exists(required))) stop("Missing input: ", paste(required[!file.exists(required)], collapse = "; "))

obj <- readRDS(obj_file)
counts <- obj$counts_clean
tax <- obj$taxonomy
meta <- obj$metadata
soil <- read.delim(soil_file, check.names = FALSE, stringsAsFactors = FALSE)
if (!all(soil$sample_id %in% meta$sample_id) || nrow(soil) != 18L) stop("Soil metadata mismatch")

genus <- tax$Genus
genus[is.na(genus) | !nzchar(genus)] <- "Unclassified_at_genus"
genus_counts <- t(rowsum(t(counts), genus, reorder = TRUE))
if (!all(rowSums(genus_counts) == rowSums(counts))) stop("Genus aggregation failed")
write.table(data.frame(sample_id = rownames(genus_counts), genus_counts, check.names = FALSE),
            file.path(out16, "objects", "16S_genus_count_table.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ma <- meta[meta$cohort == "cultivation_paired", ]
bmeta <- ma[ma$condition == "before_cultivation", ]
ameta <- ma[ma$condition == "after_cultivation", ]
bmeta <- bmeta[match(ameta$pair_id, bmeta$pair_id), ]
if (nrow(ameta) != 9L || !identical(bmeta$pair_id, ameta$pair_id)) stop("M/A pairing mismatch")
bids <- bmeta$sample_id; aids <- ameta$sample_id
maids <- c(bids, aids)
total_reads <- colSums(genus_counts)
prev_ma <- colSums(genus_counts[maids, , drop = FALSE] > 0)

clr_transform <- function(x, pc) {
  z <- log2(x + pc)
  z - rowMeans(z)
}
pcs <- c(0.1, 0.5, 1.0)
clr_list <- lapply(pcs, function(pc) clr_transform(genus_counts, pc)); names(clr_list) <- as.character(pcs)

signflip_p <- function(delta) {
  n <- length(delta)
  signs <- as.matrix(expand.grid(rep(list(c(-1, 1)), n)))
  null <- as.vector(signs %*% delta / n)
  mean(abs(null) >= abs(mean(delta)) - 1e-12)
}

stage_rows <- vector("list", ncol(genus_counts))
gh_rows <- list()
for (j in seq_len(ncol(genus_counts))) {
  deltas <- lapply(clr_list, function(z) z[aids, j] - z[bids, j])
  names(deltas) <- names(clr_list)
  d <- deltas[["0.5"]]
  gh <- aggregate(d, list(greenhouse_id = ameta$greenhouse_id), mean)
  names(gh)[2] <- "mean_delta_clr_pc05"
  gh$genus <- colnames(genus_counts)[j]
  gh_rows[[j]] <- gh[, c("genus", "greenhouse_id", "mean_delta_clr_pc05")]
  wp <- tryCatch(wilcox.test(d, mu = 0, paired = FALSE, exact = FALSE)$p.value, error = function(e) NA_real_)
  stage_rows[[j]] <- data.frame(
    genus = colnames(genus_counts)[j], total_reads_24 = total_reads[j], prevalence_MA = prev_ma[j],
    mean_relative_M = mean(genus_counts[bids, j] / rowSums(genus_counts)[bids]),
    mean_relative_A = mean(genus_counts[aids, j] / rowSums(genus_counts)[aids]),
    median_delta_clr_pc01 = median(deltas[["0.1"]]), median_delta_clr_pc05 = median(d),
    median_delta_clr_pc10 = median(deltas[["1"]]),
    pseudocount_direction_stable = length(unique(sign(c(median(deltas[["0.1"]]), median(d), median(deltas[["1"]]))))) == 1L,
    positive_pairs = sum(d > 0), negative_pairs = sum(d < 0),
    pair_direction_consistency = max(sum(d > 0), sum(d < 0)) / 9,
    greenhouse_direction_consistent = all(gh$mean_delta_clr_pc05 > 0) || all(gh$mean_delta_clr_pc05 < 0),
    paired_wilcoxon_p_exploratory = wp, exact_pair_signflip_p_exploratory = signflip_p(d),
    stringsAsFactors = FALSE)
}
stage <- do.call(rbind, stage_rows)
stage$paired_wilcoxon_BH <- p.adjust(stage$paired_wilcoxon_p_exploratory, "BH")
stage$exact_pair_signflip_BH <- p.adjust(stage$exact_pair_signflip_p_exploratory, "BH")
write.table(do.call(rbind, gh_rows), file.path(out16, "tables", "MA_genus_greenhouse_effects.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Dirichlet Monte-Carlo CLR: propagates count uncertainty; this is a custom paired implementation, not ALDEx2.
mc_n <- 128L
mc_effect <- matrix(NA_real_, nrow = mc_n, ncol = ncol(genus_counts), dimnames = list(NULL, colnames(genus_counts)))
for (k in seq_len(mc_n)) {
  prop <- matrix(rgamma(length(genus_counts), shape = as.vector(genus_counts) + 0.5, rate = 1),
                 nrow = nrow(genus_counts), ncol = ncol(genus_counts), dimnames = dimnames(genus_counts))
  prop <- prop / rowSums(prop)
  z <- log2(prop); z <- z - rowMeans(z)
  mc_effect[k, ] <- apply(z[aids, , drop = FALSE] - z[bids, , drop = FALSE], 2, median)
}
stage$mc_median_paired_effect <- apply(mc_effect, 2, median)[stage$genus]
stage$mc_effect_q025 <- apply(mc_effect, 2, quantile, 0.025)[stage$genus]
stage$mc_effect_q975 <- apply(mc_effect, 2, quantile, 0.975)[stage$genus]
stage$mc_direction_probability <- pmax(colMeans(mc_effect > 0), colMeans(mc_effect < 0))[stage$genus]
stage$mc_interval_excludes_zero <- (stage$mc_effect_q025 > 0 & stage$mc_effect_q975 > 0) |
  (stage$mc_effect_q025 < 0 & stage$mc_effect_q975 < 0)

eligible <- stage$genus != "Unclassified_at_genus" & stage$prevalence_MA >= 6L & stage$total_reads_24 >= 100L
stage$priority <- "not_prioritized"
moderate <- eligible & stage$pseudocount_direction_stable & stage$greenhouse_direction_consistent &
  stage$pair_direction_consistency >= 6/9 & abs(stage$median_delta_clr_pc05) >= 0.5 &
  stage$mc_direction_probability >= 0.90
high <- eligible & stage$pseudocount_direction_stable & stage$greenhouse_direction_consistent &
  stage$pair_direction_consistency >= 7/9 & abs(stage$median_delta_clr_pc05) >= 1 &
  stage$mc_direction_probability >= 0.95 & stage$mc_interval_excludes_zero
stage$priority[moderate] <- "moderate"
stage$priority[high] <- "high"
stage$direction <- ifelse(stage$median_delta_clr_pc05 > 0, "higher_relative_share_A", "lower_relative_share_A")
stage <- stage[order(factor(stage$priority, c("high", "moderate", "not_prioritized")),
                     -abs(stage$median_delta_clr_pc05)), ]
write.table(stage, file.path(out16, "tables", "MA_genus_compositional_screening.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(stage[stage$priority != "not_prioritized", ],
            file.path(out16, "tables", "MA_prioritized_genera_robust.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

prev_sens <- do.call(rbind, lapply(c(4L, 6L, 9L), function(th) {
  z <- stage[stage$genus != "Unclassified_at_genus" & stage$prevalence_MA >= th & stage$total_reads_24 >= 100L, ]
  data.frame(prevalence_threshold_of_18 = th, genera_testable = nrow(z),
             high_candidates = sum(z$priority == "high"), moderate_candidates = sum(z$priority == "moderate"),
             candidate_names = paste(z$genus[z$priority != "not_prioritized"], collapse = ";"), stringsAsFactors = FALSE)
}))
write.table(prev_sens, file.path(out16, "qc", "prevalence_threshold_sensitivity.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# H/D remains descriptive: two points in one greenhouse.
hids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "healthy"]
dids <- meta$sample_id[meta$cohort == "white_mold_case_example" & meta$condition == "diseased"]
z05 <- clr_list[["0.5"]]
hd <- data.frame(genus = colnames(genus_counts), total_reads_24 = total_reads,
                 prevalence_HD = colSums(genus_counts[c(hids, dids), , drop = FALSE] > 0),
                 mean_clr_H = colMeans(z05[hids, , drop = FALSE]), mean_clr_D = colMeans(z05[dids, , drop = FALSE]),
                 delta_clr_D_minus_H = colMeans(z05[dids, , drop = FALSE]) - colMeans(z05[hids, , drop = FALSE]),
                 stringsAsFactors = FALSE)
hd$descriptive_candidate <- hd$genus != "Unclassified_at_genus" & hd$prevalence_HD >= 2 &
  hd$total_reads_24 >= 100 & abs(hd$delta_clr_D_minus_H) >= 1
hd <- hd[order(-abs(hd$delta_clr_D_minus_H)), ]
write.table(hd, file.path(out16, "tables", "HD_genus_descriptive.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Soil association uses paired changes only, avoiding raw cross-phase correlation.
soil <- soil[match(c(bids, aids), soil$sample_id), ]
soil_metrics <- c("pH", "total_nitrogen", "total_phosphorus", "total_potassium",
                  "available_nitrogen", "available_phosphorus", "available_potassium", "organic_matter")
sb <- soil[match(bids, soil$sample_id), ]; sa <- soil[match(aids, soil$sample_id), ]
soil_delta <- as.data.frame(as.matrix(sa[, soil_metrics]) - as.matrix(sb[, soil_metrics]))
soil_delta$pair_id <- ameta$pair_id; soil_delta$greenhouse_id <- ameta$greenhouse_id
write.table(soil_delta, file.path(out16, "tables", "soil_paired_changes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

candidates <- stage[stage$priority != "not_prioritized", ]
soil_rows <- list(); idx <- 1L
for (g in candidates$genus) {
  delta_g <- z05[aids, g] - z05[bids, g]
  for (sm in soil_metrics) {
    ct <- suppressWarnings(cor.test(delta_g, soil_delta[[sm]], method = "spearman", exact = FALSE))
    soil_rows[[idx]] <- data.frame(genus = g, genus_priority = candidates$priority[candidates$genus == g],
                                   stage_direction = candidates$direction[candidates$genus == g],
                                   soil_metric = sm, n_spatial_pairs = 9L,
                                   spearman_rho = unname(ct$estimate), p_value_exploratory = ct$p.value,
                                   stringsAsFactors = FALSE)
    idx <- idx + 1L
  }
}
soil_assoc <- do.call(rbind, soil_rows)
soil_assoc$BH_q_exploratory <- p.adjust(soil_assoc$p_value_exploratory, "BH")
soil_assoc$screening_signal <- abs(soil_assoc$spearman_rho) >= 0.70 & soil_assoc$p_value_exploratory < 0.05
soil_assoc$limitation <- "9_spatial_pairs_nested_in_3_greenhouses;stage_batch_confounded;candidate_constrained_screen"
soil_assoc <- soil_assoc[order(soil_assoc$BH_q_exploratory, -abs(soil_assoc$spearman_rho)), ]
write.table(soil_assoc, file.path(out16, "tables", "MA_candidate_genus_soil_paired_spearman.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(soil_assoc[soil_assoc$screening_signal, ], file.path(out16, "tables", "MA_genus_soil_screening_signals.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

best_soil <- do.call(rbind, lapply(split(soil_assoc, soil_assoc$genus), function(z) z[which.max(abs(z$spearman_rho)), ]))
stage_best <- merge(candidates, best_soil[, c("genus", "soil_metric", "spearman_rho", "p_value_exploratory",
                                              "BH_q_exploratory", "screening_signal")], by = "genus", all.x = TRUE)
hdmatch <- hd[match(stage_best$genus, hd$genus), ]
stage_best$HD_delta_clr_D_minus_H <- hdmatch$delta_clr_D_minus_H
stage_best$HD_descriptive_candidate <- hdmatch$descriptive_candidate
stage_best$evidence_score <- ifelse(stage_best$priority == "high", 3L, 2L) +
  ifelse(stage_best$screening_signal, 1L, 0L) + ifelse(stage_best$BH_q_exploratory < 0.05, 2L, 0L) +
  ifelse(stage_best$HD_descriptive_candidate, 1L, 0L)
write.table(stage_best[order(-stage_best$evidence_score, -abs(stage_best$median_delta_clr_pc05)), ],
            file.path(out16, "tables", "16S_candidate_evidence_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Integrate with existing ITS evidence without altering frozen/current ITS outputs.
its_stage <- read.delim(its_stage_file, check.names = FALSE, stringsAsFactors = FALSE)
its_stage <- its_stage[its_stage$priority != "not_prioritized", ]
its_soil <- read.delim(its_soil_file, check.names = FALSE, stringsAsFactors = FALSE)
its_hd <- read.delim(its_hd_file, check.names = FALSE, stringsAsFactors = FALSE)
its_best <- do.call(rbind, lapply(split(its_soil[its_soil$genus %in% its_stage$genus, ],
                                       its_soil$genus[its_soil$genus %in% its_stage$genus]),
                                  function(z) z[which.max(abs(z$spearman_rho)), ]))

int16 <- data.frame(marker = "16S", taxon = stage_best$genus, stage_priority = stage_best$priority,
                    stage_effect = stage_best$median_delta_clr_pc05, stage_direction = stage_best$direction,
                    pair_direction_consistency = stage_best$pair_direction_consistency,
                    greenhouse_direction_consistent = stage_best$greenhouse_direction_consistent,
                    mc_direction_probability = stage_best$mc_direction_probability,
                    best_soil_metric = stage_best$soil_metric, soil_rho = stage_best$spearman_rho,
                    soil_p = stage_best$p_value_exploratory, soil_q = stage_best$BH_q_exploratory,
                    soil_screening_signal = stage_best$screening_signal,
                    HD_effect = stage_best$HD_delta_clr_D_minus_H,
                    HD_descriptive_candidate = stage_best$HD_descriptive_candidate,
                    evidence_score = stage_best$evidence_score, stringsAsFactors = FALSE)
im <- match(its_stage$genus, its_best$genus); ih <- match(its_stage$genus, its_hd$genus)
intits <- data.frame(marker = "ITS", taxon = its_stage$genus, stage_priority = its_stage$priority,
                    stage_effect = its_stage$median_paired_clr_change,
                    stage_direction = ifelse(its_stage$median_paired_clr_change > 0, "higher_relative_share_A", "lower_relative_share_A"),
                    pair_direction_consistency = its_stage$pair_direction_consistency,
                    greenhouse_direction_consistent = its_stage$greenhouse_direction_consistent,
                    mc_direction_probability = NA_real_, best_soil_metric = its_best$soil_metric[im],
                    soil_rho = its_best$spearman_rho[im], soil_p = its_best$p_value_exploratory[im],
                    soil_q = its_best$BH_q_exploratory[im],
                    soil_screening_signal = abs(its_best$spearman_rho[im]) >= 0.70 & its_best$p_value_exploratory[im] < 0.05,
                    HD_effect = its_hd$log2_mean_relative_ratio_D_vs_H[ih],
                    HD_descriptive_candidate = its_hd$priority[ih] != "not_prioritized",
                    evidence_score = ifelse(its_stage$priority == "high", 3L, 2L) +
                      ifelse(abs(its_best$spearman_rho[im]) >= 0.70 & its_best$p_value_exploratory[im] < 0.05, 1L, 0L) +
                      ifelse(its_best$BH_q_exploratory[im] < 0.05, 2L, 0L) +
                      ifelse(its_hd$priority[ih] != "not_prioritized", 1L, 0L), stringsAsFactors = FALSE)
integrated <- rbind(int16, intits)
integrated$soil_q_below_005 <- !is.na(integrated$soil_q) & integrated$soil_q < 0.05
integrated$interpretation <- "exploratory_relative_abundance_evidence_not_absolute_or_causal"
integrated <- integrated[order(-integrated$evidence_score, integrated$marker, -abs(integrated$stage_effect)), ]
write.table(integrated, file.path(outjoint, "tables", "integrated_16S_ITS_soil_candidate_evidence.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

metric_summary <- do.call(rbind, lapply(split(integrated, list(integrated$marker, integrated$best_soil_metric), drop = TRUE), function(z) {
  data.frame(marker = unique(z$marker), soil_metric = unique(z$best_soil_metric), candidate_count = nrow(z),
             positive_rho = sum(z$soil_rho > 0, na.rm = TRUE), negative_rho = sum(z$soil_rho < 0, na.rm = TRUE),
             screening_signals = sum(z$soil_screening_signal, na.rm = TRUE), stringsAsFactors = FALSE)
}))
write.table(metric_summary, file.path(outjoint, "tables", "candidate_best_soil_metric_summary.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Figures
theme_pub <- theme_bw(base_size = 10) + theme(panel.grid.minor = element_blank())
plot_stage <- candidates[order(candidates$median_delta_clr_pc05), ]
plot_stage$genus <- factor(plot_stage$genus, levels = plot_stage$genus)
p1 <- ggplot(plot_stage, aes(median_delta_clr_pc05, genus, color = priority)) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_errorbarh(aes(xmin = mc_effect_q025, xmax = mc_effect_q975), height = 0.2, alpha = 0.7) +
  geom_point(aes(size = pair_direction_consistency)) + theme_pub +
  labs(x = "Median paired CLR change (A - M)", y = NULL, color = "Priority", size = "Pair consistency")
ggsave(file.path(out16, "figures", "16S_MA_compositional_candidate_effects.pdf"), p1, width = 9, height = 9)

sig <- soil_assoc[soil_assoc$screening_signal, ]
if (nrow(sig) > 0) {
  sig$genus <- factor(sig$genus, levels = unique(sig$genus[order(sig$spearman_rho)]))
  p2 <- ggplot(sig, aes(soil_metric, genus, fill = spearman_rho)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 2.8) +
    scale_fill_gradient2(low = "#3B6FB6", mid = "white", high = "#C43C39", limits = c(-1, 1)) +
    theme_pub + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(x = NULL, y = NULL, fill = "Spearman rho")
  ggsave(file.path(out16, "figures", "16S_candidate_soil_screening_heatmap.pdf"), p2,
         width = 9, height = max(5, 0.32 * length(unique(sig$genus)) + 2))
}

topint <- do.call(rbind, lapply(split(integrated, integrated$marker), function(z) head(z, 18L)))
topint$label <- paste(topint$marker, topint$taxon, sep = ":")
topint$label <- factor(topint$label, levels = rev(topint$label))
p3 <- ggplot(topint, aes(stage_effect, label, color = marker, shape = soil_screening_signal)) +
  geom_vline(xintercept = 0, color = "grey70") + geom_point(aes(size = evidence_score), alpha = 0.85) +
  theme_pub + labs(x = "Stage-associated CLR effect (A - M)", y = NULL,
                   size = "Evidence score", shape = "Soil screen", color = "Marker")
ggsave(file.path(outjoint, "figures", "integrated_candidate_evidence.pdf"), p3, width = 10, height = 10)

qc <- data.frame(metric = c("16S_genera_total", "16S_testable_main", "16S_high", "16S_moderate",
                            "16S_soil_tests", "16S_soil_screening_signals", "16S_soil_BH_q_lt_005",
                            "ITS_candidates_integrated", "integrated_candidates", "yield_data_present",
                            "numeric_disease_outcome_present"),
                 value = c(ncol(genus_counts), sum(eligible), sum(stage$priority == "high"),
                           sum(stage$priority == "moderate"), nrow(soil_assoc), sum(soil_assoc$screening_signal),
                           sum(soil_assoc$BH_q_exploratory < 0.05), nrow(intits), nrow(integrated), 0, 0),
                 stringsAsFactors = FALSE)
write.table(qc, file.path(outjoint, "qc", "integration_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c("stage_effects_are_relative_not_absolute", "M_A_stage_fully_confounded_with_sequencing_batch",
             "cropping_year_fully_confounded_with_greenhouse", "three_spatial_subsamples_are_nested_not_independent",
             "HD_is_two_points_in_one_greenhouse", "soil_tests_use_9_paired_deltas",
             "no_yield_or_numeric_disease_outcome"), file.path(outjoint, "qc", "design_constraints.txt"))
writeLines(capture.output(sessionInfo()), file.path(out16, "sessionInfo.txt"))
writeLines(capture.output(sessionInfo()), file.path(outjoint, "sessionInfo.txt"))
cat("MORCHELLA_16S_COMPOSITIONAL_SOIL_JOINT_V1_OK\n")

