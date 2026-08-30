#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({ library(vegan); library(permute) })
project <- "."
root <- file.path(project, "analysis", "downstream", "16S_v1_silva_diversity")
beta_out <- file.path(root, "beta", "MA_beta_rarefaction_sensitivity.tsv")
alpha_out <- file.path(root, "alpha", "MA_alpha_rarefaction_sensitivity_summary.tsv")
qc_out <- file.path(root, "qc", "validation_checks.tsv")
if (any(file.exists(c(beta_out, alpha_out, qc_out)))) stop("Refusing to overwrite validation outputs")

obj <- readRDS(file.path(root, "objects", "16S_analysis_object_v1.rds"))
counts <- obj$counts_clean
meta <- obj$metadata
results_beta <- list(); results_alpha <- list()

for (d in c(1500L, 1800L, 2500L)) {
  eligible <- rownames(counts)[rowSums(counts) >= d]
  ma <- meta[meta$cohort == "cultivation_paired" & meta$sample_id %in% eligible, ]
  complete_pairs <- names(which(table(ma$pair_id) == 2L))
  ma <- ma[ma$pair_id %in% complete_pairs, ]
  ids <- ma$sample_id
  set.seed(if (d == 1500L) 20260822L else 20260822L + d)
  rr <- rrarefy(counts[ids, , drop = FALSE], d)
  ma <- ma[match(rownames(rr), ma$sample_id), ]
  ma$condition <- factor(ma$condition, c("before_cultivation", "after_cultivation"))
  for (metric in c("Bray_Curtis", "Jaccard_binary")) {
    dd <- if (metric == "Bray_Curtis") vegdist(rr, "bray") else vegdist(rr, "jaccard", binary = TRUE)
    ctrl <- how(nperm = 9999, blocks = factor(ma$pair_id))
    fit <- adonis2(dd ~ condition, data = ma, permutations = ctrl)
    disp <- permutest(betadisper(dd, ma$condition), permutations = 9999)
    results_beta[[paste(d, metric)]] <- data.frame(
      depth = d, metric = metric, retained_MA_samples = nrow(ma), complete_pairs = length(complete_pairs),
      removed_pair_ids = paste(setdiff(unique(meta$pair_id[meta$cohort == "cultivation_paired"]), complete_pairs), collapse = ";"),
      R2 = fit$R2[1], pseudo_F = fit$F[1], permanova_p = fit$`Pr(>F)`[1],
      dispersion_F = disp$tab[1, "F"], dispersion_p = disp$tab[1, "Pr(>F)"],
      stringsAsFactors = FALSE)
  }
  av <- data.frame(sample_id = rownames(rr), observed = specnumber(rr), shannon = diversity(rr, "shannon"),
                   pair_id = ma$pair_id, condition = ma$condition, stringsAsFactors = FALSE)
  b <- av[av$condition == "before_cultivation", ]; a <- av[av$condition == "after_cultivation", ]
  b <- b[match(a$pair_id, b$pair_id), ]
  results_alpha[[as.character(d)]] <- data.frame(
    depth = d, complete_pairs = nrow(a),
    observed_mean_change = mean(a$observed - b$observed),
    observed_positive_pairs = sum(a$observed > b$observed),
    shannon_mean_change = mean(a$shannon - b$shannon),
    shannon_positive_pairs = sum(a$shannon > b$shannon), stringsAsFactors = FALSE)
}

write.table(do.call(rbind, results_beta), beta_out, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, results_alpha), alpha_out, sep = "\t", quote = FALSE, row.names = FALSE)

tax <- obj$taxonomy
checks <- data.frame(
  check = c("samples_24", "clean_asvs_9477", "all_depths_positive", "counts_integer",
            "taxonomy_order_matches", "rarefied_primary_depth_1500", "metadata_order_matches"),
  pass = c(nrow(counts) == 24L, ncol(counts) == 9477L, all(rowSums(counts) > 0),
           all(counts == round(counts)), identical(colnames(counts), tax$ASV_ID),
           all(rowSums(obj$counts_rarefied_depth1500) == 1500L),
           identical(rownames(counts), meta$sample_id)), stringsAsFactors = FALSE)
write.table(checks, qc_out, sep = "\t", quote = FALSE, row.names = FALSE)
if (!all(checks$pass)) stop("Validation check failed")
cat("MORCHELLA_16S_VALIDATION_V1_OK\n")

