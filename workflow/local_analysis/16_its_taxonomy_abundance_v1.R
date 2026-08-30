#!/usr/bin/env Rscript

.libPaths(c("C:/path/to/R/library", .libPaths()))
suppressPackageStartupMessages({
  library(vegan)
  library(ape)
  library(ggplot2)
})

set.seed(20260821L)
project <- "."
count_file <- file.path(project, "analysis", "hpc_results", "ITS_dada2_v1_corrected_v1",
                        "results", "asv_count_table.tsv")
tax_file <- file.path(project, "analysis", "hpc_results", "ITS_UNITE_v2_job119642037",
                      "its_unite_v2", "its_asv_taxonomy.tsv")
metadata_file <- file.path(project, "analysis", "metadata", "sample_metadata.tsv")
result_dir <- file.path(project, "analysis", "downstream", "ITS_v2_taxonomy")
if (dir.exists(result_dir)) stop("Refusing to overwrite existing ITS_v2_taxonomy directory")
for (d in file.path(result_dir, c("qc", "tables", "diversity", "figures"))) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

count_tab <- read.delim(count_file, check.names = FALSE, stringsAsFactors = FALSE)
tax <- read.delim(tax_file, check.names = FALSE, stringsAsFactors = FALSE,
                  na.strings = character())
meta_all <- read.delim(metadata_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta_all[meta_all$marker == "ITS", ]
stopifnot(nrow(count_tab) == 8056L, nrow(tax) == 8056L, nrow(meta) == 24L)
stopifnot(!anyDuplicated(count_tab$ASV_ID), !anyDuplicated(tax$ASV_ID),
          setequal(count_tab$ASV_ID, tax$ASV_ID))

counts <- t(as.matrix(count_tab[, -1L, drop = FALSE]))
storage.mode(counts) <- "numeric"
colnames(counts) <- count_tab$ASV_ID
tax <- tax[match(colnames(counts), tax$ASV_ID), ]
stopifnot(identical(colnames(counts), tax$ASV_ID))
stopifnot(setequal(rownames(counts), meta$biological_sample_id))
meta <- meta[match(rownames(counts), meta$biological_sample_id), ]
meta$sample_id <- meta$biological_sample_id

has_kingdom <- nzchar(tax$Kingdom)
is_fungi <- tax$Kingdom == "Fungi"
has_genus <- nzchar(tax$Genus)
has_species <- nzchar(tax$Species)
high_species <- tax$confidence == "high_species_SH_candidate"
genus_candidate <- tax$confidence == "genus_candidate"

sample_depth <- rowSums(counts)
coverage <- data.frame(
  sample_id = rownames(counts),
  total_reads = sample_depth,
  kingdom_annotated_reads = rowSums(counts[, has_kingdom, drop = FALSE]),
  fungal_reads = rowSums(counts[, is_fungi, drop = FALSE]),
  genus_label_reads = rowSums(counts[, has_genus, drop = FALSE]),
  species_label_reads = rowSums(counts[, has_species, drop = FALSE]),
  high_species_candidate_reads = rowSums(counts[, high_species, drop = FALSE]),
  genus_candidate_reads = rowSums(counts[, genus_candidate, drop = FALSE]),
  stringsAsFactors = FALSE
)
for (name in names(coverage)[3:ncol(coverage)]) {
  coverage[[sub("_reads$", "_fraction", name)]] <- coverage[[name]] / coverage$total_reads
}
coverage <- cbind(coverage, meta[, c("cohort", "condition", "continuous_cropping_years",
                                    "greenhouse_id", "spatial_point", "pair_id",
                                    "independent_unit")])
write.table(coverage, file.path(result_dir, "qc", "sample_taxonomy_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

total_asv_reads <- colSums(counts)
overall <- data.frame(
  metric = c("total", "kingdom_annotated", "fungal", "genus_label", "species_label",
             "high_species_candidate", "genus_candidate", "unclassified"),
  asv_count = c(ncol(counts), sum(has_kingdom), sum(is_fungi), sum(has_genus),
                sum(has_species), sum(high_species), sum(genus_candidate), sum(!has_kingdom)),
  read_count = c(sum(counts), sum(total_asv_reads[has_kingdom]), sum(total_asv_reads[is_fungi]),
                 sum(total_asv_reads[has_genus]), sum(total_asv_reads[has_species]),
                 sum(total_asv_reads[high_species]), sum(total_asv_reads[genus_candidate]),
                 sum(total_asv_reads[!has_kingdom])),
  stringsAsFactors = FALSE
)
overall$asv_fraction <- overall$asv_count / ncol(counts)
overall$read_fraction <- overall$read_count / sum(counts)
write.table(overall, file.path(result_dir, "qc", "overall_taxonomy_coverage.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fungal_counts <- counts[, is_fungi, drop = FALSE]
fungal_tax <- tax[is_fungi, ]
fungal_count_out <- data.frame(ASV_ID = colnames(fungal_counts), t(fungal_counts),
                               check.names = FALSE)
write.table(fungal_count_out, file.path(result_dir, "tables", "fungal_asv_count_table.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(fungal_tax, file.path(result_dir, "tables", "fungal_asv_taxonomy.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

collapse_rank <- function(rank, top_n = 15L) {
  labels <- tax[[rank]]
  labels[!nzchar(labels)] <- "Unclassified_or_higher_rank"
  lev <- unique(labels)
  agg <- sapply(lev, function(x) rowSums(counts[, labels == x, drop = FALSE]))
  if (is.null(dim(agg))) agg <- matrix(agg, ncol = 1L, dimnames = list(rownames(counts), lev))
  totals <- colSums(agg)
  top <- names(sort(totals, decreasing = TRUE))[seq_len(min(top_n, length(totals)))]
  display <- ifelse(colnames(agg) %in% top, colnames(agg), "Other")
  collapsed <- sapply(unique(display), function(x) rowSums(agg[, display == x, drop = FALSE]))
  if (is.null(dim(collapsed))) collapsed <- matrix(collapsed, ncol = 1L,
                                                   dimnames = list(rownames(counts), unique(display)))
  rel <- collapsed / rowSums(collapsed)
  long <- data.frame(sample_id = rep(rownames(rel), times = ncol(rel)),
                     taxon = rep(colnames(rel), each = nrow(rel)),
                     relative_abundance = as.vector(rel), stringsAsFactors = FALSE)
  long <- merge(long, meta[, c("sample_id", "cohort", "condition", "greenhouse_id")],
                by = "sample_id", sort = FALSE)
  write.table(long, file.path(result_dir, "tables", paste0("composition_", tolower(rank), ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  long
}

phylum_long <- collapse_rank("Phylum", 12L)
genus_long <- collapse_rank("Genus", 15L)

alpha_metrics <- function(x, prefix) {
  data.frame(sample_id = rownames(x), subset = prefix, depth = rowSums(x),
             observed_asv = specnumber(x), shannon = diversity(x, "shannon"),
             simpson = diversity(x, "simpson"), stringsAsFactors = FALSE)
}
alpha <- rbind(alpha_metrics(counts, "all_ASVs"),
               alpha_metrics(fungal_counts, "UNITE_kingdom_Fungi"))
alpha <- merge(alpha, meta[, c("sample_id", "cohort", "condition", "greenhouse_id",
                               "pair_id", "independent_unit")], by = "sample_id", sort = FALSE)
write.table(alpha, file.path(result_dir, "diversity", "alpha_all_vs_fungal.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

fungal_relative <- fungal_counts / rowSums(fungal_counts)
fungal_bray <- vegdist(fungal_relative, method = "bray")
fit <- cmdscale(fungal_bray, k = 2, eig = TRUE, add = TRUE)
explained <- 100 * fit$eig[1:2] / sum(fit$eig[fit$eig > 0])
pcoa <- data.frame(sample_id = rownames(fit$points), axis1 = fit$points[, 1],
                   axis2 = fit$points[, 2], axis1_percent = explained[1],
                   axis2_percent = explained[2], stringsAsFactors = FALSE)
pcoa <- merge(pcoa, meta[, c("sample_id", "cohort", "condition", "greenhouse_id",
                             "pair_id", "independent_unit")], by = "sample_id", sort = FALSE)
write.table(pcoa, file.path(result_dir, "diversity", "fungal_bray_pcoa.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

plot_composition <- function(dat, rank) {
  ggplot(dat, aes(sample_id, relative_abundance, fill = taxon)) +
    geom_col(width = 0.9) + facet_grid(~cohort, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    theme_bw(base_size = 9) + labs(x = NULL, y = "Relative abundance", fill = rank) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1), panel.grid.major.x = element_blank())
}
ggsave(file.path(result_dir, "figures", "phylum_composition_all_reads.pdf"),
       plot_composition(phylum_long, "Phylum"), width = 12, height = 6)
ggsave(file.path(result_dir, "figures", "genus_composition_all_reads.pdf"),
       plot_composition(genus_long, "Genus"), width = 12, height = 7)

p_cov <- ggplot(coverage, aes(condition, fungal_fraction, color = greenhouse_id)) +
  geom_point(size = 2.5) + facet_wrap(~cohort, scales = "free_x") +
  scale_y_continuous(labels = function(x) paste0(round(100 * x), "%"), limits = c(0, 1)) +
  theme_bw(base_size = 10) + labs(x = NULL, y = "Reads assigned to Kingdom Fungi", color = "Greenhouse") +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))
ggsave(file.path(result_dir, "figures", "sample_fungal_annotation_coverage.pdf"),
       p_cov, width = 9, height = 5)

p_beta <- ggplot(pcoa, aes(axis1, axis2, color = condition, shape = greenhouse_id)) +
  geom_point(size = 3) + theme_bw(base_size = 11) +
  labs(x = sprintf("PCoA1 (%.1f%%)", explained[1]),
       y = sprintf("PCoA2 (%.1f%%)", explained[2]), color = "Condition", shape = "Greenhouse")
ggsave(file.path(result_dir, "figures", "fungal_bray_pcoa.pdf"), p_beta, width = 8, height = 6)

fungal_row <- overall[overall$metric == "fungal", ]
report <- c(
  "# ITS UNITE annotation and abundance QC", "",
  paste0("- Samples: ", nrow(counts), "; ASVs: ", ncol(counts), "; total reads: ", sum(counts), "."),
  paste0("- ASVs annotated at Kingdom: ", sum(has_kingdom), " (", sprintf("%.2f%%", 100 * mean(has_kingdom)), ")."),
  paste0("- Kingdom Fungi ASVs: ", sum(is_fungi), "; reads: ", fungal_row$read_count,
         " (", sprintf("%.2f%%", 100 * fungal_row$read_fraction), ")."),
  paste0("- Per-sample fungal read fraction range: ",
         sprintf("%.2f%%", 100 * min(coverage$fungal_fraction)), " to ",
         sprintf("%.2f%%", 100 * max(coverage$fungal_fraction)), "."),
  paste0("- ASVs with a Species label: ", sum(has_species),
         "; strict high-confidence species plus SH candidates: ", sum(high_species), "."),
  "", "## Interpretation constraints", "",
  "- Species names are candidate labels unless they pass the strict confidence rule.",
  "- Composition tables retain Unclassified_or_higher_rank reads.",
  "- Diversity is reported for all ASVs and the Kingdom Fungi subset as a sensitivity analysis.",
  "- M/A is paired, but spatial subsamples are nested within greenhouses; year and greenhouse are confounded.",
  "- H/D is an exploratory two-point comparison within one greenhouse.",
  "- No sample is deleted and no confirmatory differential-abundance test is run here."
)
writeLines(report, file.path(result_dir, "QC_AND_METHODS.md"), useBytes = TRUE)
writeLines(capture.output(sessionInfo()), file.path(result_dir, "sessionInfo.txt"))
cat("ITS_TAXONOMY_ABUNDANCE_V1_OK\n")

