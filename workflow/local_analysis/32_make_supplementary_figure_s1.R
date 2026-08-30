options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(svglite)
  library(ragg)
  library(grid)
})

args <- commandArgs(trailingOnly = TRUE)
project <- if (length(args)) args[1] else "."
down <- file.path(project, "analysis", "downstream")
out <- file.path(project, "analysis", "manuscript_v1", "figures_v2")
src <- file.path(out, "source_data")
dir.create(src, recursive = TRUE, showWarnings = FALSE)

soil <- read.delim(file.path(down, "ITS_v5_soil_integration_v2", "tables",
                            "soil_paired_changes.tsv"), check.names = FALSE)
evidence <- read.delim(file.path(down, "joint_16S_ITS_soil_v2", "tables",
                                "integrated_16S_ITS_soil_candidate_evidence.tsv"),
                       check.names = FALSE)
alpha16 <- read.delim(file.path(down, "16S_v1_silva_diversity", "alpha",
                               "alpha_diversity_depth1500.tsv"), check.names = FALSE)
alphaITS <- read.delim(file.path(down, "ITS_v2_taxonomy", "diversity",
                                "alpha_all_vs_fungal.tsv"), check.names = FALSE)

metric_key <- c(
  pH = "pH", total_nitrogen = "Total N", total_phosphorus = "Total P",
  total_potassium = "Total K", available_nitrogen = "Available N",
  available_phosphorus = "Available P", available_potassium = "Available K",
  organic_matter = "Organic matter"
)

change_cols <- paste0(names(metric_key), "_change")
soil_long <- do.call(rbind, lapply(seq_along(change_cols), function(i) {
  z <- soil[[change_cols[i]]]
  data.frame(pair_id = soil$pair_id, greenhouse_id = soil$greenhouse_id,
             metric = unname(metric_key[i]), raw_change = z,
             standardized_change = as.numeric(scale(z)))
}))
soil_long$metric <- factor(soil_long$metric, levels = unname(metric_key))

frozen <- c("Terrimonas", "Chitinophaga", "Gemmata", "Pseudarthrobacter",
            "Morchella", "Mortierella", "Alternaria", "Botryotrichum")
corr <- evidence[evidence$taxon %in% frozen,
                 c("marker", "taxon", "best_soil_metric", "soil_rho", "soil_q")]
corr$metric <- unname(metric_key[corr$best_soil_metric])
corr$metric <- factor(corr$metric, levels = unname(metric_key))
corr$taxon <- factor(corr$taxon, levels = rev(frozen))

# Expand the complete frozen-candidate x soil-metric matrix from marker-specific tables.
c16 <- read.delim(file.path(down, "16S_v2_compositional_soil", "tables",
                           "MA_candidate_genus_soil_paired_spearman.tsv"), check.names = FALSE)
cit <- read.delim(file.path(down, "ITS_v5_soil_integration_v2", "tables",
                           "genus_soil_spearman.tsv"), check.names = FALSE)
c16 <- c16[c16$genus %in% frozen, c("genus", "soil_metric", "spearman_rho", "BH_q_exploratory")]
cit <- cit[cit$genus %in% frozen, c("genus", "soil_metric", "spearman_rho", "BH_q_exploratory")]
corr_all <- rbind(c16, cit)
corr_all$metric <- factor(unname(metric_key[corr_all$soil_metric]), levels = unname(metric_key))
corr_all$genus <- factor(corr_all$genus, levels = rev(frozen))
stopifnot(all(corr_all$BH_q_exploratory > 0.05, na.rm = TRUE))

names(alpha16)[1] <- "sample_id"
a16 <- alpha16[alpha16$sample_id %in% c("H-1", "H-2", "H-3", "D-1", "D-2", "D-3"), ]
a16$condition <- ifelse(substr(a16$sample_id, 1, 1) == "H", "Healthy point", "Diseased point")
a16_long <- rbind(
  data.frame(sample_id = a16$sample_id, condition = a16$condition,
             metric = "Observed ASVs", value = a16$observed_asv),
  data.frame(sample_id = a16$sample_id, condition = a16$condition,
             metric = "Shannon index", value = a16$shannon)
)
a16_long$condition <- factor(a16_long$condition, levels = c("Healthy point", "Diseased point"))

names(alphaITS)[1] <- "sample_id"
subset_col <- if ("subset" %in% names(alphaITS)) "subset" else names(alphaITS)[2]
aits <- alphaITS[alphaITS$sample_id %in% c("H-1", "H-2", "H-3", "D-1", "D-2", "D-3") &
                   alphaITS[[subset_col]] == "UNITE_kingdom_Fungi", ]
aits$condition <- ifelse(substr(aits$sample_id, 1, 1) == "H", "Healthy point", "Diseased point")
obs_col <- intersect(c("observed_asv", "observed", "Observed"), names(aits))[1]
sha_col <- intersect(c("shannon", "Shannon"), names(aits))[1]
aITS_long <- rbind(
  data.frame(sample_id = aits$sample_id, condition = aits$condition,
             metric = "Observed ASVs", value = aits[[obs_col]]),
  data.frame(sample_id = aits$sample_id, condition = aits$condition,
             metric = "Shannon index", value = aits[[sha_col]])
)
aITS_long$condition <- factor(aITS_long$condition, levels = c("Healthy point", "Diseased point"))

write.table(soil_long, file.path(src, "SupplementaryFigureS1a_soil_paired_changes.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(corr_all, file.path(src, "SupplementaryFigureS1b_candidate_soil_correlations.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(a16_long, file.path(src, "SupplementaryFigureS1c_16S_HD_alpha.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
write.table(aITS_long, file.path(src, "SupplementaryFigureS1d_ITS_HD_alpha.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

pal_gh <- c(GH2 = "#4C78A8", GH3 = "#59A14F", GH4 = "#F28E2B")
pal_hd <- c("Healthy point" = "#4C78A8", "Diseased point" = "#D65F5F")

theme_set(theme_classic(base_size = 7, base_family = "Arial") +
            theme(axis.line = element_line(linewidth = 0.35),
                  axis.ticks = element_line(linewidth = 0.35),
                  strip.background = element_blank(),
                  strip.text = element_text(face = "bold", size = 6.5),
                  legend.title = element_text(size = 6.5),
                  legend.text = element_text(size = 6.2),
                  plot.title = element_text(face = "bold", size = 7.2),
                  plot.margin = margin(3, 4, 3, 4)))

p_a <- ggplot(soil_long, aes(standardized_change, metric, colour = greenhouse_id)) +
  geom_vline(xintercept = 0, colour = "grey70", linewidth = 0.35) +
  geom_point(size = 1.6, alpha = 0.9) +
  scale_colour_manual(values = pal_gh, name = "Greenhouse") +
  labs(x = "Standardized paired change (after - before)", y = NULL,
       title = "a  Paired soil-property changes") +
  theme(legend.position = "bottom")

p_b <- ggplot(corr_all, aes(metric, genus, fill = spearman_rho)) +
  geom_tile(colour = "white", linewidth = 0.35) +
  geom_text(aes(label = sprintf("%.2f", spearman_rho)), size = 1.7) +
  scale_fill_gradient2(low = "#4C78A8", mid = "white", high = "#D65F5F",
                       midpoint = 0, limits = c(-1, 1), name = "Spearman rho") +
  labs(x = NULL, y = NULL, title = "b  Frozen candidates and soil changes") +
  annotate("text", x = Inf, y = -Inf, label = "All BH q > 0.05",
           hjust = 1.05, vjust = -0.6, size = 2.1, colour = "grey25") +
  theme(axis.text.x = element_text(angle = 35, hjust = 1),
        legend.position = "bottom")

alpha_plot <- function(dat, ttl) {
  ggplot(dat, aes(condition, value, colour = condition)) +
    geom_point(position = position_jitter(width = 0.08, height = 0), size = 1.7) +
    stat_summary(fun = mean, geom = "crossbar", width = 0.42, linewidth = 0.45,
                 colour = "black") +
    facet_wrap(~metric, scales = "free_y", nrow = 1) +
    scale_colour_manual(values = pal_hd, guide = "none") +
    labs(x = NULL, y = NULL, title = ttl) +
    theme(axis.text.x = element_text(angle = 25, hjust = 1))
}

p_c <- alpha_plot(a16_long, "c  16S alpha diversity")
p_d <- alpha_plot(aITS_long, "d  ITS fungal alpha diversity")

draw_figure <- function() {
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(2, 2,
    heights = unit(c(0.60, 0.40), "npc"), widths = unit(c(0.48, 0.52), "npc"))))
  print(p_a, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_b, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  print(p_c, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
  print(p_d, vp = viewport(layout.pos.row = 2, layout.pos.col = 2))
  popViewport()
}

base <- file.path(out, "SupplementaryFigureS1_soil_and_white_mold_context")
w <- 183 / 25.4
h <- 142 / 25.4
svglite(paste0(base, ".svg"), width = w, height = h); draw_figure(); dev.off()
cairo_pdf(paste0(base, ".pdf"), width = w, height = h, family = "Arial"); draw_figure(); dev.off()
tiff(paste0(base, ".tiff"), width = w, height = h, units = "in", res = 600,
     compression = "lzw", type = "cairo"); draw_figure(); dev.off()
png(paste0(base, ".png"), width = w, height = h, units = "in", res = 300,
    type = "cairo"); draw_figure(); dev.off()

writeLines(c(
  "Supplementary Figure S1 QA notes",
  "Core conclusion: soil and white-mold point contrasts provide exploratory context only.",
  "Archetype: quantitative grid; R-only rendering; final size 183 x 142 mm.",
  "Panel a: 9 spatial pairs nested in 3 greenhouses; standardized only for display; raw changes retained in Source Data.",
  "Panel b: 8 frozen candidates x 8 soil metrics; Spearman n=9 spatial pairs; all BH q>0.05.",
  "Panels c-d: one healthy and one diseased spatial point in one greenhouse; three values per point are spatial subsamples, not independent disease replicates.",
  "No inferential disease-state test was performed. Points and means are shown; no error bars imply population uncertainty.",
  "No rows or observations were excluded from the planned panels."
), file.path(out, "SUPPLEMENTARY_FIGURE_S1_QA_REPORT_v1.md"))

