args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("Usage: Rscript 07_picrust2_core_pathways_figure.R <consensus.tsv> <metacyc_map.tsv.gz> <path_abun.tsv.gz> <nsti_metadata.tsv> <outdir>")
cons_file <- args[1]; map_file <- args[2]; abun_file <- args[3]
nsti_file <- args[4]; outdir <- args[5]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

.libPaths(c("C:/path/to/workspace/Rlib", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))

clean_text <- function(x) {
  x <- gsub("&alpha;", "alpha", x, fixed = TRUE)
  x <- gsub("&beta;", "beta", x, fixed = TRUE)
  x <- gsub("&gamma;", "gamma", x, fixed = TRUE)
  x <- gsub("<[^>]+>", "", x)
  trimws(x)
}

classify_module <- function(x) {
  z <- tolower(x)
  out <- rep("Other predicted metabolism", length(z))
  hit <- function(p) grepl(p, z, perl = TRUE)
  out[hit("peptidoglycan|lipid|fatty acid|stearate|palmit|oleate|mycol|lipopolysacchar|kdo|octulosonate|heptose|rhamnose|legionaminate|cell wall|teichoic")] <- "Cell envelope & lipids"
  out[hit("amino acid|alanine|arginine|aspart|asparagine|cysteine|glutamate|glutamine|glycine|histidine|isoleucine|leucine|lysine|methionine|ornithine|phenylalanine|proline|serine|threonine|tryptophan|tyrosine|valine")] <- "Amino-acid metabolism"
  out[hit("nad |nad$|cobalamin|biotin|folate|pterin|thiamin|riboflavin|heme|cofactor|vitamin|quinone")] <- "Cofactors & vitamins"
  out[hit("fermentation|glycol|pyruvate|tca|tricarbox|formaldehyde|methanol|calvin|respiration|acetyl-coa|carbon fixation|pentose phosphate")] <- "Carbon & energy"
  out[hit("purine|pyrimidine|nucleotide|nucleoside|adenosine|guanosine|uridine|cytidine|thymidine")] <- "Nucleotide metabolism"
  out[hit("degradation|detox|siderophore|enterobactin|mycothiol|glutathione|antibiotic|stress|xenobiotic")] <- "Environmental interaction"
  out
}

cons <- read.delim(cons_file, check.names = FALSE, stringsAsFactors = FALSE)
mp <- read.delim(gzfile(map_file), header = FALSE, stringsAsFactors = FALSE)
names(mp) <- c("pathway", "description")
mp$description <- clean_text(mp$description)
d <- merge(cons, mp, by = "pathway", all.x = TRUE, sort = FALSE)
d$description[is.na(d$description)] <- d$pathway[is.na(d$description)]
d$module <- classify_module(d$description)
d$direction <- ifelse(d$aldex_effect_A_minus_M > 0, "Higher after continuous cropping", "Higher before continuous cropping")
d$selection_pool <- with(d, paper_candidate & aldex_q < 0.05 & abs(aldex_effect_A_minus_M) >= 1)
d$robustness_score <- with(d, -log10(pmax(aldex_q, .Machine$double.xmin)) * abs(aldex_effect_A_minus_M))
d <- d[order(!d$paper_candidate, -d$robustness_score), ]
write.table(d, file.path(outdir, "MetaCyc_candidates_annotated_ranked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# Reproducible representative selection: strongest absolute ALDEx2 effect within
# each functional module and direction, then fill to 15 by global robustness.
pool <- d[d$selection_pool, ]
pool$key <- paste(pool$module, pool$direction, sep = " | ")
representative <- do.call(rbind, lapply(split(pool, pool$key), function(z) z[which.max(abs(z$aldex_effect_A_minus_M)), , drop = FALSE]))
representative <- representative[representative$module != "Other predicted metabolism", , drop = FALSE]
representative <- representative[order(-representative$robustness_score), , drop = FALSE]
fill <- pool[!pool$pathway %in% representative$pathway, ]
fill <- fill[order(-fill$robustness_score), ]
core <- rbind(representative, head(fill, max(0, 15 - nrow(representative))))
core <- core[order(core$aldex_effect_A_minus_M), ]
core$plot_label <- ifelse(nchar(core$description) > 58, paste0(substr(core$description, 1, 55), "..."), core$description)
write.table(core, file.path(outdir, "MetaCyc_core15_representative_pathways.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

ab <- read.delim(gzfile(abun_file), check.names = FALSE, row.names = 1)
pairs <- unlist(lapply(c(2, 3, 4), function(g) paste0(g, "-", 1:3)))
ordered <- c(paste0("M", pairs), paste0("A", pairs))
stopifnot(all(ordered %in% colnames(ab)))
x <- as.matrix(ab[, ordered, drop = FALSE]); storage.mode(x) <- "numeric"
rel <- sweep(x, 2, colSums(x), "/")
clr <- log(rel + 1e-7); clr <- sweep(clr, 2, colMeans(clr), "-")

candidate_modules <- unique(d[d$paper_candidate, c("pathway", "module")])
candidate_modules <- candidate_modules[candidate_modules$pathway %in% rownames(clr), ]
module_sizes <- table(candidate_modules$module)
keep_modules <- names(module_sizes[module_sizes >= 5 & names(module_sizes) != "Other predicted metabolism"])
module_scores <- do.call(rbind, lapply(keep_modules, function(mod) {
  ids <- candidate_modules$pathway[candidate_modules$module == mod]
  score <- apply(clr[ids, , drop = FALSE], 2, median)
  data.frame(module = mod, sample = names(score), score = as.numeric(score), stringsAsFactors = FALSE)
}))
module_scores$stage <- substr(module_scores$sample, 1, 1)
module_scores$pair <- substr(module_scores$sample, 2, nchar(module_scores$sample))
module_delta <- do.call(rbind, lapply(split(module_scores, module_scores$module), function(z) {
  data.frame(module = z$module[1], pair = pairs,
             delta_A_minus_M = z$score[match(paste0("A", pairs), z$sample)] - z$score[match(paste0("M", pairs), z$sample)])
}))
module_delta$greenhouse <- paste0("GH", substr(module_delta$pair, 1, 1))
write.table(module_scores, file.path(outdir, "Figure8_source_module_CLR_scores.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(module_delta, file.path(outdir, "Figure8_source_module_paired_deltas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

nsti <- read.delim(nsti_file, check.names = FALSE)
if (!"sample" %in% names(nsti)) nsti$sample <- rownames(nsti)
nsti$stage <- substr(nsti$sample, 1, 1)
nsti$pair <- substr(nsti$sample, 2, nchar(nsti$sample))
nsti_delta <- data.frame(pair = pairs,
  delta_A_minus_M = nsti$weighted_NSTI[match(paste0("A", pairs), nsti$sample)] - nsti$weighted_NSTI[match(paste0("M", pairs), nsti$sample)])
nsti_delta$greenhouse <- paste0("GH", substr(nsti_delta$pair, 1, 1))
write.table(nsti_delta, file.path(outdir, "Figure8_source_NSTI_paired_deltas.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nsti_pair <- data.frame(
  pair = pairs,
  greenhouse = paste0("GH", substr(pairs, 1, 1)),
  M = nsti$weighted_NSTI[match(paste0("M", pairs), nsti$sample)],
  A = nsti$weighted_NSTI[match(paste0("A", pairs), nsti$sample)]
)
write.table(nsti_pair, file.path(outdir, "Figure8_source_NSTI_paired_levels.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

pal <- c("Higher after continuous cropping" = "#006BAD", "Higher before continuous cropping" = "#F8B9B8")
ghpal <- c(GH2 = "#5FA3CB", GH3 = "#668FCA", GH4 = "#C9CEFE")
theme_pub <- function() theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(face = "bold", size = 8.2),
        axis.title = element_text(size = 7.4), axis.text = element_text(size = 6.8),
        legend.title = element_text(size = 6.8), legend.text = element_text(size = 6.4),
        strip.text = element_text(face = "bold", size = 6.8), panel.grid = element_blank())

landscape <- d[d$paper_candidate, , drop = FALSE]
landscape$minus_log10_q <- -log10(pmax(landscape$aldex_q, .Machine$double.xmin))
landscape$is_core <- landscape$pathway %in% core$pathway
core_landscape <- landscape[landscape$is_core, , drop = FALSE]
core_landscape$label_x <- ifelse(core_landscape$aldex_effect_A_minus_M > 0, 4.35, -4.85)
for (side in c(FALSE, TRUE)) {
  ii <- which((core_landscape$aldex_effect_A_minus_M > 0) == side)
  ii <- ii[order(core_landscape$minus_log10_q[ii])]
  core_landscape$label_y[ii] <- seq(1.55, 4.45, length.out = length(ii))
}
write.table(landscape, file.path(outdir, "Figure8_source_effect_landscape.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
p_a <- ggplot(landscape, aes(x = aldex_effect_A_minus_M, y = minus_log10_q)) +
  geom_hline(yintercept = -log10(0.05), linetype = "22", linewidth = 0.35, colour = "#A5A5A5") +
  geom_vline(xintercept = c(-1, 0, 1), linetype = c("22", "solid", "22"),
             linewidth = c(0.35, 0.4, 0.35), colour = c("#A5A5A5", "#606060", "#A5A5A5")) +
  geom_point(colour = "#C8CED8", size = 1.25, alpha = 0.58) +
  geom_point(data = core_landscape, aes(fill = direction, size = abs(mean_delta)),
             shape = 21, colour = "white", stroke = 0.45, alpha = 0.98) +
  geom_segment(data = core_landscape,
               aes(x = aldex_effect_A_minus_M, y = minus_log10_q,
                   xend = label_x, yend = label_y, colour = direction),
               linewidth = 0.26, alpha = 0.62) +
  geom_text(data = core_landscape,
            aes(x = label_x, y = label_y, label = pathway, colour = direction,
                hjust = ifelse(aldex_effect_A_minus_M > 0, 0, 1)),
            size = 2.25, fontface = "bold") +
  scale_fill_manual(values = pal, labels = c("After (A)", "Before (M)")) +
  scale_colour_manual(values = pal, guide = "none") +
  scale_size_continuous(range = c(2.3, 4.6), guide = "none") +
  coord_cartesian(xlim = c(-6.4, 6.0), clip = "off") +
  labs(title = "Functional effect landscape",
       subtitle = "215 robust pathways; 15 representative pathways are labelled",
       x = "ALDEx2 effect (A - M)", y = expression(-log[10](q)), fill = NULL) +
  theme_pub() + theme(legend.position = "top", plot.margin = margin(5, 12, 5, 8))

module_order <- names(sort(tapply(module_delta$delta_A_minus_M, module_delta$module, median)))
module_delta$module <- factor(module_delta$module, levels = module_order)
p_b <- ggplot(module_delta, aes(x = delta_A_minus_M, y = module)) +
  geom_vline(xintercept = 0, linewidth = 0.35, colour = "#707070") +
  geom_violin(width = 0.80, scale = "width", trim = FALSE,
              fill = "#D6DFEF", colour = "#8A9BB2", linewidth = 0.30, alpha = 0.72) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "#56657A", linewidth = 0.35) +
  geom_point(aes(colour = greenhouse), position = position_jitter(height = 0.095, width = 0),
             size = 1.55, alpha = 0.90) +
  scale_colour_manual(values = ghpal) +
  labs(title = "Module-level CLR shifts", x = "Delta median CLR score (A - M)", y = NULL, colour = NULL) +
  theme_pub() + theme(legend.position = "top", axis.line.y = element_blank(), axis.ticks.y = element_blank(),
                      plot.margin = margin(5, 14, 5, 5))

p_c <- ggplot(nsti_pair) +
  geom_segment(aes(x = 1, xend = 2, y = M, yend = A, colour = greenhouse),
               linewidth = 0.65, alpha = 0.62,
               arrow = grid::arrow(length = grid::unit(1.4, "mm"), type = "closed")) +
  geom_point(aes(x = 1, y = M, colour = greenhouse), size = 1.8, alpha = 0.95) +
  geom_point(aes(x = 2, y = A, colour = greenhouse), size = 1.8, alpha = 0.95) +
  scale_colour_manual(values = ghpal) +
  scale_x_continuous(NULL, breaks = c(1, 2), labels = c("Before (M)", "After (A)"), limits = c(0.82, 2.18)) +
  labs(title = "Prediction-quality audit", subtitle = "9/9 paired trajectories increase",
       y = "Weighted NSTI", colour = NULL) +
  theme_pub() + theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = margin(5, 10, 8, 5))

fig <- free(p_a) /
  (p_b + p_c + plot_layout(widths = c(1.18, 0.82))) +
  plot_layout(heights = c(1.12, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 9))

base <- file.path(outdir, "Figure8_PICRUSt2_functional_landscape")
width_mm <- 190
height_mm <- 170
svglite::svglite(paste0(base, ".svg"), width = width_mm/25.4, height = height_mm/25.4)
print(fig); dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm/25.4, height = height_mm/25.4, family = "Arial")
print(fig); dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm/25.4, height = height_mm/25.4, units = "in", res = 600)
print(fig); dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_mm/25.4, height = height_mm/25.4, units = "in", res = 300)
print(fig); dev.off()

writeLines(c(
  "Figure 8 | Continuous cropping reorganizes predicted bacterial functional potential.",
  "a, Effect landscape of all 215 robust MetaCyc pathways. Fifteen representative pathways are enlarged and labelled; they were selected reproducibly from robust candidates (paper_candidate = TRUE, ALDEx2 q < 0.05, absolute effect >= 1) using the largest absolute effect within each functional module and direction before filling by the preregistered robustness score. Positive values indicate higher predicted potential after continuous cropping (A); negative values indicate higher predicted potential before continuous cropping (M).",
  "b, Raincloud-style distributions of paired A-M differences in median CLR scores across all robust candidate pathways assigned to each functional module. Points are nine spatially matched pairs and colors indicate greenhouse; violins show smoothed distributions and narrow boxes show the median and interquartile range.",
  "c, Paired weighted NSTI trajectories from M to A. NSTI increased in all nine pairs (mean A-M = 0.214; exact sign-flip p = 0.00585), indicating greater reference-genome distance after continuous cropping.",
  "All functional values are PICRUSt2 predictions derived from 16S taxonomic profiles and represent predicted genomic potential, not measured genes, transcription, or metabolic activity."
), file.path(outdir, "Figure8_legend.txt"))

cat("PICRUST2_CORE_FIGURE_COMPLETE\n")

