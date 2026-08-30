options(stringsAsFactors=FALSE)
.libPaths(c("C:/path/to/workspace/Rlib", .libPaths()))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=2) stop("Usage: script <source_dir> <outdir>")
src <- args[1]; outdir <- args[2]; dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- read.delim(file.path(src,"MetaCyc_candidates_annotated_ranked.tsv"),check.names=FALSE)
core <- read.delim(file.path(src,"MetaCyc_core15_representative_pathways.tsv"),check.names=FALSE)
module_delta <- read.delim(file.path(src,"Figure8_source_module_paired_deltas.tsv"),check.names=FALSE)
nsti_pair <- read.delim(file.path(src,"Figure8_source_NSTI_paired_levels.tsv"),check.names=FALSE)
pal <- c("Higher after continuous cropping" = "#006BAD", "Higher before continuous cropping" = "#F8B9B8")
ghpal <- c(GH2 = "#5FA3CB", GH3 = "#668FCA", GH4 = "#C9CEFE")
theme_pub <- function() theme_classic(base_size = 6.3, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(face = "bold", size = 7.3),
        axis.title = element_text(size = 6.3), axis.text = element_text(size = 6.3),
        legend.title = element_text(size = 6.3), legend.text = element_text(size = 5.9),
        strip.text = element_text(face = "bold", size = 6.3), panel.grid = element_blank())

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
  labs(title = "Predicted pathway effects",
       subtitle = "215 robust pathways; 15 labelled representatives",
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
  labs(title = "Functional-module shifts", x = "Delta median CLR score (A - M)", y = NULL, colour = NULL) +
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
  labs(title = "NSTI sensitivity", subtitle = "9/9 pairs increased",
       y = "Weighted NSTI", colour = NULL) +
  theme_pub() + theme(legend.position = "none", axis.line.x = element_blank(), axis.ticks.x = element_blank(),
                      plot.margin = margin(5, 10, 8, 5))

fig <- free(p_a) /
  (p_b + p_c + plot_layout(widths = c(1.18, 0.82))) +
  plot_layout(heights = c(1.02, 1)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = "bold", size = 8))

base <- file.path(outdir, "Figure8_PICRUSt2_academic")
width_mm <- 169
height_mm <- 155
svglite::svglite(paste0(base, ".svg"), width = width_mm/25.4, height = height_mm/25.4)
print(fig); dev.off()
grDevices::cairo_pdf(paste0(base, ".pdf"), width = width_mm/25.4, height = height_mm/25.4, family = "Arial")
print(fig); dev.off()
ragg::agg_tiff(paste0(base, ".tiff"), width = width_mm/25.4, height = height_mm/25.4, units = "in", res = 600)
print(fig); dev.off()
ragg::agg_png(paste0(base, ".png"), width = width_mm/25.4, height = height_mm/25.4, units = "in", res = 300)
print(fig); dev.off()


