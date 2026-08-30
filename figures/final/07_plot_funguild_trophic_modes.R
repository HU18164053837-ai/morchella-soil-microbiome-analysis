options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 07_plot_funguild_trophic_modes.R <data_dir> <out_prefix>")
data_dir <- args[1]
out_prefix <- args[2]

ab <- read.delim(file.path(data_dir, "strict_trophic_mode_sample_abundance.tsv"), check.names = FALSE)
st <- read.delim(file.path(data_dir, "FUNGuild_paired_function_statistics.tsv"), check.names = FALSE)
names(ab)[names(ab) == "function"] <- "trophic_mode"
names(st)[names(st) == "function"] <- "trophic_mode"

# Panel a uses all nine spatial pairs and the predeclared strict annotation set.
a <- subset(ab, denominator == "annotated_reads")
a$stage <- factor(a$stage, levels = c("M", "A"), labels = c("Before", "After"))
a$trophic_mode <- factor(a$trophic_mode, levels = c("Pathotroph", "Saprotroph", "Symbiotroph"))

ann <- subset(st, confidence_set == "strict" & ontology == "trophic_mode" & denominator == "annotated_reads")
ann$trophic_mode <- factor(ann$trophic_mode, levels = levels(a$trophic_mode))
ann$label <- sprintf("Delta %+.1f pp\nP = %.3f; q = %.3f",
                     100 * ann$mean_delta, ann$exact_p, ann$BH_q)

pal <- c("Before" = "#F8B9B8", "After" = "#5FA3CB")
p_a <- ggplot(a, aes(stage, 100 * relative_abundance, group = pair_id)) +
  geom_line(colour = "#B8C2CC", linewidth = 0.42, alpha = 0.72) +
  geom_point(aes(fill = stage), shape = 21, colour = "white", stroke = 0.35, size = 2.25) +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#273746", linewidth = 0.85) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", shape = 23, fill = "white",
               colour = "#273746", stroke = 0.6, size = 2.6) +
  geom_text(data = ann, aes(x = 1.5, y = Inf, label = label, group = NULL),
            inherit.aes = FALSE, vjust = 1.25, size = 2.15, colour = "#273746", lineheight = 0.92) +
  facet_wrap(~trophic_mode, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.23))) +
  labs(x = NULL, y = "Relative abundance among annotated reads (%)") +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        strip.background = element_rect(fill = "#EEF4F8", colour = NA),
        strip.text = element_text(face = "bold", size = 7.2, colour = "#273746"),
        panel.spacing.x = grid::unit(3.5, "mm"), plot.margin = margin(7, 4, 3, 5))

# Panel b shows sensitivity to confidence threshold, denominator, and greenhouse omission.
b <- subset(st, ontology == "trophic_mode")
b$confidence_set <- factor(b$confidence_set, levels = c("strict", "inclusive"), labels = c("Strict", "Inclusive"))
b$denominator <- factor(b$denominator, levels = c("annotated_reads", "all_reads"),
                        labels = c("Annotated reads", "All reads"))
b$trophic_mode <- factor(b$trophic_mode, levels = c("Symbiotroph", "Pathotroph", "Saprotroph"))
b$lo <- 100 * apply(b[, c("loo_GH2", "loo_GH3", "loo_GH4")], 1, min)
b$hi <- 100 * apply(b[, c("loo_GH2", "loo_GH3", "loo_GH4")], 1, max)
b$effect <- 100 * b$mean_delta
b$signal <- ifelse(b$effect >= 0, "Increase", "Decrease")

p_b <- ggplot(b, aes(effect, trophic_mode, colour = signal, shape = confidence_set)) +
  geom_vline(xintercept = 0, colour = "#9AA5B1", linewidth = 0.4, linetype = 2) +
  geom_segment(aes(x = lo, xend = hi, yend = trophic_mode), linewidth = 0.8, alpha = 0.75) +
  geom_point(fill = "white", size = 2.45, stroke = 0.75) +
  facet_grid(denominator ~ ., scales = "free_x") +
  scale_colour_manual(values = c("Increase" = "#006BAD", "Decrease" = "#F8B9B8")) +
  scale_shape_manual(values = c("Strict" = 16, "Inclusive" = 1)) +
  labs(x = "After - before effect (percentage points)", y = NULL,
       colour = "Direction", shape = "Annotation") +
  theme_classic(base_size = 7.2, base_family = "Arial") +
  theme(axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        strip.background = element_rect(fill = "#EEF4F8", colour = NA),
        strip.text = element_text(face = "bold", size = 6.8, colour = "#273746"),
        legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0,0,0,0),
        legend.key.width = grid::unit(3, "mm"), plot.margin = margin(7, 5, 3, 5))

make_page <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(nrow = 1, ncol = 2, widths = grid::unit(c(0.66, 0.34), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_a, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_b, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grid::upViewport()
  grid::grid.text("a", x = grid::unit(1.6, "mm"), y = grid::unit(109, "mm"),
                  just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9))
  grid::grid.text("b", x = grid::unit(122, "mm"), y = grid::unit(109, "mm"),
                  just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 9))
}

width_in <- 183 / 25.4
height_in <- 112 / 25.4
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
svglite::svglite(paste0(out_prefix, ".svg"), width = width_in, height = height_in)
make_page(); dev.off()
grDevices::cairo_pdf(paste0(out_prefix, ".pdf"), width = width_in, height = height_in, family = "Arial")
make_page(); dev.off()
ragg::agg_tiff(paste0(out_prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, background = "white")
make_page(); dev.off()
ragg::agg_png(paste0(out_prefix, ".png"), width = width_in, height = height_in, units = "in", res = 300, background = "white")
make_page(); dev.off()

