options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("Usage: Rscript 07b_plot_funguild_trophic_modes_matrix.R <data_dir> <out_prefix>")
data_dir <- args[1]
out_prefix <- args[2]

ab <- read.delim(file.path(data_dir, "strict_trophic_mode_sample_abundance.tsv"), check.names = FALSE)
st <- read.delim(file.path(data_dir, "FUNGuild_paired_function_statistics.tsv"), check.names = FALSE)
names(ab)[names(ab) == "function"] <- "trophic_mode"
names(st)[names(st) == "function"] <- "trophic_mode"
modes <- c("Pathotroph", "Saprotroph", "Symbiotroph")
mode_cols <- c("Pathotroph" = "#F8B9B8", "Saprotroph" = "#5FA3CB", "Symbiotroph" = "#C9CEFE")

theme_set(theme_classic(base_size = 6.8, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
        plot.title = element_text(face = "bold", size = 7.0, colour = "#273746")))

# a, mean-composition river: category widths encode average relative abundance.
a0 <- subset(ab, denominator == "annotated_reads")
means <- aggregate(relative_abundance ~ stage + trophic_mode, a0, mean)
means$stage <- factor(means$stage, levels = c("M", "A"))
means$trophic_mode <- factor(means$trophic_mode, levels = modes)
means <- means[order(means$stage, means$trophic_mode), ]
means$ymax <- ave(means$relative_abundance, means$stage, FUN = cumsum)
means$ymin <- means$ymax - means$relative_abundance

smoothstep <- function(t) 3*t^2 - 2*t^3
ribbons <- do.call(rbind, lapply(modes, function(md) {
  l <- means[means$stage == "M" & means$trophic_mode == md, ]
  r <- means[means$stage == "A" & means$trophic_mode == md, ]
  tt <- seq(0, 1, length.out = 80); ss <- smoothstep(tt)
  lo <- l$ymin + (r$ymin - l$ymin) * ss
  hi <- l$ymax + (r$ymax - l$ymax) * ss
  data.frame(trophic_mode = md, x = c(0.16 + 0.68*tt, rev(0.16 + 0.68*tt)),
             y = c(lo, rev(hi)))
}))

p_a <- ggplot() +
  geom_polygon(data = ribbons, aes(x, y, group = trophic_mode, fill = trophic_mode),
               alpha = 0.62, colour = "white", linewidth = 0.25) +
  geom_rect(data = means, aes(xmin = ifelse(stage == "M", 0.02, 0.86),
                              xmax = ifelse(stage == "M", 0.14, 0.98), ymin = ymin, ymax = ymax,
                              fill = trophic_mode), colour = "white", linewidth = 0.35) +
  geom_text(data = means, aes(x = ifelse(stage == "M", -0.01, 1.01), y = (ymin+ymax)/2,
                              label = sprintf("%.0f%%", 100*relative_abundance),
                              hjust = ifelse(stage == "M", 1, 0)), size = 2.2, colour = "#273746") +
  annotate("text", x = c(0.08, 0.92), y = 1.055, label = c("Before", "After"),
           fontface = "bold", size = 2.55, colour = "#273746") +
  scale_fill_manual(values = mode_cols, name = NULL) +
  coord_cartesian(xlim = c(-0.14, 1.14), ylim = c(0, 1.08), clip = "off") +
  labs(title = "Trophic-mode composition", x = NULL, y = NULL) +
  theme_void(base_family = "Arial", base_size = 6.8) +
  theme(plot.title = element_text(face = "bold", size = 7.0, colour = "#273746", hjust = 0, margin = margin(l = 10)),
        legend.position = "bottom", legend.text = element_text(size = 5.8),
        legend.key.size = grid::unit(3, "mm"), plot.margin = margin(8, 7, 5, 7))

# b, distribution silhouettes summarize all nine paired changes without a point-line panel.
wide_m <- subset(a0, stage == "M", select = c(pair_id, greenhouse_id, trophic_mode, relative_abundance))
wide_a <- subset(a0, stage == "A", select = c(pair_id, greenhouse_id, trophic_mode, relative_abundance))
names(wide_m)[4] <- "before"; names(wide_a)[4] <- "after"
delta <- merge(wide_m, wide_a, by = c("pair_id", "greenhouse_id", "trophic_mode"))
delta$change <- 100 * (delta$after - delta$before)
delta$pair_label <- paste0(delta$greenhouse_id, " | ", sub("^[234]-", "P", delta$pair_id))
delta$pair_label <- factor(delta$pair_label, levels = rev(unique(delta$pair_label[order(delta$greenhouse_id, delta$pair_id)])))
delta$trophic_mode <- factor(delta$trophic_mode, levels = modes)
direction_summary <- aggregate(change ~ trophic_mode, delta, function(z) paste0(sum(z > 0), "/9 increase"))

p_b <- ggplot(delta, aes(trophic_mode, change, fill = trophic_mode)) +
  geom_hline(yintercept = 0, colour = "#9AA5B1", linewidth = 0.4, linetype = 2) +
  geom_violin(width = 0.82, trim = FALSE, alpha = 0.72, colour = "white", linewidth = 0.45,
              adjust = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA, fill = "white", colour = "#273746",
               linewidth = 0.4, alpha = 0.9) +
  geom_text(data = direction_summary, aes(x = trophic_mode, y = 59, label = change),
            inherit.aes = FALSE, size = 2.1, fontface = "bold", colour = "#273746") +
  scale_fill_manual(values = mode_cols, guide = "none") +
  scale_y_continuous(breaks = c(-60, -30, 0, 30, 60), expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(ylim = c(-76, 64)) +
  labs(title = "Paired changes", x = NULL, y = "After - before (percentage points)") +
  theme_classic(base_size = 6.8, base_family = "Arial") +
  theme(panel.grid.major.y = element_line(colour = "#E7EBEF", linewidth = 0.3),
        panel.grid.minor = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 5.8, angle = 22, hjust = 1),
        axis.text.y = element_text(size = 5.8, colour = "#273746"), axis.title.y = element_text(size = 5.8),
        plot.title = element_text(face = "bold", size = 7.0, colour = "#273746"),
        plot.margin = margin(8, 5, 5, 5))

# c, a structured effect matrix compares confidence sets within each denominator.
c0 <- subset(st, ontology == "trophic_mode")
c0$trophic_mode <- factor(c0$trophic_mode, levels = modes)
c0$row_label <- paste(ifelse(c0$confidence_set == "strict", "Strict", "Inclusive"),
                      ifelse(c0$denominator == "annotated_reads", "annotated", "all reads"), sep = " | ")
c0$row_label <- factor(c0$row_label,
                       levels = rev(c("Strict | annotated", "Inclusive | annotated", "Strict | all reads", "Inclusive | all reads")))
c0$effect <- 100 * c0$mean_delta
c0$cell_label <- sprintf("%+.1f", c0$effect)
is_sap_ann <- c0$trophic_mode == "Saprotroph" & c0$denominator == "annotated_reads"
c0$cell_label[is_sap_ann] <- sprintf("%+.1f\nq=%.3f", c0$effect[is_sap_ann], c0$BH_q[is_sap_ann])
c0$confidence <- factor(ifelse(c0$confidence_set == "strict", "Strict", "Inclusive"),
                        levels = c("Inclusive", "Strict"))
c0$denominator_label <- factor(ifelse(c0$denominator == "annotated_reads", "Annotated reads", "All reads"),
                               levels = c("Annotated reads", "All reads"))
c0$display <- sprintf("%+.1f", c0$effect)
focus <- c0$trophic_mode == "Saprotroph" & c0$denominator == "annotated_reads"
c0$display[focus] <- sprintf("%+.1f\nq = %.3f", c0$effect[focus], c0$BH_q[focus])

p_c <- ggplot(c0, aes(trophic_mode, confidence, fill = effect)) +
  geom_tile(colour = "white", linewidth = 1.0, width = 0.96, height = 0.92) +
  geom_text(aes(label = display), size = 2.15, lineheight = 0.9, colour = "#273746") +
  facet_grid(denominator_label ~ ., switch = "y") +
  scale_fill_gradient2(low = "#F8B9B8", mid = "#F7F9FB", high = "#006BAD", midpoint = 0,
                       limits = c(-15, 15), oob = scales::squish,
                       name = "After - before\n(percentage points)") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5,
                               barwidth = grid::unit(45, "mm"), barheight = grid::unit(3, "mm"))) +
  labs(title = "Sensitivity analysis", x = NULL, y = NULL) +
  theme_minimal(base_size = 6.8, base_family = "Arial") +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(size = 5.8, angle = 25, hjust = 1, colour = "#273746"),
        axis.text.y = element_text(size = 5.8, colour = "#273746"),
        strip.placement = "outside", strip.background = element_rect(fill = "#EEF4F8", colour = NA),
        strip.text.y.left = element_text(angle = 90, face = "bold", size = 5.8, colour = "#273746"),
        legend.position = "bottom", legend.key.width = grid::unit(18, "mm"),
        legend.title = element_text(size = 5.8), legend.text = element_text(size = 5.6),
        plot.title = element_text(face = "bold", size = 7.0, colour = "#273746", margin = margin(l = 10)),
        plot.margin = margin(8, 5, 5, 5))

p_a <- p_a + theme(plot.title=element_blank(), plot.subtitle=element_blank())
p_b <- p_b + theme(plot.title=element_blank(), plot.subtitle=element_blank())
p_c <- p_c + theme(plot.title=element_blank(), plot.subtitle=element_blank())

make_page <- function() {
  grid::grid.newpage()
  lay <- grid::grid.layout(1, 3, widths = grid::unit(c(0.28, 0.31, 0.41), "null"))
  grid::pushViewport(grid::viewport(layout = lay))
  print(p_a, vp = grid::viewport(layout.pos.col = 1))
  print(p_b, vp = grid::viewport(layout.pos.col = 2))
  print(p_c, vp = grid::viewport(layout.pos.col = 3))
  grid::upViewport()
  for (z in list(c("a", 1.5), c("b", 47.5), c("c", 100.0)))
    grid::grid.text(z[1], x = grid::unit(as.numeric(z[2]), "mm"), y = grid::unit(117, "mm"),
                    just = c("left", "top"), gp = grid::gpar(fontfamily = "Arial", fontface = "bold", fontsize = 8))
}

width_in <- 169/25.4; height_in <- 120/25.4
dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
svglite::svglite(paste0(out_prefix, ".svg"), width = width_in, height = height_in); make_page(); dev.off()
grDevices::cairo_pdf(paste0(out_prefix, ".pdf"), width = width_in, height = height_in, family = "Arial"); make_page(); dev.off()
ragg::agg_tiff(paste0(out_prefix, ".tiff"), width = width_in, height = height_in, units = "in", res = 600, background = "white"); make_page(); dev.off()
ragg::agg_png(paste0(out_prefix, ".png"), width = width_in, height = height_in, units = "in", res = 300, background = "white"); make_page(); dev.off()



