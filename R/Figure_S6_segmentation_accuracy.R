# Introduction -----------------------------------------------------------

# Supplementary Figure S6: segmentation accuracy per image, across carbon
# sources and time points.
#   (a) automated vs manual cell count per image (line: 1:1)
#   (b) count error and (c) recall per image over time; the grey band marks
#       the 3-7 h window used for the biofilm growth rate
#   Filled circles: with CAA; open circles: without CAA.
# R version of the panel made by validation/final/figures_final.py (same
# data and layout).
#
# Input: validation/final/results/partA_per_image.csv (see
#        Table_S3_segmentation_accuracy.R)
# Output: output/Figure_S6_segmentation_accuracy.png

source("R/00_setup.R")

part_a <- read_csv(
  "validation/final/results/partA_per_image.csv",
  show_col_types = FALSE
) |>
  mutate(
    carbon = str_to_title(carbon),
    caa = if_else(medium == "caa", "+CAA", "-CAA"),
    # small horizontal offset so carbon sources at the same time do not overlap
    time_jit = timepoint_h +
      c(Glucose = -0.15, Glycerol = 0, Succinate = 0.15)[carbon]
  )

fig_theme <- theme_journal(theme_name = "classic") +
  theme(
    axis.text.x = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 9)
  )

point_layers <- list(
  geom_point(
    aes(colour = carbon, fill = interaction(carbon, caa)),
    shape = 21,
    size = 1.8,
    stroke = 0.5
  ),
  scale_colour_manual(values = carbon_colours),
  # filled with CAA, white without
  scale_fill_manual(
    values = c(
      setNames(carbon_colours, paste0(names(carbon_colours), ".+CAA")),
      setNames(rep("white", 3), paste0(names(carbon_colours), ".-CAA"))
    ),
    guide = "none"
  )
)

lim <- c(1, max(part_a$manual_cells, part_a$automated_objects) * 1.3)

s6a <- ggplot(part_a, aes(manual_cells, automated_objects)) +
  geom_abline(colour = "grey50", linewidth = 0.4) +
  point_layers +
  scale_x_log10(limits = lim) +
  scale_y_log10(limits = lim) +
  coord_equal() +
  labs(x = "Manual count", y = "Automated count", colour = NULL) +
  guides(
    colour = guide_legend(override.aes = list(fill = unname(carbon_colours)))
  ) +
  fig_theme +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt")
  )

over_time <- function(var, ylab) {
  ggplot(part_a, aes(time_jit, .data[[var]])) +
    annotate(
      "rect",
      xmin = 2.7,
      xmax = 7.3,
      ymin = -Inf,
      ymax = Inf,
      fill = "grey93"
    ) +
    geom_line(
      data = part_a |> filter(caa == "+CAA"),
      aes(colour = carbon),
      linewidth = 0.4,
      alpha = 0.7
    ) +
    point_layers +
    scale_x_continuous(breaks = 2:8) +
    labs(x = "Time (h)", y = ylab) +
    fig_theme
}

s6b <- over_time("count_error_pct", "Count error (%)") +
  geom_hline(yintercept = 0, colour = "grey50", linewidth = 0.4)
s6c <- over_time("recall", "Recall") + coord_cartesian(ylim = c(0, 1.05))

fig_s6 <- (s6a |
  s6b |
  s6c) +
  plot_annotation(tag_levels = "a", tag_suffix = ')')

ggsave(
  file.path(
    out_dir,
    "supplementary_figures",
    "Figure_S6_segmentation_accuracy.png"
  ),
  plot = fig_s6,
  width = 180,
  height = 64,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S6 to", out_dir, "\n")
