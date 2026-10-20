# Introduction -----------------------------------------------------------

# Figure 3: biofilm growth curves (cells per colony) with and without CAA
# (Reviewer Comment 12).
#   - points: every colony (10 per biological replicate x 3 replicates)
#   - points outside the 3-7 h window are faded
#   - line: drawn over 3-7 h only, from the mean intercept and mean slope of
#     the three replicate fits, so its slope is the reported mu
#   - label: mu (mean +/- SD of 3 replicates) and doubling time td
# Values are the same as in Table 1.
#
# Output: output/Figure_3_biofilm_growth.png

source("R/00_setup.R")

# Data and fits ------------------------------------------------------------

colonies <- read_biofilm_colonies() |>
  mutate(lifestyle = "Biofilm") |>
  left_join(windows, by = c("lifestyle", "carbon", "caa")) |>
  mutate(in_window = time >= start & time <= end)

rates <- replicate_means() |>
  filter(lifestyle == "Biofilm") |>
  fit_growth_rates()

fit_summary <- rates |>
  summarise(
    intercept = mean(intercept),
    mu_sd = sd(mu),
    mu = mean(mu),
    .by = c(carbon, caa, start, end)
  ) |>
  mutate(
    td = log(2) / mu,
    label = sprintf(
      "µ = %.2f ± %s h⁻¹\ntd = %.2f h",
      mu,
      if_else(mu_sd < 0.005, "<0.01", sprintf("%.2f", mu_sd)),
      td
    )
  )

fit_lines <- fit_summary |>
  mutate(time = map2(start, end, \(a, z) seq(a, z, length.out = 30))) |>
  unnest(time) |>
  mutate(y = exp(intercept + mu * time))

# Panels -----------------------------------------------------------------------

growth_panel <- function(cb, cond, show_y, show_x) {
  d <- colonies |> filter(carbon == cb, caa == cond)
  l <- fit_lines |> filter(carbon == cb, caa == cond)
  s <- fit_summary |> filter(carbon == cb, caa == cond)
  col <- carbon_colours[[cb]]

  ggplot(d, aes(time, ncells)) +
    geom_point(
      aes(alpha = in_window),
      position = position_jitter(width = 0.1, height = 0, seed = 42),
      shape = 21,
      stroke = 0.4,
      size = 1.1,
      colour = col,
      fill = lighten(col, 0.5)
    ) +
    geom_line(
      data = l,
      aes(y = y),
      linewidth = 0.6,
      colour = "black",
      lineend = "round"
    ) +
    annotate(
      "text",
      x = 1.7,
      y = Inf,
      label = s$label,
      hjust = 0,
      vjust = 1.2,
      size = 2.4,
      lineheight = 0.95,
      family = "Arial"
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    scale_y_log10(breaks = c(1, 10, 100, 1000), labels = scales::comma) +
    coord_cartesian(ylim = c(1, 4000), clip = "off") +
    scale_x_continuous(breaks = 0:8, limits = c(1.6, 8.4)) +
    labs(
      title = paste(cb, if (cond == "+CAA") "with CAA" else "without CAA"),
      x = if (show_x) "Time (h)" else NULL,
      y = if (show_y) "Cells per colony" else NULL
    ) +
    theme_journal(theme_name = "minimal") +
    theme(
      axis.text.x = element_text(colour = "black", size = 8),
      plot.title = element_text(size = 8)
    )
}

# Top row with CAA, bottom row without; y label on the left column, x label
# on the bottom row
panels <- expand_grid(
  caa = c("+CAA", "-CAA"),
  carbon = names(carbon_colours)
) |>
  pmap(\(caa, carbon) {
    growth_panel(
      carbon,
      caa,
      show_y = carbon == "Glucose",
      show_x = caa == "-CAA"
    )
  })

fig3 <- wrap_plots(panels, nrow = 2) +
  plot_annotation(tag_levels = "a", tag_suffix = ")") &
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial"),
    plot.tag = element_text(margin = margin(r = 5), family = "Arial")
  )

ggsave(
  file.path(out_dir, "figures", "Figure_3_biofilm_growth.png"),
  plot = fig3,
  width = 180,
  height = 150,
  dpi = 300,
  units = "mm",
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure 3 to", out_dir, "\n")
