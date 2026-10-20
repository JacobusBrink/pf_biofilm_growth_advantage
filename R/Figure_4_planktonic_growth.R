# Introduction -----------------------------------------------------------

# Figure 4: planktonic growth curves (OD600) with and without CAA
# (Reviewer Comment 12).
#   - points: every flask (3 per biological replicate x 3 replicates),
#     OD600 x dilution factor
#   - points outside the exponential window are faded (with CAA: glucose and
#     glycerol 2-5 h, succinate 2-6 h; without CAA: 3-7 h)
#   - line: drawn over the window only, from the mean intercept and mean
#     slope of the three replicate fits, so its slope is the reported mu
#   - label: mu (mean +/- SD of 3 replicates) and doubling time td
# Values are the same as in Table 2. No planktonic data were recorded for
# glycerol without CAA, so that panel is left empty.
#
# Output: output/Figure_4_planktonic_growth.png

source("R/00_setup.R")

# Data and fits ------------------------------------------------------------

flasks <- read_planktonic_flasks() |>
  mutate(lifestyle = "Planktonic") |>
  left_join(windows, by = c("lifestyle", "carbon", "caa")) |>
  mutate(in_window = time >= start & time <= end)

rates <- replicate_means() |>
  filter(lifestyle == "Planktonic") |>
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
  d <- flasks |> filter(carbon == cb, caa == cond)
  if (nrow(d) == 0) {
    return(plot_spacer()) # glycerol without CAA
  }
  l <- fit_lines |> filter(carbon == cb, caa == cond)
  s <- fit_summary |> filter(carbon == cb, caa == cond)
  col <- carbon_colours[[cb]]

  ggplot(d, aes(time, od)) +
    geom_point(
      aes(alpha = in_window),
      position = position_jitter(width = 0.1, height = 0, seed = 42),
      shape = 21,
      stroke = 0.4,
      size = 1.4,
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
      x = -0.3,
      y = Inf,
      label = s$label,
      hjust = 0,
      vjust = 1.2,
      size = 2.4,
      lineheight = 0.95,
      family = "Arial"
    ) +
    scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.3), guide = "none") +
    scale_y_log10(breaks = c(0.01, 0.1, 1)) +
    coord_cartesian(ylim = c(0.004, 4), clip = "off") +
    scale_x_continuous(breaks = 0:8, limits = c(-0.4, 8.4)) +
    labs(
      title = paste(cb, if (cond == "+CAA") "with CAA" else "without CAA"),
      x = if (show_x) "Time (h)" else NULL,
      y = if (show_y) expression("OD"[600]) else NULL
    ) +
    theme_journal(theme_name = "minimal") +
    theme(
      axis.text.x = element_text(colour = "black", size = 8),
      plot.title = element_text(size = 8)
    )
}

# Top row with CAA, bottom row without. Glycerol with CAA also gets an x
# label because the panel below it is empty.
panels <- expand_grid(
  caa = c("+CAA", "-CAA"),
  carbon = names(carbon_colours)
) |>
  pmap(\(caa, carbon) {
    growth_panel(
      carbon,
      caa,
      show_y = carbon == "Glucose",
      show_x = caa == "-CAA" | carbon == "Glycerol"
    )
  })

fig4 <- wrap_plots(panels, nrow = 2) +
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
  file.path(out_dir, "figures", "Figure_4_planktonic_growth.png"),
  plot = fig4,
  width = 180,
  height = 150,
  dpi = 300,
  units = "mm",
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure 4 to", out_dir, "\n")
