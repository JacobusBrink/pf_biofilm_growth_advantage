# Introduction -----------------------------------------------------------

# Supplementary Figure S1: planktonic growth curve of every biological
# replicate (columns) and condition (rows), showing which timepoints were
# used to estimate the specific growth rate.
#   - points: mean OD600 (x dilution factor) of 3 flasks +/- SD
#   - filled points: inside the exponential window, used for mu
#   - open points: outside the window
#   - line: log-linear regression over the window; window, mu and R2 printed
#
# Output: output/Figure_S1_planktonic_replicate_fits.png

source("R/00_setup.R")

# Data and fits ------------------------------------------------------------

growth <- replicate_means() |>
  filter(lifestyle == "Planktonic") |>
  left_join(windows, by = c("lifestyle", "carbon", "caa")) |>
  mutate(in_window = time >= start & time <= end)

fits <- fit_growth_rates(growth |> select(-start, -end))

fit_lines <- fits |>
  mutate(time = map2(start, end, \(a, z) seq(a, z, length.out = 20))) |>
  unnest(time) |>
  mutate(y = exp(intercept + mu * time))

labels <- fits |>
  mutate(
    label = sprintf("%d–%d h\nµ = %.2f h⁻¹\nR² = %.3f", start, end, mu, r2)
  )

# Figure --------------------------------------------------------------------

fig_theme <- theme_journal(theme_name = "bw") +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

fig_s1 <- ggplot(growth, aes(time, y, colour = carbon)) +
  # error bars floored at y/10 so that they stay visible on the log axis
  geom_errorbar(
    aes(ymin = pmax(y - y_sd, y / 10), ymax = y + y_sd),
    width = 0.2,
    linewidth = 0.3
  ) +
  geom_point(aes(shape = in_window), size = 1.4, fill = "white") +
  geom_line(data = fit_lines, linewidth = 0.5) +
  geom_text(
    data = labels,
    aes(x = min(growth$time) - 0.25, y = Inf, label = label),
    hjust = 0,
    vjust = 1.15,
    size = 2.1,
    lineheight = 0.9,
    colour = "black",
    family = "Arial"
  ) +
  facet_grid(condition ~ paste("Biological replicate", bio_rep)) +
  scale_y_log10() +
  scale_x_continuous(breaks = 0:8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
  scale_colour_manual(values = carbon_colours) +
  labs(x = "Time (h)", y = expression("OD"[600])) +
  fig_theme

n_rows <- n_distinct(growth$condition)

ggsave(
  file.path(
    out_dir,
    "supplementary_figures",
    "Figure_S1_planktonic_replicate_fits.png"
  ),
  plot = fig_s1,
  width = 170,
  height = 20 + 42 * n_rows,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S1 to", out_dir, "\n")
