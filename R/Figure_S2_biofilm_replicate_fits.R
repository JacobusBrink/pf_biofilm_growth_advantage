# Introduction -----------------------------------------------------------

# Supplementary Figure S2: biofilm growth curve of every biological
# replicate (columns) and condition (rows), showing which timepoints were
# used to estimate the specific growth rate (Reviewer Comment 4).
#   - points: mean number of segmented cells per colony (10 colonies) +/- SD
#     (the SD is colony-to-colony spread, not measurement error)
#   - filled points (3-7 h): used for mu; open points: not used
#   - line: log-linear regression over 3-7 h; window, mu and R2 printed
#
# Output: output/Figure_S2_biofilm_replicate_fits.png

source("R/00_setup.R")

# Data and fits ------------------------------------------------------------

growth <- replicate_means() |>
  filter(lifestyle == "Biofilm") |>
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

fig_s2 <- ggplot(growth, aes(time, y, colour = carbon)) +
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
  scale_y_log10(labels = scales::label_comma()) +
  scale_x_continuous(breaks = 0:8) +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 21)) +
  scale_colour_manual(values = carbon_colours) +
  labs(x = "Time (h)", y = "Cells per colony (mean of 10)") +
  fig_theme

n_rows <- n_distinct(growth$condition)

ggsave(
  file.path(
    out_dir,
    "supplementary_figures",
    "Figure_S2_biofilm_replicate_fits.png"
  ),
  plot = fig_s2,
  width = 170,
  height = 20 + 42 * n_rows,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S2 to", out_dir, "\n")
