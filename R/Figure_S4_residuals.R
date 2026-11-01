# Introduction -----------------------------------------------------------

# Supplementary Figure S4: residuals of the log-linear regressions used to
# estimate the specific growth rate.
#   Upper rows: the selected windows (00_setup.R).
#   Lower rows: the original 2-8 h window of the first submission.
# Residual = observed ln(y) minus the fitted line, per biological replicate.
# Points scatter around 0 when growth is exponential over the window; an
# arch (positive in the middle, negative at the ends) means growth slowed
# within the window.
#
# Output: output/Figure_S4_residuals.png

source("R/00_setup.R")

# Residuals for both windows --------------------------------------------------

growth <- replicate_means()

window_sets <- bind_rows(
  windows |> mutate(scenario = "Selected window"),
  windows |> mutate(start = 2, end = 8, scenario = "Original 2-8 h window")
)

residuals_long <- window_sets |>
  inner_join(
    growth,
    by = c("lifestyle", "carbon", "caa"),
    relationship = "many-to-many"
  ) |>
  filter(time >= start, time <= end) |>
  mutate(
    residual = resid(lm(log(y) ~ time)),
    .by = c(scenario, lifestyle, carbon, caa, bio_rep)
  ) |>
  mutate(
    scenario = factor(
      scenario,
      levels = c("Selected window", "Original 2-8 h window")
    )
  )

# Figure ------------------------------------------------------------------------

fig_theme <- theme_journal(theme_name = "bw") +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.title = element_text(size = 8, face = "bold")
  )

residual_plot <- function(life) {
  d <- residuals_long |> filter(lifestyle == life)
  m <- d |>
    summarise(
      residual = mean(residual),
      .by = c(scenario, condition, carbon, time)
    )
  ggplot(d, aes(time, residual, colour = carbon)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_point(aes(shape = factor(bio_rep)), size = 1.2, alpha = 0.7) +
    geom_line(data = m, linewidth = 0.6) +
    facet_grid(scenario ~ condition) +
    scale_x_continuous(breaks = 2:8, limits = c(1.7, 8.3)) +
    scale_colour_manual(values = carbon_colours) +
    labs(title = life, x = "Time (h)", y = "Residual of ln(y)") +
    fig_theme
}

fig_s4 <- residual_plot("Biofilm") / residual_plot("Planktonic")

ggsave(
  file.path(out_dir, "supplementary_figures", "Figure_S4_residuals.png"),
  plot = fig_s4,
  width = 180,
  height = 170,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S4 to", out_dir, "\n")
