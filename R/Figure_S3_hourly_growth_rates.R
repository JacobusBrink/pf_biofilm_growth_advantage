# Introduction -----------------------------------------------------------

# Supplementary Figure S3: hour-to-hour specific growth rate, the basis for
# choosing the exponential-phase windows (Reviewer Comment 4).
#   hour-to-hour rate = ln(y_t / y_t-1) per hour, plotted at the midpoint of
#   each interval; y = OD600 x dilution factor (planktonic) or cells per
#   colony (biofilm), mean of the technical replicates.
#   - thin lines: biological replicates; thick lines: their mean
#   - grey shading: the window used to estimate mu
# Top block with CAA, bottom block without.
#
# Outputs: output/Figure_S3_hourly_growth_rates.png
#          output/Figure_S3_hourly_growth_rates.xlsx (the plotted rates)

source("R/00_setup.R")

# Hour-to-hour growth rates -----------------------------------------------------

growth <- replicate_means() |>
  left_join(windows, by = c("lifestyle", "carbon", "caa"))

hourly <- growth |>
  arrange(lifestyle, condition, bio_rep, time) |>
  mutate(
    mu_hourly = log(y / lag(y)) / (time - lag(time)),
    mid = time - 0.5,
    .by = c(lifestyle, condition, bio_rep)
  ) |>
  filter(!is.na(mu_hourly), is.finite(mu_hourly))

hourly_mean <- hourly |>
  summarise(
    mu_hourly = mean(mu_hourly),
    .by = c(lifestyle, carbon, caa, condition, mid, start, end)
  )

# Mean hour-to-hour rate per interval, printed as a table
options(width = 200)
print(
  hourly_mean |>
    arrange(mid) |>
    mutate(
      interval = paste0(mid - 0.5, "-", mid + 0.5, " h"),
      mu_hourly = round(mu_hourly, 2)
    ) |>
    select(lifestyle, condition, interval, mu_hourly) |>
    pivot_wider(names_from = interval, values_from = mu_hourly) |>
    arrange(desc(lifestyle), condition),
  n = Inf
)

write_xlsx(
  list(
    per_replicate = hourly |>
      transmute(lifestyle, carbon, caa, bio_rep, interval_mid = mid, mu_hourly),
    mean = hourly_mean |>
      transmute(
        lifestyle,
        carbon,
        caa,
        interval_mid = mid,
        mu_hourly,
        window_start = start,
        window_end = end
      )
  ),
  file.path(out_dir, "stats", "Figure_S3_hourly_growth_rates.xlsx")
)

# Figure ------------------------------------------------------------------------

fig_theme <- theme_journal(theme_name = "bw") +
  theme(
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA),
    plot.title = element_text(size = 8, face = "bold")
  )

hourly_plot <- function(cond_caa) {
  h <- hourly |> filter(caa == cond_caa)
  hm <- hourly_mean |> filter(caa == cond_caa)
  shade <- windows |>
    filter(caa == cond_caa) |>
    semi_join(h, by = c("lifestyle", "carbon"))

  ggplot(h, aes(mid, mu_hourly, colour = carbon)) +
    geom_rect(
      data = shade,
      aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
      inherit.aes = FALSE,
      fill = "grey88"
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    geom_line(aes(group = bio_rep), linewidth = 0.3, alpha = 0.6) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_line(data = hm, linewidth = 0.9) +
    facet_grid(
      factor(lifestyle, levels = c("Planktonic", "Biofilm")) ~ carbon,
      drop = FALSE
    ) +
    scale_x_continuous(breaks = 0:8, limits = c(0, 8)) +
    scale_colour_manual(values = carbon_colours) +
    labs(
      title = if (cond_caa == "+CAA") {
        "With Casamino Acids"
      } else {
        "Without Casamino Acids"
      },
      x = "Time (h)",
      y = expression("Hour-to-hour specific growth rate" ~ (h^-1))
    ) +
    fig_theme
}

fig_s3 <- hourly_plot("+CAA") / hourly_plot("-CAA")

ggsave(
  file.path(
    out_dir,
    "supplementary_figures",
    "Figure_S3_hourly_growth_rates.png"
  ),
  plot = fig_s3,
  width = 170,
  height = 200,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S3 to", out_dir, "\n")
