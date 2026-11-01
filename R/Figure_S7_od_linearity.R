# Introduction -----------------------------------------------------------

# Supplementary Figure S7: linearity of the OD600 readings (Libra S11)
# Also the source of the Results numbers on OD600
# linearity (log-log slope; deviation from proportionality at OD 0.3-0.4).
#
# Experiment (Methods, OD600 Linearity and Viable Counts): one culture per
# carbon source (with CAA) was serially diluted two-fold in fresh medium
# (8 concentrations, 3 independent dilution series, 3 readings each).
#   - log-log slope b: OD600 is proportional to concentration^b; b = 1 means
#     perfectly proportional. Fitted per series over the OD600 range used for
#     growth rates, then averaged (mean and 95% CI over the 9 series).
#   - deviation from proportionality: each reading x its dilution fold,
#     compared with the proportional line set by the three most dilute
#     points of the same series.
#   (a) OD600 against relative concentration (dashed: proportional line;
#       solid: calibration curve, a quadratic in log-log space)
#   (b) deviation of each reading from proportionality; grey shading marks
#       the OD600 range used for growth rate estimation
#
# Outputs: output/Figure_S7_od_linearity.png
#          output/Figure_S7_od_linearity.xlsx

source("R/00_setup.R")

lin_file <- "input/od_vs_cfu/od_linearity_template.xlsx"

# OD600 range of the readings used in the planktonic growth-rate fits
work_lo <- 0.024
work_hi <- 0.687

# Import ------------------------------------------------------------------

# Mean of the 3 readings of each dilution
series <- read_excel(lin_file, sheet = "dilution_series") |>
  mutate(carbon = str_to_title(carbon)) |>
  summarise(
    od = mean(od),
    .by = c(carbon, series, dilution_fold, relative_conc)
  ) |>
  mutate(back_calc = od * dilution_fold)

# Log-log slope ------------------------------------------------------------

mean_ci <- function(x) {
  n <- length(x)
  h <- qt(0.975, n - 1) * sd(x) / sqrt(n)
  tibble(
    b = mean(x),
    b_sd = sd(x),
    b_lo = mean(x) - h,
    b_hi = mean(x) + h,
    n = n
  )
}

loglog_slopes <- series |>
  filter(od >= work_lo, od <= work_hi) |>
  summarise(
    b = coef(lm(log(od) ~ log(relative_conc)))[2],
    .by = c(carbon, series)
  )

loglog_summary <- bind_rows(
  loglog_slopes |> reframe(mean_ci(b), .by = carbon),
  loglog_slopes |> reframe(mean_ci(b)) |> mutate(carbon = "All")
)

# Deviation from proportionality ------------------------------------------------

deviation <- series |>
  mutate(
    ref = mean(back_calc[dilution_fold >= 32]),
    deviation_pct = (back_calc / ref - 1) * 100,
    .by = c(carbon, series)
  )

deviation_summary <- deviation |>
  summarise(
    od = mean(od),
    deviation_pct = mean(deviation_pct),
    .by = c(carbon, dilution_fold)
  )

options(width = 200, pillar.sigfig = 3)
cat("\n== Log-log slope b (OD600", work_lo, "-", work_hi, ") ==\n")
print(loglog_summary)
cat("\n== Deviation from proportionality (mean of 3 series) ==\n")
print(
  deviation_summary |>
    pivot_wider(
      id_cols = dilution_fold,
      names_from = carbon,
      values_from = c(od, deviation_pct)
    ),
  n = Inf
)

write_xlsx(
  list(
    loglog_slope_summary = loglog_summary,
    loglog_slope_per_series = loglog_slopes,
    deviation = deviation
  ),
  file.path(out_dir, "stats", "Figure_S7_od_linearity.xlsx")
)

# Figure -----------------------------------------------------------------------

# Calibration curve per carbon source (quadratic in log-log space), drawn in (a)
calib_fits <- series |>
  nest(.by = carbon) |>
  mutate(
    fit = map(data, \(d) {
      lm(log(od) ~ log(relative_conc) + I(log(relative_conc)^2), data = d)
    })
  )

curve_pred <- calib_fits |>
  mutate(
    pred = map(fit, \(f) {
      tibble(relative_conc = exp(seq(log(1 / 128), 0, length.out = 100))) |>
        mutate(od = exp(predict(f, newdata = pick(everything()))))
    })
  ) |>
  select(carbon, pred) |>
  unnest(pred)

prop_line <- deviation |>
  summarise(ref = mean(ref), .by = carbon) |>
  cross_join(tibble(relative_conc = c(1 / 128, 1))) |>
  mutate(od = ref * relative_conc)

lin_theme <- theme_journal(theme_name = "bw") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

s7a <- ggplot(series, aes(relative_conc, od, colour = carbon)) +
  geom_line(data = prop_line, linetype = "dashed", linewidth = 0.35) +
  geom_line(data = curve_pred, linewidth = 0.4) +
  geom_point(size = 1.1, alpha = 0.8) +
  scale_x_log10(breaks = 2^-(0:7), labels = c("1", paste0("1/", 2^(1:7)))) +
  scale_y_log10() +
  scale_colour_manual(values = carbon_colours) +
  labs(
    x = "Relative concentration (2-fold serial dilution)",
    y = expression(OD[600] ~ "reading")
  ) +
  lin_theme

s7b <- ggplot(deviation, aes(od, deviation_pct, colour = carbon)) +
  annotate(
    "rect",
    xmin = work_lo,
    xmax = work_hi,
    ymin = -Inf,
    ymax = Inf,
    fill = "grey90",
    alpha = 0.6
  ) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_line(
    aes(group = interaction(carbon, series)),
    linewidth = 0.3,
    alpha = 0.6
  ) +
  geom_point(size = 1.1, alpha = 0.8) +
  scale_x_log10() +
  scale_colour_manual(values = carbon_colours) +
  labs(
    x = expression(OD[600] ~ "reading"),
    y = "Deviation from \nproportionality (%)"
  ) +
  lin_theme

fig_s7 <- ((s7a +
  guides(colour = guide_legend(override.aes = list(linetype = 0, alpha = 1)))) |
  (s7b + guides(colour = "none"))) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a", tag_suffix = ')') &
  theme(legend.position = "bottom")

ggsave(
  file.path(out_dir, "supplementary_figures", "Figure_S7_od_linearity.png"),
  plot = fig_s7,
  width = 180,
  height = 75,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S7 to", out_dir, "\n")
