# Introduction -----------------------------------------------------------

# Supplementary Figure S8: OD600 and colony-forming units (CFU) in growing
# planktonic cultures with CAA (Reviewer Comment 6). Also the source of the
# Results numbers on CFU per OD600 and on mu from CFU vs mu from OD600.
#
# Experiment (Methods, OD600 Linearity and Viable Counts): 3 biological
# replicate cultures per carbon source, sampled at 0 h and hourly 2-6 h;
# OD600 and CFU (Miles and Misra drop plating, 3 x 20 uL drops) measured on
# the same sample.
#   CFU/mL = mean colonies per drop x 50 x 10^(dilution)
#   CFU per OD600: does it change over 2-6 h? ln(CFU/OD600) ~ culture + time,
#     one common slope; % change per hour = (exp(slope) - 1) x 100
#   mu from CFU and from OD600: log-linear slope per culture over the same
#     window as the main analysis (glucose and glycerol 2-5 h, succinate
#     2-6 h); paired difference mu_CFU - mu_OD, mean and 95% CI over the 9
#     cultures
#   (a) CFU and OD600 of each culture relative to its 2 h value
#   (b) mu_CFU against mu_OD per culture (error bars: SE of the slope;
#       dashed line: equality)
#
# Outputs: output/Figure_S8_od_vs_cfu.png
#          output/Figure_S8_od_vs_cfu.xlsx

source("R/00_setup.R")

cfu_file <- "input/od_vs_cfu/cfu_od_template.xlsx"

cfu_windows <- tibble(
  carbon = c("Glucose", "Glycerol", "Succinate"),
  start = c(2, 2, 2),
  end = c(5, 5, 6)
)

# Import ------------------------------------------------------------------

od <- read_excel(cfu_file, sheet = "od") |>
  transmute(
    carbon = str_to_title(carbon),
    bio_rep,
    time_h,
    od = od * od_dilution_factor
  )

cfu <- read_excel(cfu_file, sheet = "cfu") |>
  mutate(carbon = str_to_title(carbon)) |>
  pivot_longer(spot1:spot3, names_to = "spot", values_to = "colonies") |>
  filter(!is.na(colonies)) |>
  summarise(
    colonies_mean = mean(colonies),
    .by = c(carbon, bio_rep, time_h, dilution_exp)
  ) |>
  mutate(cfu_per_ml = colonies_mean * 50 * 10^abs(dilution_exp))

paired <- od |>
  left_join(cfu, by = c("carbon", "bio_rep", "time_h")) |>
  mutate(cfu_per_od = cfu_per_ml / od)

# CFU per OD600 over 2-6 h ----------------------------------------------------------

cfu_od_model <- lm(
  log(cfu_per_od) ~ interaction(carbon, bio_rep) + time_h,
  data = paired |> filter(time_h >= 2, time_h <= 6)
)
cfu_od_trend <- summary(cfu_od_model)$coefficients["time_h", ]
cfu_od_trend <- tibble(
  pct_per_hour = (exp(cfu_od_trend[["Estimate"]]) - 1) * 100,
  p = cfu_od_trend[["Pr(>|t|)"]]
)

# mu from CFU and from OD600 ------------------------------------------------------

fit_slope <- function(time, y) {
  s <- summary(lm(log(y) ~ time))$coefficients
  tibble(mu = s[2, 1], se = s[2, 2])
}

paired_mu <- paired |>
  inner_join(cfu_windows, by = "carbon") |>
  filter(time_h >= start, time_h <= end) |>
  reframe(
    cfu = fit_slope(time_h, cfu_per_ml),
    od = fit_slope(time_h, od),
    .by = c(carbon, bio_rep, start, end)
  ) |>
  unnest(c(cfu, od), names_sep = "_") |>
  mutate(diff = cfu_mu - od_mu)

paired_difference <- paired_mu |>
  summarise(
    mean = mean(diff),
    lo = mean(diff) - qt(0.975, n() - 1) * sd(diff) / sqrt(n()),
    hi = mean(diff) + qt(0.975, n() - 1) * sd(diff) / sqrt(n()),
    n = n()
  )

options(width = 200, pillar.sigfig = 3)
cat("\n== CFU per OD600 over 2-6 h (% change per hour) ==\n")
print(cfu_od_trend)
cat("\n== mu_CFU and mu_OD per culture ==\n")
print(paired_mu, n = Inf)
cat("\n== Paired difference mu_CFU - mu_OD (h-1), 95% CI ==\n")
print(paired_difference)

write_xlsx(
  list(
    cfu_per_od_trend = cfu_od_trend,
    paired_mu = paired_mu,
    paired_difference = paired_difference,
    timecourse = paired
  ),
  file.path(out_dir, "stats", "Figure_S8_od_vs_cfu.xlsx")
)

# Figure -----------------------------------------------------------------------

fig_theme <- theme_journal(theme_name = "bw") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    axis.text.x = element_text(colour = "black"),
    axis.text.y = element_text(colour = "black"),
    strip.background = element_rect(fill = "grey95", colour = NA)
  )

# (a) each measure relative to its own 2 h value (0 h not shown: not used
# in any analysis)
norm <- paired |>
  select(carbon, bio_rep, time_h, `OD600` = od, CFU = cfu_per_ml) |>
  pivot_longer(c(`OD600`, CFU), names_to = "measure", values_to = "y") |>
  mutate(y_norm = y / y[time_h == 2], .by = c(carbon, bio_rep, measure)) |>
  filter(time_h >= 2)

s8a <- ggplot(
  norm,
  aes(time_h, y_norm, colour = measure, group = interaction(measure, bio_rep))
) +
  geom_line(linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 1) +
  facet_wrap(~carbon) +
  scale_y_log10() +
  scale_x_continuous(breaks = 2:6) +
  scale_colour_manual(
    values = c("OD600" = "grey30", "CFU" = palette_OkabeIto[[5]]),
    breaks = c("CFU", "OD600"),
    labels = c("CFU", expression(OD[600]))
  ) +
  labs(x = "Time (h)", y = "Fold change relative to 2 h") +
  fig_theme +
  theme(
    legend.position = c(0.1, 0.9),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt")
  )

# (b) square plot with equal axes
lims <- range(c(
  paired_mu$od_mu - paired_mu$od_se,
  paired_mu$od_mu + paired_mu$od_se,
  paired_mu$cfu_mu - paired_mu$cfu_se,
  paired_mu$cfu_mu + paired_mu$cfu_se
))

s8b <- ggplot(paired_mu, aes(od_mu, cfu_mu, colour = carbon)) +
  geom_abline(linetype = "dashed", linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = cfu_mu - cfu_se, ymax = cfu_mu + cfu_se),
    width = 0,
    linewidth = 0.3
  ) +
  geom_errorbar(
    aes(xmin = od_mu - od_se, xmax = od_mu + od_se),
    width = 0,
    linewidth = 0.3,
    orientation = "y"
  ) +
  geom_point(size = 1.6) +
  coord_equal(xlim = lims, ylim = lims) +
  scale_colour_manual(values = carbon_colours) +
  labs(x = expression(mu[OD] ~ (h^-1)), y = expression(mu[CFU] ~ (h^-1))) +
  fig_theme +
  theme(
    legend.position = c(0.25, 0.85),
    legend.background = element_blank(),
    legend.key.spacing.y = unit(0, "pt"),
    legend.key.height = unit(10, "pt")
  )

fig_s8 <- (s8a | s8b) +
  plot_layout(widths = c(2.2, 1)) +
  plot_annotation(tag_levels = "a", tag_suffix = ')')

ggsave(
  file.path(out_dir, "supplementary_figures", "Figure_S8_od_vs_cfu.png"),
  plot = fig_s8,
  width = 180,
  height = 80,
  units = "mm",
  dpi = 300,
  bg = "white",
  device = ragg::agg_png
)

cat("\nSaved Figure S8 to", out_dir, "\n")
