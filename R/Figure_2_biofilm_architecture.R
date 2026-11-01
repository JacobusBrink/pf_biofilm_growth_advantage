# Introduction -----------------------------------------------------------

# Figure 2: biofilm architecture at 8 h with CAA (segmented images).
# The section at the end adds the mean +/- SD and contrast estimates quoted
# in the Results.
#
# Outputs: output/Figure_2_biofilm_architecture.png
#          output/Figure_2_biofilm_architecture_stats.xlsx

# This document offers a straightforward analysis of the biofilm architecture
# experiments performed with M9 minimal media (glucose, glycerol, and
# sodium succinate as carbon sources) using *P. fluorescens* wild-type.

# Each metric has 10
# colonies (technical replicates) imaged per well, all grown from the same
# overnight culture and media batch, nested within 3 independent days
# (biological replicates). Treating all 30 colonies per medium as
# independent observations ignores that nesting and overstates the sample
# size (pseudoreplication). A mixed model with `media` as a fixed effect and
# `biological_replicate` as a random intercept keeps all 30 colony-level
# observations (so statistical power isn't lost the way it would be by
# collapsing to n = 3 per medium and rerunning a rank test), while
# attributing the shared within-day variation to the random effect instead
# of letting it inflate the apparent evidence for a media difference.

# Load required packages -------------------------------------------------

library(tidyverse) # data manipulation and plotting
library(colorspace) # colour manipulation
library(patchwork) # combine multiple plots
library(ggpubr) # significance brackets via stat_pvalue_manual()
library(rstatix) # add_significance()
library(lme4) # mixed-effects models
library(emmeans) # estimated marginal means and pairwise contrasts
library(colorblindr) # palette_OkabeIto (colour-blind safe colours)

# Load data --------------------------------------------------------------

# The dataset has already been loaded and cleaned in the previous step of this pipeline.

# Let's load the clean version.
data <- readRDS(file = "input/clean_biofilm_data.rds")

# Theme and colour palette setup ------------------------------------------

# Load colour palettes, the custom theme_journal() theme, and geom defaults
source("R/00_setup.R")

# Subset data ------------------------------------------------------------

# Recode M9 media labels to short names and keep only the timepoint 8 (endpoint) measurements
biofilm_data <- data |>
  mutate(
    media = case_when(
      media == "Glucose (CAA 120 mM)" ~ "Glucose",
      media == "Glycerol (CAA 120 mM)" ~ "Glycerol",
      media == "Succinate (CAA 120 mM)" ~ "Succinate",
      TRUE ~ media
    )
  ) |>
  filter(
    timepoint == "8",
    media %in%
      c(
        "Glucose",
        "Glycerol",
        "Succinate"
      )
  ) |>
  mutate(
    media = as.factor(media),
    biological_replicate = as.factor(biological_replicate),
    biofilm_height = as.numeric(biofilm_height),
    biofilm_volume = as.numeric(biofilm_volume),
    biofilm_base_area = as.numeric(biofilm_base_area),
    distance_to_nearest_neighbor_mean = as.numeric(
      distance_to_nearest_neighbor_mean
    ),
    architecture_nematic_order_parameter_range30_mean = as.numeric(
      architecture_nematic_order_parameter_range30_mean
    ),
    distance_inter_cell_spacing_mean_range20_mean = as.numeric(
      distance_inter_cell_spacing_mean_range20_mean
    )
  )

# Statistical comparisons -------------------------------------------------

# All pairwise carbon-source comparisons, tested on every metric below via a
# linear mixed-effects model (family of 3 comparisons per metric,
# Bonferroni-adjusted)
carbon_source_comparisons <- list(
  c("Glucose", "Glycerol"),
  c("Glucose", "Succinate"),
  c("Glycerol", "Succinate")
)

# Fit media ~ (1 | biological_replicate) per metric, then get
# Bonferroni-adjusted pairwise contrasts between media from the fitted
# model's estimated marginal means. Output columns (group1, group2, p.adj,
# p.adj.signif) are shaped to match rstatix's wilcox_test() output so the
# same stat_pvalue_manual() plotting code and reporting table below still
# work unchanged.
run_mixed <- function(metric) {
  model <- lmer(
    as.formula(paste(metric, "~ media + (1 | biological_replicate)")),
    data = biofilm_data
  )

  emmeans(model, ~media) |>
    pairs(adjust = "bonferroni") |>
    as.data.frame() |>
    as_tibble() |>
    separate(contrast, into = c("group1", "group2"), sep = " - ") |>
    mutate(
      group1 = str_trim(group1),
      group2 = str_trim(group2),
      p.adj = p.value
    ) |>
    add_significance("p.adj") |>
    select(group1, group2, p.adj, p.adj.signif)
}

stats_height <- run_mixed("biofilm_height")
stats_volume <- run_mixed("biofilm_volume")
stats_base_area <- run_mixed("biofilm_base_area")
stats_nearest_neighbor <- run_mixed("distance_to_nearest_neighbor_mean")
stats_inter_cell_spacing <- run_mixed(
  "distance_inter_cell_spacing_mean_range20_mean"
)

# Biofilm architecture metrics by growth medium (timepoint 8) ------------

p1 <- ggplot(
  data = biofilm_data,
  aes(
    x = media,
    y = biofilm_volume,
    color = media
  )
) +

  geom_boxplot(
    aes(color = media),
    fill = NA,
    outlier.alpha = 0
  ) +

  # Plot individual data points with jitter
  geom_point(
    aes(fill = stage(media, after_scale = lighten(fill, .5))),
    position = position_jitterdodge(jitter.width = .2, seed = 42),
    shape = 21,
    size = 1,
    stroke = 0.5
  ) +

  # Axis labels (medium on x, metric on y)
  labs(
    x = NULL,
    y = expression(Volume ~ (mu * m^3))
  ) +

  # Constrain y-axis range and disable clipping
  coord_cartesian(
    ylim = c(0, 6500),
    clip = "off"
  ) +

  # Format y-axis labels with comma thousand separators
  scale_y_continuous(
    breaks = c(0, 2000, 4000, 6000),
    labels = scales::comma
  ) +

  # Fill points with a lightened shade of each medium's palette colour
  # scale_fill_manual(values = pal_bold) +
  scale_fill_OkabeIto() +

  # Colour boxplots, outlines and points using the Bold palette (one colour per medium)
  # scale_color_manual(values = pal_bold) +
  scale_color_OkabeIto() +

  # Pairwise mixed-model contrasts between carbon sources (Bonferroni-adjusted),
  # drawn as brackets with significance asterisks derived from p.adj
  stat_pvalue_manual(
    stats_volume,
    label = "p.adj.signif",
    y.position = c(5600, 6150, 6600),
    tip.length = 0.01,
    color = "grey30",
    size = 3,
    family = "Arial"
  ) +

  # Apply the custom minimal journal theme
  theme_journal(theme_name = "minimal")

p1

p2 <- ggplot(
  data = biofilm_data,
  aes(
    x = media,
    y = biofilm_base_area,
    color = media
  )
) +

  geom_boxplot(
    aes(color = media),
    fill = NA,
    outlier.alpha = 0
  ) +

  # Plot individual data points with jitter
  geom_point(
    aes(fill = stage(media, after_scale = lighten(fill, .5))),
    position = position_jitterdodge(jitter.width = .2, seed = 42),
    shape = 21,
    size = 1,
    stroke = 0.5
  ) +

  # Axis labels (medium on x, metric on y)
  labs(
    x = NULL,
    y = expression("Base Area" ~ (mu * m^2))
  ) +

  # Constrain y-axis range and disable clipping
  coord_cartesian(
    ylim = c(1500, 10500),
    clip = "off"
  ) +

  # Format y-axis labels with comma thousand separators
  scale_y_continuous(
    breaks = c(2500, 5000, 7500, 10000),
    labels = scales::comma
  ) +

  # Fill points with a lightened shade of each medium's palette colour
  # scale_fill_manual(values = pal_bold) +
  scale_fill_OkabeIto() +

  # Colour boxplots, outlines and points using the Bold palette (one colour per medium)
  # scale_color_manual(values = pal_bold) +
  scale_color_OkabeIto() +

  # Pairwise mixed-model contrasts between carbon sources (Bonferroni-adjusted),
  # drawn as brackets with significance asterisks derived from p.adj
  stat_pvalue_manual(
    stats_base_area,
    label = "p.adj.signif",
    y.position = c(8950, 9550, 10250),
    tip.length = 0.01,
    color = "grey30",
    size = 3,
    family = "Arial"
  ) +

  # Apply the custom minimal journal theme
  theme_journal(theme_name = "minimal")

p2

# Boxplot of biofilm height by medium, with jittered individual data points
p3 <- ggplot(
  data = biofilm_data,
  aes(
    x = media,
    y = biofilm_height,
    color = media
  )
) +

  geom_boxplot(
    aes(color = media),
    fill = NA,
    outlier.alpha = 0
  ) +

  # Plot individual data points with jitter
  geom_point(
    aes(fill = stage(media, after_scale = lighten(fill, .5))),
    position = position_jitterdodge(jitter.width = .2, seed = 42),
    shape = 21,
    size = 1,
    stroke = 0.5
  ) +

  # Axis labels (medium on x, metric on y)
  labs(
    x = NULL,
    y = expression("Height" ~ (mu * m))
  ) +

  # Add headroom above the data for the three stacked significance brackets
  # Constrain y-axis range and disable clipping
  coord_cartesian(
    ylim = c(5, 17),
    clip = "off"
  ) +

  # Set y-axis tick breaks
  scale_y_continuous(
    breaks = c(6, 9, 12, 15)
  ) +

  # Fill points with a lightened shade of each medium's palette colour
  # scale_fill_manual(values = pal_bold) +
  scale_fill_OkabeIto() +

  # Colour boxplots, outlines and points using the Bold palette (one colour per medium)
  # scale_color_manual(values = pal_bold) +
  scale_color_OkabeIto() +

  # Pairwise mixed-model contrasts between carbon sources (Bonferroni-adjusted),
  # drawn as brackets with significance asterisks derived from p.adj
  stat_pvalue_manual(
    stats_height,
    label = "p.adj.signif",
    y.position = c(14.5, 15.5, 16.5),
    tip.length = 0.01,
    color = "grey30",
    size = 3,
    family = "Arial"
  ) +

  # Apply the custom minimal journal theme
  theme_journal(theme_name = "minimal")

# Preview panel interactively (Positron's Plots pane handles the Arial font
# family fine)
p3

p4 <- ggplot(
  data = biofilm_data,
  aes(
    x = media,
    y = distance_to_nearest_neighbor_mean,
    color = media
  )
) +

  geom_boxplot(
    aes(color = media),
    fill = NA,
    outlier.alpha = 0
  ) +

  # Plot individual data points with jitter
  geom_point(
    aes(fill = stage(media, after_scale = lighten(fill, .5))),
    position = position_jitterdodge(jitter.width = .2, seed = 42),
    shape = 21,
    size = 1,
    stroke = 0.5
  ) +

  # Axis labels (medium on x, metric on y)
  labs(
    x = NULL,
    y = expression("Nearest Neighbour" ~ (mu * m))
  ) +

  # Constrain y-axis range and disable clipping
  coord_cartesian(
    ylim = c(1, 1.85),
    clip = "off"
  ) +

  # Set y-axis tick breaks
  scale_y_continuous(
    breaks = c(1, 1.25, 1.5, 1.75)
  ) +

  # Fill points with a lightened shade of each medium's palette colour
  # scale_fill_manual(values = pal_bold) +
  scale_fill_OkabeIto() +

  # Colour boxplots, outlines and points using the Bold palette (one colour per medium)
  # scale_color_manual(values = pal_bold) +
  scale_color_OkabeIto() +

  # Pairwise mixed-model contrasts between carbon sources (Bonferroni-adjusted),
  # drawn as brackets with significance asterisks derived from p.adj
  stat_pvalue_manual(
    stats_nearest_neighbor,
    label = "p.adj.signif",
    y.position = c(1.63, 1.7, 1.78),
    tip.length = 0.01,
    color = "grey30",
    size = 3,
    family = "Arial"
  ) +

  # Apply the custom minimal journal theme
  theme_journal(theme_name = "minimal")

p4

p5 <- ggplot(
  data = biofilm_data,
  aes(
    x = media,
    y = distance_inter_cell_spacing_mean_range20_mean,
    color = media
  )
) +

  geom_boxplot(
    aes(color = media),
    fill = NA,
    outlier.alpha = 0
  ) +

  # Plot individual data points with jitter
  geom_point(
    aes(fill = stage(media, after_scale = lighten(fill, .5))),
    position = position_jitterdodge(jitter.width = .2, seed = 42),
    shape = 21,
    size = 1,
    stroke = 0.5
  ) +

  # Axis labels (medium on x, metric on y)
  labs(
    x = NULL,
    y = expression("Inter-Cell Spacing" ~ (mu * m))
  ) +

  # Add headroom above the data for the three stacked significance brackets
  # Constrain y-axis range and disable clipping
  coord_cartesian(
    ylim = c(0.53, 0.635),
    clip = "off"
  ) +

  # Set y-axis tick breaks
  scale_y_continuous(
    breaks = c(0.54, 0.58, 0.62)
  ) +

  # Fill points with a lightened shade of each medium's palette colour
  # scale_fill_manual(values = pal_bold) +
  scale_fill_OkabeIto() +

  # Colour boxplots, outlines and points using the Bold palette (one colour per medium)
  # scale_color_manual(values = pal_bold) +
  scale_color_OkabeIto() +

  # Pairwise mixed-model contrasts between carbon sources (Bonferroni-adjusted),
  # drawn as brackets with significance asterisks derived from p.adj
  stat_pvalue_manual(
    stats_inter_cell_spacing,
    label = "p.adj.signif",
    y.position = c(0.615, 0.625, 0.634),
    tip.length = 0.01,
    color = "grey30",
    size = 3,
    family = "Arial"
  ) +

  # Apply the custom minimal journal theme
  theme_journal(theme_name = "minimal")

p5

# Combine plots with patchwork -------------------------------------------

# Arrange the six plots in a 2-row grid, collect a single shared legend,
# add lower-case letter tags, and apply a shared transparent theme via &
p7 <- wrap_plots(p1, p2, p3, p4, p5, nrow = 2) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = 'a', tag_suffix = ')') &
  theme(
    text = element_text(family = "Arial"),
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background = element_rect(fill = "transparent", colour = NA),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(family = "Arial"),
    # Add space below the a)/b)/c)... tags so they don't sit on top of long
    # rotated y-axis titles like "Inter-Cell Spacing (µm)" in panel e)
    plot.tag = element_text(margin = margin(r = 5), family = "Arial")
  )

p7 # draw combined plot

# Save the combined figure as a white-background PNG
ggsave(
  file.path(out_dir, "figures", "Figure_2_biofilm_architecture.png"),
  plot = p7,
  width = 180,
  height = 150,
  dpi = 300,
  units = "mm",
  bg = 'white',
  device = ragg::agg_png
)

# Statistical test results (for manuscript reporting) --------------------

# Start a character vector that will hold the six column names we want to
# test/summarise, one per panel a)-f); the c() function combines them into
# a single vector
architecture_metrics <- c(
  "biofilm_height",
  "biofilm_volume",
  "biofilm_base_area",
  "distance_to_nearest_neighbor_mean",
  "distance_inter_cell_spacing_mean_range20_mean"
)

# Stack the five mixed-model contrast tables computed earlier (stats_height,
# stats_volume, ...) into one table, in the same order as
# architecture_metrics, and add a "metric" column so each row records which
# column it tested.
stat_results <- bind_rows(
  stats_height |> mutate(metric = "biofilm_height", .before = 1),
  stats_volume |> mutate(metric = "biofilm_volume", .before = 1),
  stats_base_area |> mutate(metric = "biofilm_base_area", .before = 1),
  stats_nearest_neighbor |>
    mutate(metric = "distance_to_nearest_neighbor_mean", .before = 1),
  stats_inter_cell_spacing |>
    mutate(
      metric = "distance_inter_cell_spacing_mean_range20_mean",
      .before = 1
    )
)

# Build a second table of medians/IQRs/sample sizes, starting from the
# biofilm_data used in the plots and piping it through several steps
stat_summary_table <- biofilm_data |>

  # Keep only the medium label column and the five metric columns
  select(media, all_of(architecture_metrics)) |>

  # Reshape from "one column per metric" to two columns, "metric" (the
  # column name) and "value" (the number); every row becomes one measurement
  pivot_longer(-media, names_to = "metric", values_to = "value") |>

  # From here on, treat each metric/medium combination as its own group
  group_by(metric, media) |>

  # summarise() collapses each group down to one row per group
  summarise(
    # Count how many non-missing measurements are in this group
    n = sum(!is.na(value)),

    # Median (middle value) for this metric/medium group
    median = median(value, na.rm = TRUE),

    # Interquartile range: spread of the middle 50% of values in this group
    iqr = IQR(value, na.rm = TRUE),

    # Drop the grouping afterwards so the result is a normal flat table
    .groups = "drop"
  )

# Print the medians/IQR table to the R console; n = Inf means "show every
# row" instead of truncating the printout
print(stat_summary_table, n = Inf)

# Print the full pairwise test-results table to the console as well
print(stat_results, n = Inf)

# Numbers for the Results text -------------------------------------------

# Mean +/- SD per carbon source (n = 30 colonies, 10 x 3 biological
# replicates), and the estimated difference for each pair of carbon sources
# from the same mixed model as the brackets (Bonferroni-adjusted p).
mean_sd_table <- biofilm_data |>
  select(media, all_of(architecture_metrics)) |>
  pivot_longer(-media, names_to = "metric", values_to = "value") |>
  summarise(
    n = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .by = c(metric, media)
  )

contrast_table <- map(architecture_metrics, \(metric) {
  model <- lmer(
    as.formula(paste(metric, "~ media + (1 | biological_replicate)")),
    data = biofilm_data
  )
  emmeans(model, ~media) |>
    pairs(adjust = "bonferroni") |>
    as_tibble() |>
    transmute(
      metric = metric,
      contrast,
      estimate,
      SE,
      df,
      p_bonferroni = p.value
    )
}) |>
  list_rbind()

print(mean_sd_table, n = Inf)
print(contrast_table, n = Inf)

writexl::write_xlsx(
  list(
    mean_sd = mean_sd_table,
    contrasts = contrast_table,
    median_iqr = stat_summary_table,
    brackets = stat_results
  ),
  file.path(out_dir, "stats", "Figure_2_biofilm_architecture_stats.xlsx")
)
