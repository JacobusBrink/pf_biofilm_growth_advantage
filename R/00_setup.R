# Introduction -----------------------------------------------------------

# Shared setup for every script in R/. Each figure or table script
# starts with source("R/00_setup.R"), so the things that must be the
# same everywhere are defined once, here:
#   - packages used by most scripts
#   - the journal theme and the carbon-source colours
#   - the exponential-phase windows used to estimate the specific growth rate
#   - functions that read the raw data (biofilm cell counts, planktonic OD600)
#   - the log-linear fit and the biofilm vs planktonic statistical model
#
# Run every script from the repository root (the folder that contains R/ and
# input/). All outputs are written to output/.

# Load required packages -------------------------------------------------

suppressMessages({
  library(tidyverse) # data manipulation and plotting
  library(readxl) # read Excel input files
  library(writexl) # write Excel output files
  library(patchwork) # combine plots into one figure
  library(colorblindr) # palette_OkabeIto (colour-blind safe colours)
  library(colorspace) # lighten()
})

# Output folder ----------------------------------------------------------

out_dir <- "output"
subdirs <- c(
  "figures",
  "tables",
  "supplementary_figures",
  "supplementary_tables",
  "stats"
)

for (d in subdirs) {
  dir.create(file.path(out_dir, d), showWarnings = FALSE, recursive = TRUE)
}

# Theme and colours ------------------------------------------------------

theme_journal <- function(
  base_size = 9,
  base_family = "Arial",
  theme_name = "classic"
) {
  base_theme <- switch(
    theme_name,
    classic = theme_classic(base_size = base_size, base_family = base_family),
    minimal = theme_minimal(base_size = base_size, base_family = base_family),
    bw = theme_bw(base_size = base_size, base_family = base_family),
    stop("Unknown theme_name")
  )
  base_theme +
    theme(
      text = element_text(
        color = "black",
        size = base_size,
        family = base_family
      ),
      axis.text = element_text(color = "black", size = base_size * 0.9),
      axis.text.x = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(size = base_size, lineheight = 1.5),
      panel.background = element_rect(fill = "transparent", colour = NA),
      plot.background = element_rect(fill = "transparent", colour = NA),
      legend.position = "none"
    )
}

update_geom_defaults("point", list(size = 0.75))
update_geom_defaults("boxplot", list(linewidth = 0.25))

# One colour per carbon source, used in every figure
carbon_colours <- c(
  "Glucose" = palette_OkabeIto[[1]],
  "Glycerol" = palette_OkabeIto[[2]],
  "Succinate" = palette_OkabeIto[[3]]
)

# Exponential-phase windows ----------------------------------------------

# Timepoints (h) used to estimate the specific growth rate (Methods,
# Estimation of Specific Growth Rates; Supplementary Figures S1-S3).
#   Planktonic with CAA: from 2 h until the hour-to-hour rate dropped below
#     80% of its maximum (glucose and glycerol 2-5 h, succinate 2-6 h).
#   Planktonic without CAA: no distinct exponential phase, 3-7 h (same as
#     the biofilm window).
#   Biofilm: 3-7 h (only 2-5 cells per colony at 2 h; all conditions slowed
#     between 7 and 8 h).
# No planktonic data were recorded for glycerol without CAA.
windows <- tribble(
  ~lifestyle   , ~carbon     , ~caa   , ~start , ~end ,
  "Planktonic" , "Glucose"   , "+CAA" ,      2 ,    5 ,
  "Planktonic" , "Glycerol"  , "+CAA" ,      2 ,    5 ,
  "Planktonic" , "Succinate" , "+CAA" ,      2 ,    6 ,
  "Planktonic" , "Glucose"   , "-CAA" ,      3 ,    7 ,
  "Planktonic" , "Succinate" , "-CAA" ,      3 ,    7 ,
  "Biofilm"    , "Glucose"   , "+CAA" ,      3 ,    7 ,
  "Biofilm"    , "Glycerol"  , "+CAA" ,      3 ,    7 ,
  "Biofilm"    , "Succinate" , "+CAA" ,      3 ,    7 ,
  "Biofilm"    , "Glucose"   , "-CAA" ,      3 ,    7 ,
  "Biofilm"    , "Glycerol"  , "-CAA" ,      3 ,    7 ,
  "Biofilm"    , "Succinate" , "-CAA" ,      3 ,    7
)

# Order of the six conditions in figures and tables
condition_levels <- c(
  paste(names(carbon_colours), "+CAA"),
  paste(names(carbon_colours), "-CAA")
)

# Read the raw data --------------------------------------------------------

# Biofilm: one BiofilmQ summary CSV per colony (segmented images).
# File names are <carbon><caa|control>_<time>h_br<bio rep>_tr<tech rep>.csv,
# row 1 holds the units and row 2 the values. Returns one row per colony.
read_biofilm_colonies <- function(dir = "input/biofilm") {
  files <- list.files(dir, pattern = "\\.csv$", full.names = TRUE)
  map(files, \(f) {
    # suppressWarnings: BiofilmQ rows have uneven lengths, which read_csv
    # reports as parsing issues; the column used here is read correctly
    x <- suppressWarnings(read_csv(
      f,
      col_types = cols(.default = "c"),
      n_max = 2,
      col_select = Ncells
    ))
    tibble(file = basename(f), ncells = as.numeric(x$Ncells[2]))
  }) |>
    list_rbind() |>
    extract(
      file,
      c("media", "time", "bio_rep", "tech_rep"),
      "^([a-z]+)_([0-9]+)h_br([0-9])_tr([0-9]+)\\.csv$",
      convert = TRUE
    ) |>
    mutate(
      carbon = str_to_title(str_remove(media, "caa|control")),
      caa = if_else(str_detect(media, "control"), "-CAA", "+CAA")
    ) |>
    select(carbon, caa, bio_rep, tech_rep, time, ncells)
}

# Planktonic: OD600 of each flask (3 flasks per biological replicate).
# corrected_od = OD600 reading x dilution factor (samples were diluted when
# OD600 approached 0.8). Returns one row per flask and timepoint.
read_planktonic_flasks <- function(
  file = "input/planktonic/planktonic_glucose_glycerol_succinate.xlsx"
) {
  read_excel(file, na = "NA") |>
    mutate(
      caa = if_else(str_detect(carbon, "control"), "-CAA", "+CAA"),
      carbon = str_to_title(str_remove(carbon, "_control"))
    ) |>
    transmute(
      carbon,
      caa,
      bio_rep = biological_replicate,
      tech_rep = technical_replicate,
      time,
      od = corrected_od
    )
}

# Biological replicate means ------------------------------------------------

# Growth rates are fitted per biological replicate, to the mean of its
# technical replicates (10 colonies or 3 flasks) at each timepoint.
# y = cells per colony (biofilm) or OD600 x dilution factor (planktonic);
# y_sd = SD between the technical replicates.
replicate_means <- function() {
  biofilm <- read_biofilm_colonies() |>
    summarise(
      y = mean(ncells),
      y_sd = sd(ncells),
      .by = c(carbon, caa, bio_rep, time)
    ) |>
    mutate(lifestyle = "Biofilm", .before = 1)
  planktonic <- read_planktonic_flasks() |>
    summarise(
      y = mean(od),
      y_sd = sd(od),
      .by = c(carbon, caa, bio_rep, time)
    ) |>
    mutate(lifestyle = "Planktonic", .before = 1)
  bind_rows(biofilm, planktonic) |>
    mutate(condition = factor(paste(carbon, caa), levels = condition_levels))
}

# Specific growth rate --------------------------------------------------------

# mu = slope of ln(y) against time (h-1), fitted by linear regression to the
# timepoints between start and end, separately for each biological replicate.
# Returns one row per replicate with the intercept (for drawing the line),
# mu and R2.
fit_growth_rates <- function(d, win = windows) {
  d |>
    inner_join(win, by = c("lifestyle", "carbon", "caa")) |>
    filter(time >= start, time <= end) |>
    summarise(
      fit = list(lm(log(y) ~ time)),
      .by = c(lifestyle, carbon, caa, condition, bio_rep, start, end)
    ) |>
    mutate(
      intercept = map_dbl(fit, \(m) coef(m)[[1]]),
      mu = map_dbl(fit, \(m) coef(m)[[2]]),
      r2 = map_dbl(fit, \(m) summary(m)$r.squared)
    ) |>
    select(-fit)
}

# Biofilm vs planktonic comparison ---------------------------------------------

# Model used for every biofilm vs planktonic comparison (Methods, Statistical
# analysis): log(mu) of each replicate, one mean per condition x lifestyle
# group, generalised least squares with a separate variance for each group,
# Satterthwaite degrees of freedom. The contrast biofilm - planktonic on the
# log scale is back-transformed to a ratio, so (ratio - 1) x 100 is the %
# difference. 95% CIs and p-values are Bonferroni-adjusted over the
# comparisons in the family.
#
# Returns the fitted model and its estimated marginal means, so that other
# contrasts (Table S2, family D) can be taken from the same model.
fit_group_model <- function(rates) {
  d <- rates |>
    mutate(
      condition = factor(paste(carbon, caa)),
      group = interaction(condition, lifestyle, drop = TRUE, sep = "."),
      log_mu = log(mu)
    )
  model <- nlme::gls(
    log_mu ~ 0 + group,
    data = d,
    weights = nlme::varIdent(form = ~ 1 | group)
  )
  emm <- emmeans::emmeans(model, ~group, mode = "satterthwaite", data = d)
  list(data = d, model = model, emm = emm)
}

# Contrast weights: a named vector with +1/-1 on the groups being compared
contrast_weights <- function(terms, group_names) {
  w <- setNames(rep(0, length(group_names)), group_names)
  w[names(terms)] <- terms
  w
}

# Run a family of contrasts and Bonferroni-adjust within it
tidy_contrasts <- function(emm, weights, family) {
  emmeans::contrast(emm, method = weights) |>
    summary(infer = TRUE, adjust = "none") |>
    as_tibble() |>
    mutate(
      family = family,
      ratio = exp(estimate),
      t_crit = qt(1 - 0.05 / (2 * n()), df),
      lower = exp(estimate - t_crit * SE),
      upper = exp(estimate + t_crit * SE),
      p_bonferroni = pmin(p.value * n(), 1)
    )
}

# Biofilm vs planktonic for every condition that has both lifestyles
lifestyle_contrasts <- function(fitted) {
  gn <- levels(fitted$data$group)
  conds <- c(
    "Glucose +CAA",
    "Glycerol +CAA",
    "Succinate +CAA",
    "Glucose -CAA",
    "Succinate -CAA"
  )
  conds <- conds[paste0(conds, ".Planktonic") %in% gn]
  weights <- map(conds, \(cond) {
    contrast_weights(
      setNames(c(1, -1), paste0(cond, c(".Biofilm", ".Planktonic"))),
      gn
    )
  }) |>
    set_names(conds)
  tidy_contrasts(fitted$emm, weights, "Biofilm vs planktonic")
}

# Formatting helpers -------------------------------------------------------------

fmt_p <- function(p) {
  if_else(p < 0.001, "<0.001", if_else(p >= 0.995, "1.00", sprintf("%.3f", p)))
}
fmt_pct_ci <- function(r, lo, hi) {
  sprintf(
    "%+.0f%% (%+.0f to %+.0f%%)",
    (r - 1) * 100,
    (lo - 1) * 100,
    (hi - 1) * 100
  )
}
