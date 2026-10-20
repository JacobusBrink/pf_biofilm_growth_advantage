# Runs every script in R/ and rebuilds output/ from the raw data.
# Run from the repository root: source("R/run_all.R")
# (or from a terminal: Rscript R/run_all.R)

scripts <- c(
  "Figure_2_biofilm_architecture.R",
  "Figure_3_biofilm_growth.R",
  "Figure_4_planktonic_growth.R",
  "Table_1_biofilm_growth_rates.R",
  "Table_2_planktonic_growth_rates.R",
  "Figure_S1_planktonic_replicate_fits.R",
  "Figure_S2_biofilm_replicate_fits.R",
  "Figure_S3_hourly_growth_rates.R",
  "Figure_S4_residuals.R",
  "Figure_S6_segmentation_accuracy.R",
  "Figure_S7_od_linearity.R",
  "Figure_S8_od_vs_cfu.R",
  "Table_S1_window_sensitivity.R",
  "Table_S2_segmentation_accuracy.R",
  "Table_S3_growth_rate_statistics.R",
  "cell_volume.R"
)

for (s in scripts) {
  cat("\n####", s, "\n")
  # The Figure 2 script prints each panel; outside RStudio that needs an open
  # graphics device, so the previews go to a temporary PNG file
  if (!interactive()) {
    png(tempfile(fileext = ".png"))
  }
  source(file.path("R", s), local = new.env())
  if (!interactive()) invisible(dev.off())
}
