# Load required packages -------------------------------------------------

library(tidyverse) # data manipulation
library(fs) # file system utilities
library(janitor) # clean column names
library(writexl) # write Excel files

# Import raw data --------------------------------------------------------

# Get list of all CSV files in the input/biofilm folder
csv_files <-
  dir_ls(
    path = "input/biofilm"
  )

# Read files, keep every second row, and clean column names
data <- read_csv(csv_files) |>
  slice(seq(2, n(), by = 2)) |> # take rows 2,4,6,... to skip metadata rows
  clean_names() # make column names snake_case

# Process and tidy data --------------------------------------------------

# Select columns, extract info from filename, and tidy data
biofilm_data <- data |>
  select(
    filename,
    ncells,
    biofilm_height,
    biofilm_base_area,
    biofilm_volume,
    distance_to_nearest_neighbor_mean,
    architecture_nematic_order_parameter_range30_mean,
    distance_inter_cell_spacing_mean_range20_mean
  ) |>

  # Split full filename into components
  separate(
    filename,
    into = c(
      "sample_id",
      "timepoint",
      "position",
      "channel",
      "frame"
    ),
    sep = "_"
  ) |>

  # Keep the following columns
  select(
    sample_id,
    ncells,
    biofilm_height,
    biofilm_base_area,
    biofilm_volume,
    distance_to_nearest_neighbor_mean,
    architecture_nematic_order_parameter_range30_mean,
    distance_inter_cell_spacing_mean_range20_mean
  ) |>

  # Extract media, timepoint, and replicates from sample_id
  extract(
    sample_id,
    into = c(
      "media",
      "timepoint",
      "biological_replicate",
      "technical_replicate"
    ),
    regex = "^([A-Za-z]+)([0-9]+h)([A-Za-z]+[0-9]+)([A-Za-z]+[0-9]+)$",
    remove = TRUE,
    convert = FALSE
  ) |>

  # Recode variables and convert numeric fields
  mutate(
    media = case_when(
      media == "glucosecaa" ~ "Glucose (CAA 120 mM)",
      media == "glycerolcaa" ~ "Glycerol (CAA 120 mM)",
      media == "succinatecaa" ~ "Succinate (CAA 120 mM)",
      media == "glucosecontrol" ~ "Glucose (no CAA 160 mM)",
      media == "glycerolcontrol" ~ "Glycerol (no CAA 160 mM)",
      media == "succinatecontrol" ~ "Succinate (no CAA 160 mM)",
      # media == "glucosetwozeromm" ~ "Glucose (no CAA 120 mM)",
      # media == "glucosethreemm" ~ "Glucose (no CAA 20 mM)",
      # media == "glucosethirteenmm" ~ "Glucose (no CAA 80 mM)",
      # media == "glucosethreetotalmm" ~ "Glucose (no CAA 3 mM)",
      TRUE ~ media
    ),
    timepoint = str_remove(timepoint, "h") |> as.numeric(), # remove 'h' and convert
    biological_replicate = str_extract(biological_replicate, "[0-9]+") |>
      as.numeric(),
    technical_replicate = str_extract(technical_replicate, "[0-9]+") |>
      as.numeric()
  ) |>

  # ensure final format is a tibble
  as_tibble()

# Export cleaned data ----------------------------------------------------

# Save in RDS format
saveRDS(biofilm_data, file = "input/clean_biofilm_data.rds")

# Write to an excel file
writexl::write_xlsx(biofilm_data, path = "input/clean_biofilm_data.xlsx")
