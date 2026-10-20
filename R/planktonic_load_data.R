# Load required packages -------------------------------------------------

library(readxl) # read Excel files easily
library(dplyr) # data manipulation and wrangling

# Import raw data --------------------------------------------------------

# Read data using readxl
data <- read_excel(
  "input/planktonic/planktonic_glucose_glycerol_succinate.xlsx",
  na = "NA"
)

# Process and tidy data --------------------------------------------------

# Clean and transform the data
clean_data <- data |>

  # convert carbon and biological_replicate to factors
  mutate(
    carbon = as.factor(carbon),
    biological_replicate = as.factor(biological_replicate)
  ) |>

  # rename carbon codes
  mutate(
    carbon = case_when(
      carbon == "glucose" ~ "M9-Glucose",
      carbon == "glycerol" ~ "M9-Glycerol",
      carbon == "succinate" ~ "M9-Succinate",
      carbon == "glucose_control" ~ "M9-Glucose control",
      carbon == "succinate_control" ~ "M9-Succinate control",
      TRUE ~ carbon
    )
  ) |>

  # remove od and dilution factor columns
  select(
    -od,
    -dilution_factor
  ) |>

  # return as tibble
  as_tibble()

# Export cleaned data ----------------------------------------------------

# Save in RDS format
saveRDS(clean_data, file = "input/clean_planktonic_carbon_data.rds")

# Write to an excel file
writexl::write_xlsx(
  clean_data,
  path = "input/clean_planktonic_carbon_data.xlsx"
)
