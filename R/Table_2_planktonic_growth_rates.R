# Introduction -----------------------------------------------------------

# Table 2: specific growth rate (mu), R2 and doubling time (td) of planktonic
# populations, from OD600.
#
# mu: slope of ln(OD600 x dilution factor) against time over the exponential
#     phase (with CAA: glucose and glycerol 2-5 h, succinate 2-6 h; without
#     CAA: 3-7 h), fitted for each biological replicate to the mean of its 3
#     flasks; the table gives the mean +/- SD of the 3 replicates.
# R2: mean R2 of the three replicate fits.
# td: doubling time, ln(2) / mean mu.
# No planktonic data were recorded for glycerol without CAA.
#
# Output: output/Table_2_planktonic_growth_rates.xlsx

source("R/00_setup.R")

# Growth rate per biological replicate ------------------------------------

rates <- replicate_means() |>
  filter(lifestyle == "Planktonic") |>
  fit_growth_rates()

# Summary per condition ----------------------------------------------------

table_2 <- rates |>
  summarise(
    window = paste0(first(start), "-", first(end), " h"),
    mu_mean = mean(mu),
    mu_sd = sd(mu),
    r2_mean = mean(r2),
    .by = condition
  ) |>
  mutate(td = log(2) / mu_mean) |>
  arrange(condition)

print(table_2 |> mutate(across(where(is.double), \(x) round(x, 3))))

# Formatted version (3 decimals for mu, as in the Results text)
table_2_formatted <- table_2 |>
  transmute(
    Condition = as.character(condition),
    Window = window,
    `µ (h⁻¹, mean ± SD)` = sprintf("%.3f ± %.3f", mu_mean, mu_sd),
    `R² (mean)` = sprintf("%.3f", r2_mean),
    `td (h)` = sprintf("%.2f", td)
  )

write_xlsx(
  list(
    table_2 = table_2_formatted,
    values = table_2 |> mutate(condition = as.character(condition)),
    per_replicate = rates |> mutate(condition = as.character(condition))
  ),
  file.path(out_dir, "tables", "Table_2_planktonic_growth_rates.xlsx")
)

cat("\nSaved Table 2 to", out_dir, "\n")
