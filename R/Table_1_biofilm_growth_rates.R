# Introduction -----------------------------------------------------------

# Table 1: specific growth rate (mu), R2 and doubling time (td) of biofilm
# populations, from segmented cell counts.
#
# mu: slope of ln(cells per colony) against time over 3-7 h, fitted for each
#     biological replicate to the mean of its 10 colonies; the table gives
#     the mean +/- SD of the 3 replicates.
# R2: mean R2 of the three replicate fits.
# td: doubling time, ln(2) / mean mu.
#
# Output: output/Table_1_biofilm_growth_rates.xlsx

source("R/00_setup.R")

# Growth rate per biological replicate ------------------------------------

rates <- replicate_means() |>
  filter(lifestyle == "Biofilm") |>
  fit_growth_rates()

# Summary per condition ----------------------------------------------------

table_1 <- rates |>
  summarise(
    window = paste0(first(start), "-", first(end), " h"),
    mu_mean = mean(mu),
    mu_sd = sd(mu),
    r2_mean = mean(r2),
    .by = condition
  ) |>
  mutate(td = log(2) / mu_mean) |>
  arrange(condition)

print(table_1 |> mutate(across(where(is.double), \(x) round(x, 3))))

# Formatted version (3 decimals for mu, as in the Results text)
table_1_formatted <- table_1 |>
  transmute(
    Condition = as.character(condition),
    Window = window,
    `µ (h⁻¹, mean ± SD)` = sprintf("%.3f ± %.3f", mu_mean, mu_sd),
    `R² (mean)` = sprintf("%.3f", r2_mean),
    `td (h)` = sprintf("%.2f", td)
  )

write_xlsx(
  list(
    table_1 = table_1_formatted,
    values = table_1 |> mutate(condition = as.character(condition)),
    per_replicate = rates |> mutate(condition = as.character(condition))
  ),
  file.path(out_dir, "tables", "Table_1_biofilm_growth_rates.xlsx")
)

cat("\nSaved Table 1 to", out_dir, "\n")
