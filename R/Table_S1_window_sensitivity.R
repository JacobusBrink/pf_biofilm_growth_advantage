# Introduction -----------------------------------------------------------

# Supplementary Table S1: sensitivity of the specific growth rates, and of
# the biofilm vs planktonic comparison, to the fitting window.
#   Scenarios: the selected windows (as in the main text), and common
#   windows of 2-6, 2-7, 3-7 and 2-8 h applied to both lifestyles.
#   For each scenario: mu (mean +/- SD of 3 replicates) per lifestyle, and
#   the biofilm vs planktonic % difference from the same model as Table S2
#   (gls on log mu, separate variance per group, Bonferroni over the five
#   comparisons within the scenario).
# The selected-window rows are identical to Table S2, section C.
#
# Outputs: output/Table_S1_window_sensitivity.xlsx (formatted table)
#          output/Table_S1_window_sensitivity_details.xlsx

source("R/00_setup.R")

suppressMessages({
  library(nlme) # gls (used inside fit_group_model)
  library(emmeans) # contrasts
})

# Growth rates for every scenario ----------------------------------------------

growth <- replicate_means()

scenarios <- tribble(
  ~scenario          , ~start , ~end ,
  "2-6 h"            ,      2 ,    6 ,
  "2-7 h"            ,      2 ,    7 ,
  "3-7 h"            ,      3 ,    7 ,
  "2-8 h (original)" ,      2 ,    8
)

rates <- bind_rows(
  fit_growth_rates(growth) |> mutate(scenario = "Selected windows"),
  pmap(scenarios, \(scenario, start, end) {
    fit_growth_rates(growth, windows |> mutate(start = !!start, end = !!end)) |>
      mutate(scenario = !!scenario)
  }) |>
    list_rbind()
)

scenario_levels <- c("Selected windows", scenarios$scenario)

# Biofilm vs planktonic within each scenario ------------------------------------

contrasts <- rates |>
  nest(.by = scenario) |>
  mutate(res = map(data, \(d) lifestyle_contrasts(fit_group_model(d)))) |>
  select(-data) |>
  unnest(res)

# Table ------------------------------------------------------------------------------

mu_cells <- rates |>
  summarise(
    text = sprintf("%.2f ± %.2f", mean(mu), sd(mu)),
    .by = c(scenario, lifestyle, condition)
  ) |>
  pivot_wider(names_from = lifestyle, values_from = text)

window_text <- rates |>
  distinct(scenario, lifestyle, condition, start, end) |>
  mutate(w = paste0(start, "-", end, " h")) |>
  select(-start, -end) |>
  pivot_wider(names_from = lifestyle, values_from = w) |>
  mutate(
    `Window (biofilm / planktonic)` = if_else(
      is.na(Planktonic) | Planktonic == Biofilm,
      Biofilm,
      paste(Biofilm, "/", Planktonic)
    )
  ) |>
  select(scenario, condition, `Window (biofilm / planktonic)`)

table_s1 <- mu_cells |>
  left_join(window_text, by = c("scenario", "condition")) |>
  left_join(
    contrasts |>
      transmute(
        scenario,
        condition = factor(contrast, levels = condition_levels),
        Difference = fmt_pct_ci(ratio, lower, upper),
        p = fmt_p(p_bonferroni)
      ),
    by = c("scenario", "condition")
  ) |>
  mutate(
    across(c(Planktonic, Difference, p), \(x) replace_na(x, "-")),
    scenario = factor(scenario, levels = scenario_levels)
  ) |>
  arrange(condition, scenario) |>
  transmute(
    Condition = as.character(condition),
    Scenario = as.character(scenario),
    `Window (biofilm / planktonic)`,
    `Biofilm µ (h⁻¹)` = Biofilm,
    `Planktonic µ (h⁻¹)` = Planktonic,
    `Difference (95% CI)` = Difference,
    `Adjusted p` = p
  )

options(width = 200)
print(table_s1, n = Inf)

write_xlsx(
  table_s1,
  file.path(out_dir, "supplementary_tables", "Table_S1_window_sensitivity.xlsx")
)

write_xlsx(
  list(
    growth_rates_per_replicate = rates |>
      mutate(condition = as.character(condition)),
    contrasts = contrasts
  ),
  file.path(out_dir, "stats", "Table_S1_window_sensitivity_details.xlsx")
)

cat("\nSaved Table S1 to", out_dir, "\n")
