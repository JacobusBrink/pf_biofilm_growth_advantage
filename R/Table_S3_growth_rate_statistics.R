# Introduction -----------------------------------------------------------

# Supplementary Table S3: statistical analysis of the specific growth rates
# (Reviewer Comment 10). Also the source of the statistics quoted in the
# Results (lifestyle comparison) and in the response to Comment 10.
#
# Data: one growth rate per biological replicate (n = 3 per group).
#   Biofilm:    cell counts, 3-7 h
#   Planktonic: OD600 x dilution factor, exponential phase (00_setup.R)
#
# All models use log(mu), so differences are ratios (% differences), with
# generalised least squares (gls) and a separate variance for each group.
#   A  With CAA, all three carbon sources: log(mu) ~ lifestyle * carbon.
#      Joint F tests: does mu differ between lifestyles, between carbon
#      sources, and does the lifestyle difference depend on carbon source?
#   B  Glucose and succinate, with and without CAA:
#      log(mu) ~ lifestyle * carbon * caa (glycerol has no planktonic data
#      without CAA).
#   C  Biofilm vs planktonic % difference: averaged over carbon sources
#      (from model A) and for each condition (Bonferroni over 5).
#   D  Change in the biofilm/planktonic ratio when CAA was removed
#      (Bonferroni over 2).
#   E  Carbon source comparisons within each lifestyle, with CAA
#      (Bonferroni over 6).
#
# Outputs: output/Table_S3_growth_rate_statistics.xlsx (formatted table)
#          output/Table_S3_growth_rate_statistics_details.xlsx (all
#          estimates, SE, df, unadjusted p)

source("R/00_setup.R")

suppressMessages({
  library(nlme) # gls with group-specific variances
  library(emmeans) # estimated marginal means, contrasts and joint tests
})

# Growth rate per biological replicate ---------------------------------------

rates <- fit_growth_rates(replicate_means())

# C and D: biofilm vs planktonic, and change on removal of CAA ------------------

fitted <- fit_group_model(rates)
gn <- levels(fitted$data$group)

family_c <- lifestyle_contrasts(fitted)

# D: (biofilm/planktonic without CAA) / (biofilm/planktonic with CAA)
weights_d <- map(c("Glucose", "Succinate"), \(cb) {
  contrast_weights(
    setNames(
      c(1, -1, -1, 1),
      paste0(
        cb,
        c(
          " -CAA.Biofilm",
          " -CAA.Planktonic",
          " +CAA.Biofilm",
          " +CAA.Planktonic"
        )
      )
    ),
    gn
  )
}) |>
  set_names(c("Glucose", "Succinate"))
family_d <- tidy_contrasts(
  fitted$emm,
  weights_d,
  "Change in biofilm advantage on removal of CAA"
)

# A: with CAA, lifestyle x carbon source -----------------------------------------

d_model <- rates |>
  mutate(
    log_mu = log(mu),
    lifestyle = factor(lifestyle, levels = c("Planktonic", "Biofilm")),
    carbon = factor(carbon, levels = names(carbon_colours)),
    caa = factor(caa, levels = c("+CAA", "-CAA")),
    group = interaction(carbon, caa, lifestyle, drop = TRUE, sep = ".")
  )

d_a <- d_model |> filter(caa == "+CAA") |> droplevels()
model_a <- gls(
  log_mu ~ lifestyle * carbon,
  data = d_a,
  weights = varIdent(form = ~ 1 | group)
)

joint_a <- joint_tests(model_a, mode = "satterthwaite", data = d_a) |>
  as_tibble() |>
  mutate(model = "A: with CAA (lifestyle x carbon source)")

# Overall biofilm vs planktonic difference, averaged over carbon sources
# (meaningful because the lifestyle x carbon interaction is not significant)
lifestyle_overall <- emmeans(
  model_a,
  ~lifestyle,
  mode = "satterthwaite",
  data = d_a
) |>
  pairs(reverse = TRUE) |>
  summary(infer = TRUE) |>
  as_tibble() |>
  mutate(ratio = exp(estimate), lower = exp(lower.CL), upper = exp(upper.CL))

# E: carbon source comparisons within each lifestyle
family_e <- emmeans(
  model_a,
  ~ carbon | lifestyle,
  mode = "satterthwaite",
  data = d_a
) |>
  pairs(reverse = FALSE) |>
  summary(infer = TRUE, adjust = "none", by = NULL) |>
  as_tibble() |>
  mutate(
    ratio = exp(estimate),
    t_crit = qt(1 - 0.05 / (2 * n()), df),
    lower = exp(estimate - t_crit * SE),
    upper = exp(estimate + t_crit * SE),
    p_bonferroni = pmin(p.value * n(), 1)
  )

# B: glucose and succinate, lifestyle x carbon x CAA ------------------------------

d_b <- d_model |> filter(carbon %in% c("Glucose", "Succinate")) |> droplevels()
model_b <- gls(
  log_mu ~ lifestyle * carbon * caa,
  data = d_b,
  weights = varIdent(form = ~ 1 | group)
)

joint_b <- joint_tests(model_b, mode = "satterthwaite", data = d_b) |>
  as_tibble() |>
  mutate(model = "B: glucose and succinate, with and without CAA")

joint <- bind_rows(joint_a, joint_b) |> rename(term = `model term`)

# Console output -------------------------------------------------------------------

options(width = 200)
cat("\n== Joint F tests ==\n")
print(
  joint |>
    mutate(
      across(c(F.ratio, df2), \(x) round(x, 2)),
      p.value = signif(p.value, 2)
    ),
  n = Inf
)
cat("\n== Biofilm vs planktonic, averaged over carbon sources (with CAA) ==\n")
print(
  lifestyle_overall |>
    transmute(
      pct = (ratio - 1) * 100,
      lower = (lower - 1) * 100,
      upper = (upper - 1) * 100,
      df,
      p.value
    )
)
cat("\n== C and D ==\n")
print(
  bind_rows(family_c, family_d) |>
    transmute(
      family,
      contrast,
      pct = round((ratio - 1) * 100, 1),
      ci = sprintf("%+.0f to %+.0f%%", (lower - 1) * 100, (upper - 1) * 100),
      df = round(df, 1),
      p_adj = signif(p_bonferroni, 2)
    ),
  n = Inf
)

# Supplementary Table S3 -------------------------------------------------------------

term_label <- c(
  "lifestyle" = "Lifestyle",
  "carbon" = "Carbon source",
  "caa" = "CAA",
  "lifestyle:carbon" = "Lifestyle × carbon source",
  "lifestyle:caa" = "Lifestyle × CAA",
  "carbon:caa" = "Carbon source × CAA",
  "lifestyle:carbon:caa" = "Lifestyle × carbon source × CAA"
)

section_c <- "C: Biofilm vs planktonic with CAA, % difference (95% CI)"
section_c2 <- "C2: Biofilm vs planktonic without CAA, % difference (95% CI)"
section_d <- "D: Change in biofilm vs planktonic difference on removal of CAA (95% CI)"
section_e <- "E: Carbon source comparisons with CAA, % difference (95% CI)"

table_s3 <- bind_rows(
  joint |>
    transmute(
      section = model,
      comparison = term_label[term],
      estimate = sprintf("F(%d, %.1f) = %.1f", df1, df2, F.ratio),
      p = fmt_p(p.value)
    ),
  lifestyle_overall |>
    transmute(
      section = section_c,
      comparison = "All carbon sources (average)",
      estimate = fmt_pct_ci(ratio, lower, upper),
      p = fmt_p(p.value)
    ),
  family_c |>
    transmute(
      section = if_else(str_detect(contrast, "-CAA"), section_c2, section_c),
      comparison = contrast,
      estimate = fmt_pct_ci(ratio, lower, upper),
      p = fmt_p(p_bonferroni)
    ),
  family_d |>
    transmute(
      section = section_d,
      comparison = contrast,
      estimate = fmt_pct_ci(ratio, lower, upper),
      p = fmt_p(p_bonferroni)
    ),
  family_e |>
    transmute(
      section = section_e,
      comparison = paste0(
        lifestyle,
        ": ",
        str_replace(contrast, " - ", " vs ")
      ),
      estimate = fmt_pct_ci(ratio, lower, upper),
      p = fmt_p(p_bonferroni)
    )
) |>
  rename(
    `Term or comparison` = comparison,
    `Test statistic or estimate` = estimate
  )

print(table_s3, n = Inf)

write_xlsx(
  table_s3,
  file.path(
    out_dir,
    "supplementary_tables",
    "Table_S3_growth_rate_statistics.xlsx"
  )
)

write_xlsx(
  list(
    growth_rates_per_replicate = rates |>
      mutate(condition = as.character(condition)),
    joint_tests = joint,
    C_overall = lifestyle_overall,
    C_D_contrasts = bind_rows(family_c, family_d),
    E_carbon_comparisons = family_e |>
      mutate(lifestyle = as.character(lifestyle))
  ),
  file.path(out_dir, "stats", "Table_S3_growth_rate_statistics_details.xlsx")
)

cat("\nSaved Table S3 to", out_dir, "\n")
