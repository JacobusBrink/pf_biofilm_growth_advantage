# Introduction -----------------------------------------------------------

# Mean volume per segmented biofilm cell.
# Source of the Results sentence "In biofilms
# with CAA, the mean volume per segmented cell was 1.4-1.6 um3 between 3 and
# 7 h and differed by less than 3% between carbon sources."
#   - per colony: BiofilmQ Shape_Volume_mean (mean volume of the segmented
#     cells; equal to Biofilm_Volume / Ncells)
#   - mean of the 10 colonies per biological replicate and timepoint, then
#     mean +/- SD of the 3 biological replicates (as for the growth rates)
#   - 3-7 h, the biofilm fitting window
#
# Output: output/Text_cell_volume.xlsx

source("R/00_setup.R")

# Data ----------------------------------------------------------------------

colonies <- list.files(
  "input/biofilm",
  pattern = "[.]csv$",
  full.names = TRUE
) |>
  map(\(f) {
    # suppressWarnings: BiofilmQ rows have uneven lengths (see 00_setup.R)
    x <- suppressWarnings(read_csv(
      f,
      col_types = cols(.default = "c"),
      n_max = 2,
      col_select = Shape_Volume_mean
    ))
    tibble(file = basename(f), cell_volume = as.numeric(x$Shape_Volume_mean[2]))
  }) |>
  list_rbind() |>
  extract(
    file,
    c("media", "time", "bio_rep", "tech_rep"),
    "^([a-z]+)_([0-9]+)h_br([0-9])_tr([0-9]+)[.]csv$",
    convert = TRUE
  ) |>
  mutate(
    carbon = str_to_title(str_remove(media, "caa|control")),
    caa = if_else(str_detect(media, "control"), "-CAA", "+CAA")
  ) |>
  filter(time >= 3, time <= 7)

# Volume per cell -------------------------------------------------------------

replicates <- colonies |>
  summarise(
    cell_volume = mean(cell_volume),
    .by = c(caa, carbon, bio_rep, time)
  )

# Per condition and timepoint (mean +/- SD of 3 biological replicates)
by_time <- replicates |>
  summarise(
    mean = mean(cell_volume),
    sd = sd(cell_volume),
    .by = c(caa, carbon, time)
  )

# Per condition over 3-7 h
by_condition <- replicates |>
  summarise(cell_volume = mean(cell_volume), .by = c(caa, carbon, bio_rep)) |>
  summarise(
    mean = mean(cell_volume),
    sd = sd(cell_volume),
    .by = c(caa, carbon)
  )

# Numbers quoted in the text: range of the timepoint means and the largest
# difference between carbon sources (%)
text_numbers <- by_time |>
  summarise(min = min(mean), max = max(mean), .by = caa) |>
  left_join(
    by_condition |>
      summarise(carbon_diff_pct = 100 * (max(mean) / min(mean) - 1), .by = caa),
    by = "caa"
  )

options(width = 200)
cat("\n== Volume per cell (um3), mean of 3 biological replicates ==\n")
print(
  by_time |>
    mutate(mean = round(mean, 2)) |>
    select(-sd) |>
    pivot_wider(names_from = time, values_from = mean),
  n = Inf
)
cat("\n== 3-7 h per condition ==\n")
print(by_condition |> mutate(across(c(mean, sd), \(x) round(x, 2))), n = Inf)
cat("\n== Text numbers ==\n")
print(text_numbers |> mutate(across(where(is.numeric), \(x) round(x, 2))))

write_xlsx(
  list(
    text_numbers = text_numbers,
    by_condition = by_condition,
    by_time = by_time,
    replicates = replicates
  ),
  file.path(out_dir, "stats", "Text_cell_volume.xlsx")
)

cat("\nSaved Text_cell_volume.xlsx to", out_dir, "\n")
