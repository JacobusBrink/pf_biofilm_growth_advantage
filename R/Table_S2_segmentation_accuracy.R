# Introduction -----------------------------------------------------------

# Supplementary Table S2: segmentation accuracy against blinded manual
# annotation (Reviewer Comment 5). Also the source of every number in
# Supplementary Section S1 and in the Results validation paragraph.
#
# Input: the scored validation (validation/final/results/), made in Python
# by validation/final/score_final.py, which compares the manual clicks and
# outlines with the StarDist label images (those need the image drive, so
# the scoring itself is not repeated here).
#   partA_per_image.csv - 30 stacks, one click per cell:
#     true positive  = object matched by one click (or the first of several)
#     false negative = click without an object, or each extra click in a
#                      merged object (two cells joined into one)
#     false positive = object without a click
#   partB_per_cell.csv  - manual outlines of cells on one 8 h slice per
#     carbon source, matched one-to-one to segmented cross-sections (IoU)
#
# Pooled over the images in each group:
#   precision = TP / (TP + FP), recall = TP / (TP + FN),
#   F1 = 2 x precision x recall / (precision + recall),
#   count error = (automated - manual) / manual x 100
#
# Outputs: output/Table_S2_segmentation_accuracy.xlsx
#          output/Section_S1_segmentation_numbers.xlsx

source("R/00_setup.R")

res_dir <- "validation/final/results"

part_a <- read_csv(
  file.path(res_dir, "partA_per_image.csv"),
  show_col_types = FALSE
) |>
  mutate(
    carbon = str_to_title(carbon),
    caa = if_else(
      medium == "caa",
      "With Casamino Acids",
      "Without Casamino Acids"
    ),
    density = if_else(
      kind == "whole",
      "Low density (whole stacks)",
      "High density (crops)"
    )
  )

# Table S2 ------------------------------------------------------------------

pool <- function(d, group) {
  d |>
    summarise(
      Images = n(),
      `Manual cells` = sum(manual_cells),
      `Automated objects` = sum(automated_objects),
      tp = sum(true_positive),
      fp = sum(false_positive),
      fn = sum(false_negative),
      missed = sum(missed_cells),
      .by = all_of(group)
    ) |>
    mutate(
      Precision = tp / (tp + fp),
      Recall = tp / (tp + fn),
      F1 = 2 * Precision * Recall / (Precision + Recall),
      `Count error (%)` = (`Automated objects` - `Manual cells`) /
        `Manual cells` *
        100
    ) |>
    rename(Group = all_of(group))
}

table_s2 <- bind_rows(
  part_a |> mutate(all = "All images") |> pool("all"),
  pool(part_a, "density") |> arrange(desc(Group)),
  pool(part_a, "carbon") |> arrange(Group),
  pool(part_a, "caa") |> arrange(Group)
)

print(
  table_s2 |>
    mutate(
      across(c(Precision, Recall, F1), \(x) round(x, 2)),
      `Count error (%)` = round(`Count error (%)`, 1)
    )
)

write_xlsx(
  table_s2 |>
    select(-tp, -fp, -fn, -missed) |>
    mutate(
      across(c(Precision, Recall, F1), \(x) sprintf("%.2f", x)),
      `Count error (%)` = sprintf("%.1f", `Count error (%)`)
    ),
  file.path(
    out_dir,
    "supplementary_tables",
    "Table_S2_segmentation_accuracy.xlsx"
  )
)

# Section S1 and Results numbers ------------------------------------------------

# Cells missed completely (click with no object at all; the other false
# negatives are merged cells)
missed <- tibble(
  missed_cells = sum(part_a$missed_cells),
  manual_cells = sum(part_a$manual_cells)
)

# Bias in the specific growth rate: if the undercount changed with time it
# would change mu. ln(automated / manual) regressed on time for the 15
# stacks with CAA in the 3-7 h window, common slope and a separate intercept
# per carbon source; the slope (h-1) is the error added to mu.
bias_data <- part_a |>
  filter(medium == "caa", timepoint_h >= 3, timepoint_h <= 7)
bias_model <- lm(
  log(automated_objects / manual_cells) ~ carbon + timepoint_h,
  data = bias_data
)
mu_bias <- tibble(
  images = nrow(bias_data),
  mu_bias = coef(bias_model)[["timepoint_h"]],
  ci_low = confint(bias_model)["timepoint_h", 1],
  ci_high = confint(bias_model)["timepoint_h", 2]
)

# Cell shape at 8 h (manual outlines): matched cells have IoU >= 0.5
part_b_images <- read_csv(
  file.path(res_dir, "partB_per_image.csv"),
  show_col_types = FALSE
)
part_b_cells <- read_csv(
  file.path(res_dir, "partB_per_cell.csv"),
  show_col_types = FALSE
) |>
  filter(iou >= 0.5)
shape <- tibble(
  outlined_cells = sum(part_b_images$manual_cells),
  matched_iou50 = nrow(part_b_cells),
  matched_pct = nrow(part_b_cells) / sum(part_b_images$manual_cells) * 100,
  mean_iou_matched = mean(part_b_cells$iou),
  median_area_ratio = median(
    part_b_cells$area_auto_um2 / part_b_cells$area_manual_um2
  ),
  median_centroid_offset_um = median(part_b_cells$centroid_offset_um)
)

# Nearest-neighbour distance in the cropped (dense) stacks
nearest_neighbour <- part_a |>
  filter(kind == "crop") |>
  summarise(
    nn_manual_um = mean(nn_manual_um),
    nn_automated_um = mean(nn_auto_um)
  )

cat("\n== Missed cells ==\n")
print(missed)
cat("\n== Bias in mu (h-1) ==\n")
print(mu_bias)
cat("\n== Cell shape at 8 h ==\n")
print(shape)
cat("\n== Nearest-neighbour distance, crops ==\n")
print(nearest_neighbour)

write_xlsx(
  list(
    missed = missed,
    mu_bias = mu_bias,
    shape_8h = shape,
    nearest_neighbour = nearest_neighbour
  ),
  file.path(
    out_dir,
    "stats",
    "Section_S1_segmentation_numbers.xlsx"
  )
)

cat("\nSaved Table S2 and Section S1 numbers to", out_dir, "\n")
