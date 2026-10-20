# R: one script per figure and table

Each script makes one figure or table of the manuscript, including its analysis and statistics, and saves the result to `output/`. Run the scripts from the repository root (the folder that contains `R/` and `input/`).

- `00_setup.R` is sourced at the top of every script. It holds what must be identical everywhere: packages, theme, colours, the exponential-phase **windows**, the functions that read the raw data, the log-linear fit, and the biofilm vs planktonic model.
- `run_all.R` runs every script and rebuilds `output/` from scratch: `source("R/run_all.R")`, or `Rscript R/run_all.R` from a terminal.
- `biofilm_growth_load_data.R` and `planktonic_load_data.R` make the cleaned data files in `input/` (`clean_biofilm_data.rds`, `clean_planktonic_carbon_data.rds`). They are not part of `run_all.R`.

## Scripts, outputs and where the numbers are used

| Item | Script | Output (in `output/`) | Also the source of |
|---|---|---|---|
| Figure 1 | none (microscopy images) | | |
| Figure 2 | `Figure_2_biofilm_architecture.R` | `figures/Figure_2_biofilm_architecture.png`, `stats/Figure_2_biofilm_architecture_stats.xlsx` | Results: structural properties (mean ± SD, contrasts) |
| Figure 3 | `Figure_3_biofilm_growth.R` | `figures/Figure_3_biofilm_growth.png` | |
| Figure 4 | `Figure_4_planktonic_growth.R` | `figures/Figure_4_planktonic_growth.png` | |
| Table 1 | `Table_1_biofilm_growth_rates.R` | `tables/Table_1_biofilm_growth_rates.xlsx` | Results: biofilm kinetics |
| Table 2 | `Table_2_planktonic_growth_rates.R` | `tables/Table_2_planktonic_growth_rates.xlsx` | Results: planktonic kinetics |
| Figure S1 | `Figure_S1_planktonic_replicate_fits.R` | `supplementary_figures/*.png` | |
| Figure S2 | `Figure_S2_biofilm_replicate_fits.R` | `supplementary_figures/*.png` | |
| Figure S3 | `Figure_S3_hourly_growth_rates.R` | `supplementary_figures/*.png`, `stats/*.xlsx` | Methods: window selection |
| Figure S4 | `Figure_S4_residuals.R` | `supplementary_figures/*.png` | |
| Figure S5 | none (overlay images, made in Fiji) | | |
| Figure S6 | `Figure_S6_segmentation_accuracy.R` | `supplementary_figures/*.png` | |
| Figure S7 | `Figure_S7_od_linearity.R` | `supplementary_figures/*.png`, `stats/*.xlsx` | Results: log-log slope 0.94, deviation at OD 0.3–0.4 |
| Figure S8 | `Figure_S8_od_vs_cfu.R` | `supplementary_figures/*.png`, `stats/*.xlsx` | Results: CFU per OD600 over time, µ_CFU − µ_OD |
| Table S1 | `Table_S1_window_sensitivity.R` | `supplementary_tables/*.xlsx`, `stats/*_details.xlsx` | Results: window sensitivity |
| Table S2 | `Table_S2_segmentation_accuracy.R` | `supplementary_tables/*.xlsx`, `stats/Section_S1_segmentation_numbers.xlsx` | Supplementary Section S1; Results: validation paragraph |
| Table S3 | `Table_S3_growth_rate_statistics.R` | `supplementary_tables/*.xlsx`, `stats/*_details.xlsx` | Results: lifestyle and carbon comparison, CAA effect |
| Table S4 | none (Python: `validation/final/table_s4.py`) | `segmentation_validation/Table_S4_segmentation_validation.csv` | |
| Text | `cell_volume.R` | `stats/Text_cell_volume.xlsx` | Results: cell volume |

## Input data

| Data | File |
|---|---|
| Biofilm cell counts (re-segmented, BiofilmQ, one CSV per colony) | `input/biofilm/*.csv` |
| Biofilm architecture at 8 h (Figure 2) | `input/clean_biofilm_data.rds` (made by `R/biofilm_growth_load_data.R`) |
| Planktonic OD600 | `input/planktonic/planktonic_glucose_glycerol_succinate.xlsx` |
| OD600 dilution series | `input/od_vs_cfu/od_linearity_template.xlsx` |
| Paired OD600 and CFU | `input/od_vs_cfu/cfu_od_template.xlsx` |
| Segmentation validation (scored in Python) | `validation/final/results/*.csv` |

## Notes

- **Growth rates.** µ is fitted per biological replicate, to the mean of its technical replicates (10 colonies or 3 flasks) at each timepoint, over the window in `00_setup.R`. Tables, figures and statistics all use the same function, so changing a window in `00_setup.R` changes everything consistently.
- **Small p-value differences.** With n = 3 and one group with almost no spread (glucose without CAA, biofilm: SD 0.001 h⁻¹), the Satterthwaite degrees of freedom are numerically fragile: rounding differences of 10⁻¹⁶ can move some p-values in the third decimal between machines. No conclusion changes.
- **Figure 2 outside RStudio.** The script prints each panel as it goes. With `Rscript` that needs a graphics device, which `run_all.R` opens for you.
- **Segmentation validation.** Matching clicks and outlines to the label images needs the original image drive and is done in Python (`validation/final/score_final.py`). The R scripts start from its results.
