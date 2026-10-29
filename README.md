<p align="center">
  <img src="assets/lab_logo.png" alt="Biofilm Engineering Laboratory Logo" width="160">
</p>

<h1 align="center">
  Welcome to the Biofilm Engineering Laboratory repository for Single-cell counting reveals a modest, nutrient-dependent growth advantage of early <i>Pseudomonas fluorescens</i> biofilms over planktonic cells
</h1>

Analysis code and data for comparing the specific growth rate of *Pseudomonas fluorescens* growing as biofilms versus planktonic cultures, on glucose, glycerol and succinate, with and without casamino acids (CAA).

## Contents

- `R/`: one script per figure and table of the manuscript. `00_setup.R` holds the shared settings (packages, theme, exponential-phase windows, fitting and statistics), and `run_all.R` rebuilds every output. See `R/README.md` for details.
- `input/`: raw data. Biofilm cell counts per colony (BiofilmQ, `input/biofilm/`), planktonic OD600 (`input/planktonic/`), and OD600 linearity and OD vs CFU calibration data (`input/od_vs_cfu/`).
- `output/`: figures, tables, supplementary material and statistics produced by the scripts.
- `validation/`: segmentation accuracy check (manual cell annotations scored against the StarDist segmentation in Python).

## Usage

Run from the repository root in R. The package versions are recorded with renv, so first install them:

```r
renv::restore()
source("R/run_all.R")
```

## License

The code (`R/` and the Python and ImageJ scripts in `validation/`) is released under the [MIT License](LICENSE). The data, images and results (`input/`, `output/` and `validation/`) are released under [CC BY 4.0](LICENSE-data).
