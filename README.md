## Description
> These codes come without technical support of any kind. The code is free to use, provided that the paper is cited properly.

Codes based on M. Marcellino & M. Pfarrhofer (2025): "[Nonparametric Mixed Frequency Monitoring Macro-at-Risk](https://doi.org/10.1016/j.econlet.2025.112498)" _Economics Letters_ *255* 112498, using state space MF-VARs estimated with precision sampling (including MF-BART). See also [mf-bavart](https://github.com/mpfarrho/mf-bavart) for related work.

Estimation is integrated with the [alfred](https://cran.r-project.org/web/packages/alfred/index.html) R-package to access real-time data from [alfred.stlouisfed.org](https://alfred.stlouisfed.org) for the United States.

## Source files
- `setup_mfvar` is the main estimation file which downloads and transforms data (requires user input)
- `est_mfvar` estimates the model and produces outputs (is sourced automatically from `setup_mfvar`)
- `utils` contains several helper functions (is sourced automatically from `est_mfvar`)

### Estimation options in `setup_mfvar`
- Data settings are described in the file
- `run_mean` refers to the conditional mean estimation and approximation, choose from:
  - `bart-hs` estimates a BART model and uses the shrinkage-based linear approximation
  - `bart-proj` estimates a BART model and uses the projection-based approximation
  - `lin` estimates a linear BVAR model with horseshoe prior
 - `outlier` refers to whether a homoskedastic (when `FALSE`) or outlier (when `TRUE`) specification is estimated

Model outputs (MCMC draws for the data) are saved in the directory `results` and associated charts for variables defined in `output_var` are collected in the directory `plots`.
