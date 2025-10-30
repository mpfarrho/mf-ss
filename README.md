## Description
Codes based on M. Marcellino & M. Pfarrhofer (2025): "[Nonparametric Mixed Frequency Monitoring Macro-at-Risk](https://doi.org/10.1016/j.econlet.2025.112498)" _Economics Letters_ *255* 112498, using state space MF-VARs estimated with precision sampling (including MF-BART). See also [mf-bavart](https://github.com/mpfarrho/mf-bavart) for related work.

Integrated with the [alfred](https://cran.r-project.org/web/packages/alfred/index.html) R-package to access real-time data from [alfred.stlouisfed.org](https://alfred.stlouisfed.org).

- > [!CAUTION]
> These codes come without technical support of any kind. The code is free to use, provided that the paper is cited properly.

## Source files
- `setup_mfvar` is the main estimation file which downloads and transforms data
- `est_mfvar` estimates the model and produces outputs
- `utils` contains several helper functions

### Estimation options in `setup_mfvar`
- Data settings are described in the file
- `run_mean` refers to the conditional mean estimation and approximation, choose from:
  - `bart-hs` estimates a BART model and uses the shrinkage-based linear approximation
  - `bart-proj` estimates a BART model and uses the projection-based approximation
  - `lin` estimates a linear BVAR model with horseshoe prior
 - `outlier` refers to whether a homoskedastic or outlier specification is estimated

Model outputs are saved in the directory `results` and associated charts are collected in the directory `plots`.
