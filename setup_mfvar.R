# -------------------------
# data settings
obs_start <- "2000-01-01" # start of the sample
obs_end <- "2025-10-01" # end of the sample
rt_date <- as.Date(Sys.time()) # real time vintage

# variabe codes and transformations (0 = level, 1 = log, 2 = diff, 3 = log-diff, 4 = annualized log-diff)
output_var <- c("GDPC1") # charts will be produced for this variable (can be multiple variables)

var_qd <- c("GDPC1" = 4) # quarterly variables and transformation
var_md <- c("INDPRO" = 3, "FEDFUNDS" = 3) # monthly variables and transformation

fhorz <- 8 # forecast horizon in months 
grid_tau <- seq(0.05,0.95,by = 0.01)
dir_main <- "results" # output directory for draws
dir_plots <- "plots" # output directory for plots

# -------------------------
# settings of the algorithm
nburn <- 1000 # burnin observations
nsave <- 3000 # number of saved draws
nthin <- 1 # thinfactor

run_mean <- "bart-hs" # c("bart-proj", "bart-hs", "lin")
outlier <- TRUE # outlier specification
p <- 12 # number of lags

source("est_mfvar.R")