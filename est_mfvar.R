library(Matrix)
library(MASS)
library(dbarts)
library(mvnfast)
library(invgamma)

library(zoo)
library(lubridate)

# data management
library(dplyr)
library(tidyr)
library(alfred) # load dataset from FRED

source("utils.R")
run_id <- paste0("rt",rt_date,"_mfvar_",run_mean,ifelse(outlier,"_o","hom"))
dir.create(dir_main, showWarnings = FALSE)
dir.create(dir_plots, showWarnings = FALSE)

# ------------------------------------------------------------------------------------
# transformations
obs_start_2 <- as.numeric(unlist(strsplit(substr(obs_start,1,7),"-")))
variables_tfm <- c(var_qd,var_md)
variables <- names(variables_tfm)
run_mbyq <- as.Date.mbyq(obs_end)

freq <- c(rep(4,length(var_qd)),rep(12,length(var_md)))

out <- lapply(variables, get_alfred_series,
              observation_start = obs_start, observation_end = obs_end,
              realtime_start = rt_date, realtime_end = rt_date)
alfred_to_ts <- function(x, freq){
  ts(x[,3],start=obs_start_2,frequency=freq)
}

mf_list <- mapply(alfred_to_ts, x = out, freq = freq)
names(mf_list) <- variables

df_raw <- data.frame(date = as.character(seq(as.Date(obs_start), as.Date(obs_end), by = "month")))
for(i in 1:length(mf_list)){
  lab_i <- names(mf_list)[i]
  y_i <- mf_list[[i]]
  if(variables_tfm[i] == 1){
    y_i <- 100 * log(y_i)
  }else if(variables_tfm[i] == 2){
    y_i <- diff(y_i)
  }else if(variables_tfm[i] == 3){
    y_i <- diff(100 * log(y_i))
  }else if(variables_tfm[i] == 4){
    y_i <- diff(freq[i] * 100 * log(y_i))
  }
  
  # adjust for matching periods
  date_tmp <- as.Date.ts(y_i)
  if(freq[i] == 4) date_tmp <- date_tmp %m+% months(2) # end of quarter observation for LF variable
  
  df_tmp <- data.frame(date = as.character(date_tmp), variable = y_i)
  colnames(df_tmp) <- c("date",lab_i)
  df_raw <- left_join(df_raw, df_tmp, by = "date")
}
Yraw <- ts(df_raw[,-1], start = obs_start_2, frequency = 12)
rownames(Yraw) <- as.character(as.Date.ts(Yraw))

datelabs <- rownames(Yraw)
varlabs <- gsub("_","",colnames(Yraw))
names(var_md) <- gsub("_","",names(var_md))
colnames(Yraw) <- varlabs

var_all <- variables
max_T <- nrow(Yraw)

# --------------------------------------------------------------------------------------------------
# other necessary adjustments
if(run_mbyq != 3){
  miss_mq <- 3 - (nrow(Yraw) %% 3)
  Yraw <- rbind(Yraw, matrix(NA, nrow = miss_mq, ncol = ncol(Yraw)))
  rownames(Yraw) <- c(datelabs,
                      as.character(seq(as.Date(max(datelabs)),as.Date("2100-12-01"), by = "month")[2:(miss_mq+1)]))
  datelabs_f <- as.character(seq(as.Date(max(rownames(Yraw))),as.Date("2100-12-01"), by = "month"))[2:(fhorz+1)]
}else{
  miss_mq <- 0
  datelabs_f <- as.character(seq(as.Date(max(rownames(Yraw))),as.Date("2100-12-01"), by = "month"))[2:(fhorz+1)]
}

n_lf <- length(var_qd)
n_hf <- length(var_md)
varlabs <- colnames(Yraw) # variable names, needed for selecting positions of temporal loadings
itr <- rep("D", n_lf) # triangular aggregation scheme

Yraw_orig <- Yraw
datelabs <- rownames(Yraw)
Ymu <- apply(Yraw, 2, mean, na.rm = TRUE)
Ysd <- apply(Yraw, 2, sd, na.rm = TRUE)
Yraw <- scale(Yraw, center = Ymu, scale = Ysd)

# dating of forecast horizons (on quarterly frequency) of interest
date_nc <- as.Date(obs_end) %m+% (months(3-run_mbyq))
date_bc <- seq(date_nc %m-% months(3), date_nc %m-% months(3), by = "quarter") # set here if more quarters are missing

if(fhorz > 0){
  date_fc <- seq(date_nc %m+% months(3), as.Date(max(datelabs_f)), by = "quarter")
  sl_fcdates_q <- data.frame("type" = paste0("hq",(-1):(fhorz/3)), 
                             "date" = c(date_bc, date_nc, date_fc))
  sl_fcdates_m <- data.frame("type" = paste0("hm",(-3):(fhorz + 3-run_mbyq )),
                             "date" = c(datelabs[(which(obs_end == datelabs) - 3):which(obs_end == datelabs)],
                                        as.character(seq(as.Date(obs_end) %m+% months(1), as.Date(max(datelabs_f)), by = "month"))))  
}

# --------------------------------------------------------------------------
# start of the main code
ntot <- nburn + nsave * nthin
mcmclabs <- paste0("mcmc",1:nsave)

set.mod.og <- run_mean
set.mean.tmp <- unlist(strsplit(set.mod.og,"-"))
set.mean <- set.mean.tmp[1]
if(length(set.mean.tmp) > 1){
  set.approx <- set.mean.tmp[2]  
}else{
  set.approx <- ""
}

id_spec <- paste0(obs_start,"_",set.mod.og,
                  "_",ifelse(outlier,"out","hom"))

dir_save <- paste0(dir_main,"/",id_spec)
if(file.exists(paste0(dir_save,".rda"))) stop("File already exists.")

# other settings
cons <- TRUE
shrink <- ifelse(set.approx == "hs", TRUE, FALSE)
n <- ncol(Yraw)

# subsample splits
datelabs_p <- datelabs[-c(1:p)]

# ------------------------------------------------------------------------------------------------
# fill up missings in raw matrix for initialization
Yraw_NA <- Yraw # original data with missings
for (nn in 1:n_lf) {
  Ytmp <- Yraw[, nn, drop = FALSE]
  init_obs <- min(which(!is.na(Ytmp)))
  Ytmp <- rep(Ytmp[seq(init_obs, nrow(Ytmp), by = 3)], each = 3) / 3
  Yraw[(nrow(Yraw) - NROW(Ytmp) + 1):nrow(Yraw), nn] <- Ytmp
}

# define matrix of lags
if (any(is.na(Yraw))) {
  for (i in 1:ncol(Yraw)) {
    yy <- Yraw[, i]
    if(any(is.na(yy))){
      yylags <- embed(yy[!is.na(yy)],2)
      y1 <- yylags[,1]
      x1 <- cbind(yylags[,-1],1)
      phic <- solve(crossprod(x1)) %*% crossprod(x1,y1)
      if(phic[1] >= 1) phic[1] <- 0.999
      
      yh <- y1[length(y1)]
      yfc_init <- rep(NA,sum(is.na(yy)))
      for(hh in 1:(sum(is.na(yy)))){
        yfc_init[hh] <- yh <- c(yh,1) %*% phic
      }
      
      Yraw[, i] <- c(yy[!is.na(yy)], yfc_init)
    }
  }
}

Y_init <- Yraw[1:p, ]
Ylags <- embed(Yraw, dim = p + 1)
rownames(Y_init) <- paste0("T", (-p + 1):0)
rownames(Ylags) <- paste0("T", 1:nrow(Ylags))
colnames(Ylags) <- paste0(rep(varlabs, p + 1), "_l", rep(0:p, each = n))
Y <- Ylags[, 1:n]

Y_full <- rbind(Y_init, Y)
colnames(Y) <- colnames(Yraw)

s2_ARp <- matrix(NA,n,1)
for(nn in 1:n){
  yy <- Yraw_NA[!is.na(Yraw_NA[,nn]),nn]
  sig2 <- try(arima(yy,order = c(p,0,0))$sigma2, silent = TRUE)
  s2_ARp[nn,1] <- sig2
}

# --------------------------------------------------------------------------
# match observed low-frequency data to lag structure
Ytmp <- Yraw_NA[(p + 1):nrow(Yraw_NA), 1:n_lf, drop = FALSE]
Ytmp_init <- Yraw_NA[1:p, 1:n_lf, drop = FALSE]
rownames(Ytmp) <- paste0("T", 1:nrow(Ytmp))
rownames(Ytmp_init) <- paste0("T", (-p + 1):0)

min_lf_obs <- min(which(apply(!is.na(Ytmp), 1, sum) == n_lf))
min_lf_init <- min(which(apply(!is.na(Ytmp_init), 1, sum) == n_lf))

Z <- Ytmp[seq(min_lf_obs, nrow(Ytmp), by = 3), , drop = FALSE]
Z_init <- Ytmp_init[seq(min_lf_init, nrow(Ytmp_init), by = 3), , drop = FALSE]
Z_full <- rbind(Z_init, Z)

ix_t_Y <- as.numeric(gsub("T", "", rownames(Y_full)))
ix_t_Z <- as.numeric(gsub("T", "", rownames(Z_full)))
ix_min_Z <- which(ix_t_Z - 5 == min(ix_t_Y))

Z_full <- Z_full[ix_min_Z:nrow(Z_full), , drop = FALSE]
Z_vec <- c(t(Z_full))
sl_na_lf <- is.na(Z_vec)
Z_vec <- Z_vec[!sl_na_lf]

T_obs_lf <- length(Z_vec)
X <- Ylags[, -c(1:n)]
X1l <- X[,1:n]

# holds other deterministic trends
if(cons){
  D <- matrix(1, nrow(X), 1)
  colnames(D) <- "cons"
  k_D <- ncol(D)
}else{
  k_D <- 0
  D <- matrix(0, nrow(X), k_D)
}

# --------------------------------------------------------------------------
# dimensions
k <- ncol(X)
T <- nrow(Y)

# some re-usable matrices
i_T <- matrix(1, T, 1)
i_n <- matrix(1, n, 1)
I_T <- Matrix(diag(T), sparse = TRUE)
I_n <- Matrix(diag(n), sparse = TRUE)
eps_small <- sqrt(.Machine$double.eps)

# some required indexing for missings
Yraw_NA_full <- Yraw_NA
Yraw_NA_full[, 1:n_lf] <- NA

id_lat <- is.na(Yraw_NA_full)
id_lat_vec <- c(t(id_lat))

Tn_lat <- sum(id_lat)
Tn_mea <- n * (T + p) - Tn_lat

Y_full_vec <- c(t(Y_full))
Yo_vec <- Y_full_vec[!id_lat_vec]

# setup for selection matrices
Slt_ls <- Smt_ls <- list()
for (tt in 1:(T + p)) {
  n_lat_t <- sum(id_lat[tt, ])
  n_mea_t <- n - n_lat_t

  Slt_tmp <- matrix(0, nrow = n, ncol = n_lat_t)
  Smt_tmp <- matrix(0, nrow = n, ncol = n_mea_t)
  lcount <- mcount <- 0
  for (nn in 1:n) {
    if (id_lat[tt, nn] == 1) {
      lcount <- lcount + 1
      Slt_tmp[nn, lcount] <- 1
    } else {
      mcount <- mcount + 1
      Smt_tmp[nn, mcount] <- 1
    }
  }
  Slt_ls[[tt]] <- Slt_tmp
  Smt_ls[[tt]] <- Smt_tmp
}

Sl <- bdiag(Slt_ls)
Sm <- bdiag(Smt_ls)

Sl_mat <- as.matrix(Sl)
Sm_mat <- as.matrix(Sm)

# labeling of positions
Ypos <- Y_full
Zpos <- Z_full
for (i in 1:n) {
  Ypos[, i] <- paste0(colnames(Y_full)[i], "_", rownames(Y_full))
}
for (i in 1:n_lf) {
  Zpos[, i] <- paste0(colnames(Z_full)[i], "_", rownames(Z_full))
}
Zpos[is.na(Z_full)] <- NA
Ypos_vec <- c(t(Ypos))
Zpos_vec <- c(t(Zpos))
Zpos_vec <- Zpos_vec[!sl_na_lf]

# loadings matrix for intertemporal restrictions
itr_Dload <- c(1 / 3, 2 / 3, 1, 2 / 3, 1 / 3) / 3
itr_Lload <- rep(1 / 3, 3)

Lambda <- matrix(0, T_obs_lf, Tn_lat)
colnames(Lambda) <- Ypos_subvec <- Ypos_vec[id_lat_vec]
rownames(Lambda) <- Zpos_vec

for (nn in 1:n_lf) {
  sl_var <- varlabs[nn]
  sl_Lambda_rows <- grepl(paste0(sl_var, "_"), Zpos_vec)

  Lambda_sub <- Lambda[sl_Lambda_rows, ]
  maxT <- as.numeric(gsub("T", "", do.call("rbind", strsplit(rownames(Lambda_sub), "_"))[, 2]))

  for (tt in seq_along(maxT)) {
    if (itr[nn] == "D") {
      Lambda_sub[tt, paste0(sl_var, "_T", maxT[tt]:(maxT[tt] - 4))] <- itr_Dload
    } else {
      Lambda_sub[tt, paste0(sl_var, "_T", maxT[tt]:(maxT[tt] - 2))] <- itr_Lload
    }
  }
  Lambda[sl_Lambda_rows, ] <- Lambda_sub
}

Lambda <- Matrix(Lambda, sparse = TRUE) # intertemporal restriction loadings
O_mat <- Matrix(1e-8 * diag(T_obs_lf), sparse = TRUE)
Oi_mat <- solve(O_mat)

# --------------------------------------------------------------------------
# priors and initialization
A_draw <- try(solve(crossprod(X)) %*% crossprod(X, Y), silent = TRUE)
if (is(A_draw, "try-error")) {
  A_draw <- solve(crossprod(X) + 1e-2 * diag(k)) %*% crossprod(X, Y)
}

a0_draw <- matrix(0,k_D,n)
FX <- fit <- X %*% A_draw

A_prior <- matrix(0, k, n)
A_prior[1:n, 1:n] <- diag(n) * 0
eps <- Y - fit

theta_a0 <- matrix(100, k_D, n) # prior on deterministic terms

# Global scalings (differentiate between own and other lags for Minnesota)
tau_A   <- rep(1, n)
zeta_A  <- rep(1, n)

# Local scalings accoriding to Minnesota moments
lam_A <- matrix(1, k, n)
nu_A  <- matrix(1, k, n)
psi_A <- theta_A  <- matrix(1, k, n)

theta_A <- t(tau_A * t(lam_A))

prior.sig <- c(10^50, 0.5)
cgm.level <- 0.95
cgm.exp <- 2
sd.mu <- 1.96
num.trees <- 250 #  number of trees

if (set.mean == "bart") {
  sig2.bart <- sig2.reg <- rep(1, n)
  control <- dbartsControl(
    verbose = FALSE, keepTrainingFits = TRUE, useQuantiles = FALSE,
    keepTrees = FALSE, n.samples = ntot,
    n.cuts = 100L, n.burn = nburn, n.trees = num.trees, n.chains = 1,
    n.threads = 1, n.thin = 1L, printEvery = 1,
    printCutoffs = 0L, rngKind = "default", rngNormalKind = "default",
    updateState = FALSE
  )
  sampler.run <- sampler.list <- list()
  for (jj in seq_len(n)) {
    sampler.list[[jj]] <- dbarts(Y[, jj] ~ X,
      control = control,
      tree.prior = cgm(cgm.exp, cgm.level), node.prior = normal(sd.mu),
      n.samples = nsave, weights = rep(1, T),
      sigma = 1, resid.prior = chisq(prior.sig[1], prior.sig[2])
    )
  }
}

# --------------------------------------------------------------------------
# prior on the covariance matrix
A_OLS <- try(solve(crossprod(X)) %*% crossprod(X,Y), silent = TRUE)
if(is(A_OLS, "try-error")){
  A_OLS <- try(solve(crossprod(X) + 0.01 * diag(ncol(X))) %*% crossprod(X,Y), silent = TRUE)  
}

Sigma_OLS <- crossprod(Y - X %*% A_OLS) / (T - n*p - cons)
sig2_OLS <- as.numeric(diag(Sigma_OLS))

# prior on the covariance matrix: avoid shrinkage of covariances
Sig_s0nu <- 2
sig2_tau <- 0.75

Sig_s0 <- Sig_s0nu + n - 1
Sig_A0 <- rep(1,n) # fixed scale parameter
for(i in 1:n){
  Sig_A0[i] <- stats::optim(par = 0.01, f = function(A){(sig2_tau - phalft(sig2_OLS[i], A, Sig_s0nu))^2},
                            method = "Brent",lower = 0.001,upper = 100)$par
}

Sig_A0[Sig_A0 < 0.1] <- 0.1
Sig_A0[Sig_A0 > 10] <- 10

Sig_a0 <- (1/Sig_A0^2) / (1/2 + 1)
Sig_S0 <- 2 * Sig_s0nu * diag(1/Sig_a0)

s1 <- Sig_s0 + T
Sigma_draw <- solve(matrix(rWishart(1, s1, solve(Sig_S0 + crossprod(Y - X %*% A_draw))), n, n))
Sigmai_draw <- solve(Sigma_draw)

# outlier detection setup
ot_draw <- rep(1, T)
ot_grid <- seq(1, 6, by = 1) # check here
ot_grid_ln <- length(ot_grid)
ot_p <- 1e-4
ot_p_grid <- c(1 - ot_p, rep(ot_p / (ot_grid_ln - 1), ot_grid_ln - 1))

ot_Ba <- 1
ot_Bb <- 50
Oi_T <- O_T <- I_T

# --------------------------------------------------------------------------
# set up the dynamic coefficients (here: standard normal)
H_init <- Matrix(cbind(kronecker(diag(p),diag(n)), matrix(0, n*p, T*n)), sparse = TRUE)
h_init <- matrix(0, n*p, 1)
Sig_init <- Matrix(kronecker(diag(p),diag(n)), sparse = TRUE)

Alag_sl <- matrix(seq(1, n * p), n)
H_ls <- list()
for (pp in 0:p) {
  if (pp == p) h10 <- 1 else h10 <- (-1)
  Htmp <- matrix(0, T, T + p)
  for (tt in 1:T) {
    Htmp[tt, tt + pp] <- h10
  }
  H_ls[[paste0("Hp", p - pp)]] <- Matrix(Htmp, sparse = TRUE)
}
H <- kronecker(H_ls[[paste0("Hp", 0)]], diag(n))
for (pp in 1:p) {
  H <- H + kronecker(H_ls[[paste0("Hp", pp)]], t(A_draw)[, Alag_sl[, pp]])
}
h <- kronecker(i_T, A_draw[k, ])

h <- rbind(h_init, h) # add prior on initial conditions
H <- rbind(H_init, H) # add prior on initial conditions
Sig_full <- bdiag(Sig_init, kronecker(I_T, solve(Sigma_draw)))

Gm <- H %*% Sm
Gl <- H %*% Sl

GltSigi <- crossprod(Gl, Sig_full)
K <- GltSigi %*% Gl
mu <- solve(K) %*% GltSigi %*% (h - Gm %*% Yo_vec)

M <- matrix(0,k,k)
M[(n+1):k,1:(n*(p-1))] <- diag(n*(p-1))

# approximation step shrinkage if required 
if(shrink){
  sig2_me <- rep(0.001, n)
  Sigma_me <- Matrix(diag(n) * sig2_me, sparse = TRUE)
}else{
  sig2_me <- rep(0, n)
  Sigma_me <- Matrix(diag(n) * 0, sparse = TRUE)
}

# tight prior on measurement errors
a0_me <- 10
b0_me <- 0.01

# Global scalings (differentiate between own and other lags for Minnesota)
tau_Aapx   <- 1
zeta_Aapx  <- 1
nu_Aapx  <- lam_Aapx <- matrix(1, k, n)
psi_Aapx <- matrix(1, k, n)
theta_Aapx <- t(tau_Aapx * t(lam_Aapx))

# --------------------------------------------------------------------------
# sampling
thin.set <- seq(nburn + 1, ntot, by = nthin)
savecount <- 0

A_store <- array(NA, dim = c(nsave, k, n))
Sig_store <- array(NA, dim = c(nsave, n, n))
ot_store <- array(NA, dim = c(nsave, T))

Y_store <- Ynl_store <- array(NA, dim = c(nsave, T + p + fhorz, n))

# naming 
dimnames(A_store) <- list("iter" = mcmclabs, "predictor" = colnames(X)[1:k], "variable" = varlabs)
dimnames(Sig_store) <- list("iter" = mcmclabs, "variable_r" = varlabs, "variable_c" = varlabs)
dimnames(Y_store) <- dimnames(Ynl_store) <- list("iter" = mcmclabs, "date" = c(datelabs,datelabs_f), "variable" = varlabs)
dimnames(ot_store) <- list("iter" = mcmclabs, "date" = datelabs_p)

# start the sampling loop
irep <- 1
pb <- txtProgressBar(min = 0, max = ntot, style = 3)
for (irep in 1:ntot) {
  # ----------------------------------------------------------------------
  # Part A: Sample the main VAR coefficients here (dynamic coefficients, covariances, etc.)
  # Step 1a: Sampling the VAR coefficients
  if (set.mean == "lin") {
    XD <- cbind(X,D)
    for(nn in 1:n){
      Sn0 <- Sigma_draw[nn,-nn,drop=FALSE] %*% solve(Sigma_draw[-nn,-nn,drop=FALSE])
      m_nn <- t(tcrossprod(Sn0,(Y - XD %*% rbind(A_draw,a0_draw))[,-nn]))
      S_nn <- as.numeric(Sigma_draw[nn,nn,drop=FALSE] - Sn0 %*% Sigma_draw[-nn,nn,drop=FALSE])
      
      YY <- (Y[,nn] - m_nn) / (sqrt(S_nn) * ot_draw)
      XX <- XD / (sqrt(S_nn) * ot_draw)
      
      VA_inv <- diag(k + k_D) / c(theta_A[,nn],theta_a0[,nn])
      VA_post <- solve(crossprod(XX) + VA_inv)
      A_post <- VA_post %*% (VA_inv %*% c(A_prior[,nn], rep(0,k_D)) + crossprod(XX,YY))
      
      Aa0_draw <- as.numeric(A_post + t(chol(VA_post)) %*% rnorm(k + k_D))
      A_draw[,nn] <- Aa0_draw[1:k]
      a0_draw[,nn] <- Aa0_draw[-c(1:k)]
    }
    
    fit <- cbind(X,D) %*% rbind(A_draw,a0_draw)
    eps <- Y - fit
  } else if (set.mean == "bart") {
    Xginv <- ginv(X)
    for (nn in 1:n){
      Sn0 <- Sigma_draw[nn,-nn,drop=FALSE] %*% solve(Sigma_draw[-nn,-nn,drop=FALSE])
      m_nn <- t(tcrossprod(Sn0,(Y - FX - D %*% a0_draw)[,-nn]))
      S_nn <- as.numeric(Sigma_draw[nn,nn,drop=FALSE] - Sn0 %*% Sigma_draw[-nn,nn,drop=FALSE])
      
      # BART part
      sampler.list[[nn]]$setResponse((Y[,nn] - (D %*% a0_draw)[,nn] - m_nn))
      sampler.list[[nn]]$setPredictor(X)
      sampler.list[[nn]]$setWeights(1 / (ot_draw^2 * S_nn))
      
      rep_nn <- sampler.list[[nn]]$run(0L, 1L) # run the BART update step
      sampler.run[[nn]] <- rep_nn
      sig2.bart[nn] <- rep_nn$sigma
      
      # FX[, nn] <- rep_nn$train
      FX[, nn] <- rep_nn$train - mean(rep_nn$train)
      
      # variables that enter linearly
      m0_nn <- t(tcrossprod(Sn0,(Y - FX - D %*% a0_draw)[,-nn]))
      YY <- (Y[,nn] - FX[, nn] - m0_nn) / (sqrt(S_nn) * ot_draw)
      XX <- D / (sqrt(S_nn) * ot_draw)
      
      if(ncol(XX) > 0){
        VA_inv <- diag(k_D) / c(theta_a0[,nn])
        VA_post <- solve(crossprod(XX) + VA_inv)
        A_post <- VA_post %*% (VA_inv %*% c(rep(0, k_D)) + crossprod(XX,YY))
        
        a0_tmp <- as.numeric(A_post + t(chol(VA_post)) %*% rnorm(k_D))
        a0_draw[,nn] <- a0_tmp[1:k_D]
      }
    }
    
    fit <- FX + D %*% a0_draw
    eps <- Y - fit
  }
  
  # approximate the conditional mean linearly
  XpX <- crossprod(X)
  Xpfi <- crossprod(X,FX)
  
  if(set.mean == "bart"){
    if(shrink){
      # F = XA + u, approximate shrinkage solution for A
      for(i in 1:n){
        A_po <- solve(XpX / sig2_me[i] + diag(k) / theta_Aapx[,i])
        a_po <- A_po %*% (Xpfi[,i] / sig2_me[i])
        Aapx_draw <- a_po + t(chol(A_po)) %*% rnorm(k)
        
        A_draw[,i] <- Aapx_draw
        sig2_me[i] <- 1/rgamma(1,a0_me + T/2, b0_me + crossprod(FX[,i] - X %*% Aapx_draw) / 2)
      }
      Sigma_me <- diag(n) * sig2_me
      
      # shrinkage
      hs_Aapx <- get.hs.min(bdraw     = c(A_draw),
                            lam.hs    = c(lam_Aapx), 
                            nu.hs     = c(nu_Aapx), 
                            tau.hs    = tau_Aapx, 
                            zeta.hs   = zeta_Aapx, 
                            update.ls = TRUE)
      psi_Aapx[] <- hs_Aapx$psi
      lam_Aapx[] <- hs_Aapx$lam; nu_Aapx[]  <- hs_Aapx$nu
      tau_Aapx[]  <- hs_Aapx$tau; zeta_Aapx[] <- hs_Aapx$zeta
      
      psi_Aapx[psi_Aapx < 1e-10] <- 1e-10
      psi_Aapx[psi_Aapx > 10] <- 10
      theta_Aapx <- psi_Aapx
    }else{
      if(k > (0.8 * T)){
        Xginv <- solve(crossprod(X) + 0.001*diag(k)) %*% t(X) # add a bit of shrinkage in case there's too little observations
      }else{
        Xginv <- ginv(X) # if rank(X) = k, then same as solve(XpX) %*% t(X), due to F \approx XA -> (X'X)^-1 X'F \approx A
      }
      
      A_draw <- Xginv %*% FX
      Sigma_me <- diag(n) * 0 
    }
  }
  
  # Step 1b: Sampling variances
  Sig_a0 <- 1 / rgamma(n, 0.5 * (Sig_s0nu + T), Sig_s0nu * diag(Sigmai_draw) + (1 / Sig_A0^2))
  Sig_S0 <- 2 * Sig_s0nu * diag(1/Sig_a0) # update the prior scaling matrix
  
  s1 <- Sig_s0 + T
  S1 <- Sig_S0 + crossprod(eps / ot_draw)
  
  Sigmai_draw <- matrix(rWishart(1, s1, solve(S1)), n, n)
  Sigma_draw <- solve(Sigmai_draw)

  if (outlier) {
    loglikmat <- matrix(0, T, ot_grid_ln)
    for (i in 1:ot_grid_ln) {
      loglikmat[, i] <- mvnfast::dmvn(eps, mu = rep(0,n), sigma = ot_grid[i]^2 * Sigma_draw, log = TRUE)
    }
    probs <- t(ot_p_grid * t(exp(loglikmat)))
    probs <- probs / apply(probs, 1, sum)
    
    # sample the outlier scaling
    for (tt in 1:T) {
      ot_draw[tt] <- sample(ot_grid, size = 1, prob = probs[tt, ])
    }
    
    ot_sum <- sum(ot_draw != 1)
    ot_p <- rbeta(1, ot_Ba + ot_sum, ot_Bb + T - ot_sum)
    ot_p_grid <- c(1 - ot_p, rep(ot_p / (ot_grid_ln - 1), ot_grid_ln - 1))
    
    diag(O_T) <- (ot_draw^2)
    diag(Oi_T) <- 1 / (ot_draw^2)
  }
  
  # Step 1c: shrinkage
  if (set.mean == "lin") {
    hs_A <- get.hs.min(bdraw     = c(A_draw - A_prior),
                       lam.hs    = c(lam_A), 
                       nu.hs     = c(nu_A), 
                       tau.hs    = tau_A[1], 
                       zeta.hs   = zeta_A[1], 
                       update.ls = TRUE)
    psi_A[] <- hs_A$psi
    lam_A[] <- hs_A$lam; nu_A[]  <- hs_A$nu
    tau_A[]  <- hs_A$tau; zeta_A[] <- hs_A$zeta
    
    psi_A[psi_A < 1e-8] <- 1e-8
    psi_A[psi_A > 10] <- 10
    theta_A <- psi_A
  }

  # ----------------------------------------------------------------------
  # Part B: sample the latent mixed-frequency states conditional on the VAR parameters
  h <- matrix(c(t(D %*% a0_draw)), ncol = 1) # set up companion full data matrices
  H <- kronecker(H_ls[[paste0("Hp", 0)]], diag(n))
  for (pp in 1:p) {
    H <- H + kronecker(H_ls[[paste0("Hp", pp)]], t(A_draw)[, Alag_sl[, pp]])
  }
  h <- rbind(h_init, h) # add prior on initial conditions
  H <- rbind(H_init, H) # add prior on initial conditions
  Sig_T <- solve(kronecker(O_T, Sigma_draw) + kronecker(I_T, Sigma_me))
  Sig_full[-c(1:(n*p)),-c(1:(n*p))] <- Sig_T
  
  Gm <- H %*% Sm
  Gl <- H %*% Sl

  GltSigi <- crossprod(Gl, Sig_full)
  Sigbar_i <- GltSigi %*% Gl # precision

  # unconditional of observed low-frequency info ---
  Sigbar_i_chol <- chol(Sigbar_i)
  mubar <- solve(Sigbar_i_chol, solve(t(Sigbar_i_chol), GltSigi %*% (h - Gm %*% Yo_vec))) # precision-based

  LambdaOi <- crossprod(Lambda, Oi_mat)
  Sigbar_const_i <- LambdaOi %*% Lambda + Sigbar_i # precision

  Sigbar_const_i_chol <- try(chol(Sigbar_const_i), silent = TRUE)
  if (is(Sigbar_const_i_chol, "try-error")) {
    diag(Sigbar_const_i) <- diag(Sigbar_const_i) + eps_small # in case of numerical issues, add a small offsetting constant
    Sigbar_const_i_chol <- chol(Sigbar_const_i)
  }

  # conditioning on low-frequency info ---
  mubar_const <- solve(Sigbar_const_i_chol, solve(t(Sigbar_const_i_chol), LambdaOi %*% Z_vec + Sigbar_i %*% mubar)) # precision-based
  Ym_vec <- mubar_const + solve(Sigbar_const_i_chol, rnorm(Tn_lat))

  # reconfigure data for the next iteration
  Y_full_vec <- Sm %*% Yo_vec + Sl %*% Ym_vec
  Y_full_raw <- matrix(Y_full_vec, T + p, n, byrow = TRUE)
  Y_full <- scale(Y_full_raw) # normalize so that each dataset has again mean 0 and sd 1

  Ylags <- embed(Y_full, dim = p + 1)
  Y <- Ylags[, 1:n]
  X <- Ylags[, -c(1:n)]

  # ----------------------------------------------------------------------
  # storage
  if (irep %in% thin.set) {
    savecount <- savecount + 1
    A_store[savecount, , ] <- A_draw[1:k,]
    Sig_store[savecount, , ] <- Sigma_draw
    ot_store[savecount, ] <- ot_draw
    Sig_chol <- t(chol(Sigma_draw))
    
    # nonlinear fit
    Ynl <- matrix(NA, T, n)
    if(set.mean == "bart"){
      for(nn in 1:n){
        Ynl[, nn] <- sampler.list[[nn]]$predict(X)
      }
      if(cons){
        Ynl <- Ynl + (matrix(1, T, 1) %*% a0_draw)
      }
      
      # sample including variation from the error term
      for(tt in 1:T){
        Ynl[tt,] <- Ynl[tt,,drop=FALSE] + t(ot_draw[tt] * Sig_chol %*% rnorm(n))
      }  
    }
    
    # compute forecasts
    if(fhorz > 0){
      Xf <- c(Y[T, ], X[T, 1:(n * (p-1))])
      Yf <- Yfmu <- matrix(0, fhorz, n)
      
      for(hh in 1:fhorz){
        # sample a shock
        if(outlier){
          ot_f <- sample(c(1, sample(ot_grid[-1], size = 1)), size = 1, replace = TRUE, prob = c(1 - ot_p, ot_p))
        }else{
          ot_f <- 1
        }
        u_th <- ot_f * Sig_chol %*% rnorm(n)
        
        if(set.mean == "bart"){
          for (nn in 1:n) {
            Yfmu[hh,nn] <- sampler.list[[nn]]$predict(Xf)
          }
        }else{
          Yfmu[hh, ] <- Xf %*% A_draw
        }
        if(cons){
          Yfmu[hh, ] <- Yfmu[hh, ] + (matrix(1, 1, 1) %*% a0_draw)
        }
        
        Yf[hh, ] <- Yfmu[hh, ] + u_th
        Xf <- c(Yf[hh, ], Xf[1:(n * (p-1))])
      }
    }
    
    # map to original scale
    Yq_rep <- Yqnl_rep <- matrix(NA, T + p + fhorz, n)
    
    for(j in 1:n){
      if(j <= n_lf){
        if (itr[j] == "D") {
          Yq_rep[, j] <- stats::filter(c(Y_full_raw[,j],Yf[,j]), 
                                       filter = itr_Dload, method = "convolution", sides = 1) * Ysd[j] + Ymu[j]
          Yqnl_rep[, j] <- stats::filter(c(Y_full_raw[1:p,j],Ynl[,j],Yf[,j]), 
                                         filter = itr_Dload, method = "convolution", sides = 1) * Ysd[j] + Ymu[j]
        } else {
          Yq_rep[, j] <- stats::filter(c(Y_full_raw[,j],Yf[,j]), 
                                       filter = itr_Lload, method = "convolution", sides = 1) * Ysd[j] + Ymu[j]
          Yqnl_rep[, j] <- stats::filter(c(Y_full_raw[1:p,j],Ynl[,j],Yf[,j]), 
                                         filter = itr_Lload, method = "convolution", sides = 1) * Ysd[j] + Ymu[j]
        } 
      }else{
        Yq_rep[, j] <- c(Y_full_raw[,j],Yf[,j]) * Ysd[j] + Ymu[j]
        Yqnl_rep[, j] <- c(Y_full_raw[1:p,j],Ynl[,j],Yf[,j]) * Ysd[j] + Ymu[j]
      }
    }
    
    Y_store[savecount, , ] <- Yq_rep
    if(set.mean == "bart"){
      Ynl_store[savecount, , ] <- Yqnl_rep
    }else{
      Ynl_store[savecount, , ] <- Yq_rep
    }
  }
  setTxtProgressBar(pb, irep)
}

Yq_store <- Y_store
save(Yq_store, file = paste0(dir_main,"/",run_id,".rda")) # save the mcmc draws (can take large amounts of storage depending on settings)

# some outputs
Yq_quant <- apply(apply(Yq_store,c(2,3),remove_outliers),c(2,3),quantile,probs=grid_tau,na.rm=TRUE)

plot_save <- paste0(dir_plots,"/",run_id)
dir.create(plot_save, showWarnings = FALSE)
for(sl_var in output_var){
  # plot of time series
  pdf(file = paste0(plot_save,"/ts_",sl_var,".pdf"), width = 10, height = 4)
  par(mfrow = c(1, 2), mar = c(4,4,1,1))
  ts.plot(ts(cbind(0, t(Yq_quant[c("5%","95%"),,sl_var])), start = obs_start_2, frequency = 12), ylab = sl_var)
  points(Yraw_orig[,sl_var])
  ts.plot(window(ts(cbind(0, t(Yq_quant[c("5%","95%"),,sl_var])), start = obs_start_2, frequency = 12), start = c(2022,1)), ylab = "zoomed")
  dev.off()
  
  # plot of distribution
  pdf(file = paste0(plot_save,"/dist_",sl_var,".pdf"), width = 10, height = 4)
  par(mfrow = c(1, 3), mar = c(4,4,1,1))
  plot(density(Yq_store[,as.character(date_bc),sl_var]), lwd = 2, 
       xlab = paste0("Backcast (",as.yearqtr(date_bc),")"), main = "")
  abline(v = c(0, Yq_quant[c("5%","50%","95%"),as.character(date_bc),sl_var]), lwd = c(1,1,2,1), col = c("red",rep("black",3)))
  med_num <- Yq_quant["50%",as.character(date_bc),sl_var]
  text(x = med_num, y = 0, labels = sprintf("%.2f", med_num), pos = 4, offset = 0.5, cex = 1)
  
  plot(density(Yq_store[,as.character(date_nc),sl_var]), lwd = 2, 
       xlab = paste0("Nowcast (",as.yearqtr(date_nc),")"), main = sl_var)
  abline(v = c(0, Yq_quant[c("5%","50%","95%"),as.character(date_nc),sl_var]), lwd = c(1,1,2,1), col = c("red",rep("black",3)))
  med_num <- Yq_quant["50%",as.character(date_nc),sl_var]
  text(x = med_num, y = 0, labels = sprintf("%.2f", med_num), pos = 4, offset = 0.5, cex = 1)
  
  plot(density(Yq_store[,as.character(date_fc[1]),sl_var]), lwd = 2, 
       xlab = paste0("Forecast (",as.yearqtr(date_fc[1]),")"), main = "")
  abline(v = c(0, Yq_quant[c("5%","50%","95%"),as.character(date_fc[1]),sl_var]), lwd = c(1,1,2,1), col = c("red",rep("black",3)))
  med_num <- Yq_quant["50%",as.character(date_fc[1]),sl_var]
  text(x = med_num, y = 0, labels = sprintf("%.2f", med_num), pos = 4, offset = 0.5, cex = 1)
  dev.off()
}






