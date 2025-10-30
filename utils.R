phalft <- function(x, A, nu){
  return(2 * stats::pt(x/A, nu) - 1)
}

get.hs.min <- function(bdraw,lam.hs,nu.hs,tau.hs,zeta.hs,update.ls = TRUE){
  k <- length(bdraw)
  # Local shrinkage scalings
  if(update.ls){
    lam.hs <- 1/rgamma(k,shape=1,rate=1/nu.hs+bdraw^2/(2*tau.hs))
    nu.hs <-  1/rgamma(k,shape=1,rate=1+1/lam.hs)
  }else{
    lam.hs <- lam.hs
    nu.hs <- nu.hs
  }
  # Global shrinkage parameter
  tau.hs  <- 1/rgamma(1,shape=(k+1)/2,rate=1/zeta.hs+sum(bdraw^2/lam.hs)/2) 
  zeta.hs <- 1/rgamma(1,shape=1,rate=1+1/tau.hs)
  
  ret <- list("psi"=(lam.hs*tau.hs),"lam"=lam.hs,"tau"=tau.hs,"nu"=nu.hs,"zeta"=zeta.hs)
  return(ret)
}

as.Date.mbyq <- function(x){
  x_date <- as.Date(x)
  x_m <- format(x_date, "%m")
  check <- formatC(matrix(seq(1,12),3,4), width = 2, flag = "0")
  return(which(check == x_m, arr.ind = TRUE)[1])
}

# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs = c(.25, .75), na.rm = na.rm, ...)
  H <- 2 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  return(y)
}
