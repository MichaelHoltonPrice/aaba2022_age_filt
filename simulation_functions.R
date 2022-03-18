# Calculate the negative log-likelihood (nll) for the three-state Usher
# illness-death model.
#
# theta		The parameter vector for the Usher illness-death model with
#		ordering [k1,k2,mort]. The mortality parameter vector, mort,
#		is [a1,b1,a2,a3,b3] if juv is TRUE (Siler mortality) or
#		[a2,a3,b3] if juv is FALSE (Gompertz-Makeham mortality.
#
# x		The vector of ages-at-death
# ill		The vector of illness indicates, which can have NA entries
# x0		Conditional starting age [default 0]
# juv		Whether to include a juvenile mortality term [default TRUE]
nll_usher3 <- function(theta,x,ill,x0=0) {

  k1 <- theta[1]
  if(k1 < 0) {
    return(Inf)
  }

  k2 <- theta[2]
  if(k2 < 0) {
    return(Inf)
  }

  if(any(theta[3:7] < 0)) {
    return(Inf)
  }
  a_siler <- theta[3:7]

  # No missing allowed in ill
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]

  rho2_ill <- usher3_rho2(x_ill,k1,k2,a_siler,x0)
  if(any(is.na(rho2_ill))) {
    return(Inf)
  }

  rho1_wll <- usher3_rho1(x_wll,k1,a_siler,x0)

  ll <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  return(-ll)
}

usher3_rho1 <- function(x,k1,a_siler,x0=0) {
  f13 <- dsiler(x,a_siler,x0)
  S12 <- exp(-k1*(x - x0))
  rho1 <- f13 * S12
  return(rho1)
}

usher3_rho2 <- function(x,k1,k2,a_siler,x0=0) {
  f13_0_x <- dsiler(x,a_siler)
  S13_0_x <- ssiler(x,a_siler)
  S12_0_x0 <- exp(-k1*x0)
  S13_0_x0 <- ssiler(x0,a_siler)

  integralTerm <- rep(NA,length(x))
  for(ii in 1:length(x)) {
    integralTerm[ii] <- tryCatch(integrate(usher3_integrand,
                                           x0,x[ii],
                                           k1=k1,k2=k2,
                                           a_siler=a_siler)$value,
                                 error=function(e){NA})
  }

  rho2 <- k1*k2 * f13_0_x * S13_0_x^(k2-1) * integralTerm / S12_0_x0 / S13_0_x0
  return(rho2)
}

usher3_integrand <- function(y,k1,k2,a_siler) {
  # _0_y  indicates that the relevant quantity is relative to the interval 0 to y
  S12_0_y <- exp(-k1*y)
  S13_0_y <- ssiler(y,a_siler)

  return(S12_0_y * (S13_0_y)^(1-k2))
}

usher3_hessian <- function(theta,x,ill,x0=0,juv=T) {
  H <- numDeriv::hessian(usher3_ll,theta,method.args=list(eps=1e-12),ageVect=x,illVect=ill,x0=x0,juv=juv)
  return(H)
}

# In part, use this because numDeriv::hessian has the named variable x
usher3_ll <- function(paramVect,ageVect,illVect,x0=0,juv=T) {
  return(nll_usher3(paramVect,ageVect,illVect,x0,juv))
}

usher3_errors <- function(theta,x,ill,x0=0,juv=T) {
  H <- usher3_hessian(theta,x,ill,x0,juv)
  if(juv) {
    against <- c(0,1,0,0,0,0,0)
    sideAdjustment <- c(1,2,1,1,1,1,1)
    varName <- c('k1','k2','a1','b1','a2','a3','b3')
  } else {
    against <- c(0,1,0,0,0)
    sideAdjustment <- c(1,2,1,1,1)
    varName <- c('k1','k2','a2','a3','b3')
  }
  seVect <- sqrt(diag(solve(H)))
  zVect <- (theta-against) / seVect
  pvalVect <- sideAdjustment * pnorm(-abs(zVect))
  outputDf <- data.frame(Estimate=theta,StandErr=seVect,z=zVect,pval=pvalVect,against=against,sideAdj=sideAdjustment)
  rownames(outputDf) <- varName
  return(outputDf)
}

calc_filtation_density <- function(xcalc, x_mid, infant_prop, discrete=T) {
  z_mid <- 1 / x_mid / (1+infant_prop)
  z_inf <- infant_prop * z_mid

  f <- rep(0,length(xcalc))
  ind_lo <- which(xcalc <= x_mid)
  f[ind_lo] <- z_inf + (z_mid - z_inf) * xcalc[ind_lo] / x_mid

  ind_hi <- which(xcalc > x_mid & (xcalc <= 2*x_mid))
  f[ind_hi] <- z_mid + (z_inf - z_mid)*(xcalc[ind_hi] - x_mid) / x_mid

  #ind_too_hi <- which(xcalc > 2*x_mid)
  if (discrete) {
    f <- f / sum(f)
  }
  return(f)
}


# Although the simulated samples could in principle be sampled from in an
# "exact" manner using the underlying densities (rho1 and rho2), it is far
# easier to discretize the distributions and sample with a fine age spacing and
# sample from the discrete probabilities.
sample_usher3 <- function(N, th, dx, xmax, x_mid=NA, infant_prop=NA) {
  filter_by_age <- !is.na(x_mid)
  k1 <- th[1]
  k2 <- th[2]
  a_siler <- th[3:7]

  xcalc <- seq(0, xmax, by=dx)
  if (filter_by_age) {
    f_filt <- calc_filtation_density(xcalc, x_mid, infant_prop, discrete=T)
  }
  rho1 <- usher3_rho1(xcalc, k1, a_siler)
  rho2 <- usher3_rho2(xcalc,k1,k2,a_siler)
  area1 <- dx*sum(rho1)
  area2 <- dx*sum(rho2)

  total_area <- area1 + area2
  if (abs(total_area - 1) >= .001) {
    stop("total_area is not 1 to acceptable tolerance")
  }

  # The probability of having an ill observation
  prob_ill <- area2 / total_area

  # The density for well observations (normalized to 1)
  f_well <- rho1 / sum(rho1)
  if (filter_by_age) {
    f_well <- f_well * f_filt
    f_well <- f_well / sum(f_well)
  }

  # The density for ill observations (normalized to 1)
  f_ill <- rho2 / sum(rho2)
  if (filter_by_age) {
    f_ill <- f_ill * f_filt
    f_ill <- f_ill / sum(f_ill)
  }

  x <- rep(NA, N)
  ill <- rep(NA, N)

  for (n in 1:N) {
    ill[n] <- runif(1) <= prob_ill
    if (ill[n]) {
      # Sample from f_ill
      i <- sample(1:length(xcalc), 1, replace=F, prob=f_ill)
    } else {
      # Sample from f_well
      i <- sample(1:length(xcalc), 1, replace=F, prob=f_well)
    }
    x[n] <- xcalc[i]
  }

  return(list(xcalc=xcalc, rho1=rho1, rho2=rho2, x=x, ill=ill))
}


#sample_usher3_extension <- function(numSamp,param,missProp=0) {
#  k0 <- param[1]
#  k1 <- param[2]
#  k2 <- param[3]
#
#  a_siler <- param[4:8]
#
#  # Transitions out of the well state (1)
#  x_dead <- rsiler(numSamp,a_siler)
#
#  if(k0 <= 0) {
#    x_ill <- rsiler(numSamp,c(k1,k0,0,0,1))
#  } else {
#    x_ill <- rsiler(numSamp,c(0,1,0,k1,k0))
#  }
#
#  x <- rep(NA,numSamp)
#  ill <- rep(T,numSamp)
#  bin <- x_ill >= x_dead
#  x[bin] <- x_dead[bin]
#  ill[bin] <- F
#
#  x_dead_from_ill <- rep(NA,sum(bin))
#  ind <- which(bin)
#  for(ii in 1:length(ind)) {
#    x_dead_from_ill[ii] <- rsiler(1,k2*a_siler,x0=x_ill[ind[ii]])
#  }
#
#  x[bin] <- x_dead_from_ill
#  return(list(x=x,ill=ill))
#}

#integrate_rho1_usher3 <- function(param) {
#  k1 <- param[1]
#  a_siler <- param[3:7]
#  output <- tryCatch(integrate(usher3_rho1,0,120,k1=k1,a_siler=a_siler,x0=0)$value,error=function(e){NA})
#  return(output)
#}
#
#integrate_rho1_leh <- function(param) {
#  output <- tryCatch(integrate(leh_rho1,0,120,th=param)$value,error=function(e){NA})
#  return(output)
#}