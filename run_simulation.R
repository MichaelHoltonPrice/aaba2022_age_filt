library(demohaz)

rm(list=ls())
source("simulation_functions.R")

# Age specific mortality
# Siler parameter vector from Gage and Dyke 1986, Table 2, Level 15.
# Gage and Dyke 1986 -- Parameterizing abridged mortality tables: the Siler
#                       three-component hazard model
a <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)

dx <- 0.01
xmax <- 120
N <- 100

th <- c(2e-2, 1.2, a)

set.seed(1000)
# A simulation with no age filtration
sim1 <- sample_usher3(N, th, dx, xmax)
# A simulation with age filtration
sim2 <- sample_usher3(N, th, dx, xmax, x_mid=50, infant_prop=0.1)

# A wrapper function with a transformation to ensure that parameters are
# positive
nll_usher3_wrapper <- function(th_bar, x, ill) {
  th <- exp(th_bar)
  return(nll_usher3(th, x, ill))
}

# Fit the two simulations
th0 <- th
th_bar0 <- log(th0)
fit1 <- demohaz::temper_and_tune(nll_usher3_wrapper, th_bar0, x=sim1$x, ill=sim1$ill)
fit2 <- demohaz::temper_and_tune(nll_usher3_wrapper, th_bar0, x=sim2$x, ill=sim2$ill)
th1 <- exp(fit1$th)
th2 <- exp(fit2$th)

xcalc <- sim1$xcalc
xcalc <- xcalc[xcalc <= 80]
haz0 <- hsiler(xcalc, a)
haz1 <- hsiler(xcalc, th1[3:7])
haz2 <- hsiler(xcalc, th2[3:7])

haz_max <- max(haz0, haz1, haz2)
pdf("hazard.pdf")
  plot(xcalc, haz0, ylim=c(0,haz_max), type="l", col="blue", lwd=3)
  lines(xcalc, haz1, col="black", lwd=3)
  lines(xcalc, haz2, col="red"  , lwd=3)
  legend("topright",
         c("Sim", "Unfiltered", "Filtered"),
         col=c("blue", "black", "red"),
         lwd=3)
dev.off()

xcalc <- seq(0,100, by=dx)
rho1_0 <- usher3_rho1(xcalc, th0[1], a)
rho2_0 <- usher3_rho2(xcalc, th0[1], th0[2], a)
rho1_1 <- usher3_rho1(xcalc, th1[1], th1[3:7])
rho2_1 <- usher3_rho2(xcalc, th1[1], th[2], th1[3:7])
rho1_2 <- usher3_rho1(xcalc, th2[1], th2[3:7])
rho2_2 <- usher3_rho2(xcalc, th2[1], th2[2], th2[3:7])

rho_max <- max(rho1_0, rho2_0, rho1_1, rho2_1, rho1_2, rho2_2)
pdf("rho.pdf")
  plot(xcalc,  rho1_0, ylim=c(0,rho_max), type="l", col="blue", lwd=3,
       xlab="Age [years]", ylab="Age Distribution")
  lines(xcalc, rho2_0, col="blue", lwd=3, lty=3)
  lines(xcalc, rho1_1, col="black", lwd=3)
  lines(xcalc, rho2_1, col="black", lwd=3, lty=3)
  lines(xcalc, rho1_2, col="red", lwd=3)
  lines(xcalc, rho2_2, col="red", lwd=3, lty=3)
  legend("topright",
         c("Sim", "Unfiltered", "Filtered"),
         col=c("blue", "black", "red"),
         lwd=3)
dev.off()

# Simulate 100 cases with and without filtration, saving the fitted parameter
# vectors
num_sim <- 100
TH_uflt <- matrix(NA,7,num_sim)
TH_filt <- matrix(NA,7,num_sim)
print("Starting simulation loop")
for (n_sim in 1:num_sim) {
  print(n_sim)
  # A simulation with no age filtration
  sim1 <- sample_usher3(N, th, dx, xmax)
  fit1 <- demohaz::temper_and_tune(nll_usher3_wrapper,
                                   th_bar0,
                                   x=sim1$x,
                                   ill=sim1$ill)
  th1 <- exp(fit1$th)
  TH_uflt[,n_sim] <- th1
  # A simulation with age filtration
  sim2 <- sample_usher3(N, th, dx, xmax, x_mid=50, infant_prop=0.1)
  fit2 <- demohaz::temper_and_tune(nll_usher3_wrapper,
                                   th_bar0,
                                   x=sim2$x,
                                   ill=sim2$ill)
  th2 <- exp(fit2$th)
  TH_filt[,n_sim] <- th2

  ind_done <- which(!is.na(TH_uflt[2,]))
  print("--")
  print(th0[2])
  print(mean(TH_uflt[2,ind_done]))
  print(sd(TH_uflt[2,ind_done]))
  print(mean(TH_filt[2,ind_done]))
  print(sd(TH_filt[2,ind_done]))
}

# Relative errors
re_k1_uflt <- (TH_uflt[1,] - th0[1])/th0[1]
re_k1_filt <- (TH_filt[1,] - th0[1])/th0[1]
re_k2_uflt <- (TH_uflt[2,] - th0[2])/th0[2]
re_k2_filt <- (TH_filt[2,] - th0[2])/th0[2]

hist_uflt_k2 <- hist(re_k2_uflt, plot = FALSE)
hist_filt_k2 <- hist(re_k2_filt, plot = FALSE)

xrange <- range(hist_uflt_k2$breaks, hist_filt_k2$breaks)
pdf("k2_hist.pdf")
plot(hist_uflt_k2,
     col = rgb(1,0,0,0.4),
     xlab = 'k2 value',
     ylab="Density",
     main=NULL,
     freq=FALSE,
     xlim=xrange)
plot(hist_filt_k2,
     xaxt = 'n',
     yaxt = 'n',
     col = rgb(0,0,1,0.4),
     add = TRUE,
     freq = FALSE)
dev.off()