
library(nimble)
library(nimbleJSextras)
library(here)

local_wd <- file.path(here(), "inst", "examples", "simulated_RPT")
setwd(local_wd)

#' -----------------------------------------------------------------------------
#' Create simulated finite mixture JS data
#' -----------------------------------------------------------------------------

set.seed(8675309)
# Number of occasions
K <- 36

# Mixture probabilities (1="resident", 2="part-time resident", 3="transient")
alpha <- c(.2, .25, .55)
# Entry probabilities
xi <- cbind(
  c(0.5,0.5*rep(1/(K-1), K-1)),
  c(0.5,0.5*rep(1/(K-1), K-1)),
  c(seq(1,10,length=18),seq(10,1,length=18))
  )
xi[,3] <- xi[,3]/sum(xi[,3])
# Exit probabilities
phi_RP <- 0.75
phi_T <- 0.01
phi_tilde <- cbind(
  rep(phi_RP, K-1),
  rep(phi_RP, K-1),
  rep(phi_T, K-1)
)^(1/12)
# capture probabilites
pRT <- rbeta(K,10,57)
delta_p <- 0.7
pP <- delta_p * pRT
p <- cbind(pRT, pP, pRT)

#' -----------------------------------------------------------------------------
#' Capture data simulation
#' -----------------------------------------------------------------------------
sim <- simulate_mix_ch(1000, alpha, xi, phi_tilde, p)

plot(colSums(sim$states==2))


save(list=ls(), file="rpt_data.RData")

