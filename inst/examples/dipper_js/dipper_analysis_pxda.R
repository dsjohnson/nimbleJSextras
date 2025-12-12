library(nimble)
library(coda)
library(marked)
library(tidyverse)

#' -----------------------------------------------------------------------------
#' Load the dipper data
#' -----------------------------------------------------------------------------
data("dipper")
DH <- strsplit(dipper$ch, "") %>% lapply(.,as.numeric) %>% do.call(rbind,.)
y.obs <- as.matrix(DH)
y.obs[is.na(y.obs)] <- 0   # ensure no NAs

n.obs <- nrow(y.obs)
K <- ncol(y.obs)

#' -----------------------------------------------------------------------------
#' Data augmentation
#' -----------------------------------------------------------------------------
n.aug <- 300  # choose augmented size (tune as needed)
M <- n.obs + n.aug
y_aug <- matrix(0, nrow = M, ncol = K)
y_aug[1:n.obs, ] <- y.obs

constants <- list(M = M, K = K)
data <- list(y = y_aug)

#'------------------------------------------------------------------------------
#' PX-DA code
#' -----------------------------------------------------------------------------
pxda_code <- nimbleCode({
  for(t in 1:K){
    logit_p[t] ~ dnorm(mu_p, sd=sig_p)
    p[t] <- expit(logit_p[t])
  }
  mu_p ~ dnorm(0,sd=1.5)
  sig_p ~ dexp(1)

  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mu_phi, sd=sig_phi)
    phi[t] <- expit(logit_phi[t])
  }
  mu_phi ~ dnorm(0,sd=1.5)
  sig_phi ~ dexp(1)

  d[1] <- 1
  beta[1] <- 1
  for(t in 2:K){
    log_f[t-1] ~ dnorm(mu_f, sd=sig_f)
    f[t-1] <- exp(log_f[t-1])
    beta[t] <- d[t-1] * f[t-1]
    d[t] <- d[t-1] * (phi[t-1]+f[t-1])
  }
  mu_f ~ dnorm(0,sd=1.5)
  sig_f ~ dexp(1)

  for(t in 1:(K-1)){
    xi[t] <- beta[t]/sum(beta[t:K])
  }
  xi[K] <- 1

  psi ~ dunif(0,1)
  # Individual processes and observations
  for(i in 1:M) {
    # First occasion
    w[i] ~ dbern(psi)
    z[i, 1] ~ dbern(xi[1])
    y[i, 1] ~ dbern(z[i, 1] * p[1] * w[i])
    r[i,1] <- 1
    # Later occasions
    for(t in 2:K) {
      r[i,t] <- r[i,(t-1)] * (1 - z[i, t-1])
      muz[i, t] <- z[i, t-1]*phi[t-1] + r[i,t]*xi[t]
      z[i, t] ~ dbern(muz[i, t])
      y[i, t] ~ dbern(z[i, t] * p[t] * w[i])
    }
  }
  # Derived quantities
  Nsuper <- sum(w[1:M])
  for(t in 1:K) {
    N[t] <- sum(z[1:M, t] * w[1:M])
  }
})

#' -----------------------------------------------------------------------------
#' Initial values
#' -----------------------------------------------------------------------------
inits <- list(
  logit_p = rep(qlogis(0.5), K),
  logit_phi = rep(qlogis(0.7), K-1),
  log_f = rep(0, K-1),
  mu_phi=0, mu_p=0, mu_f=0,
  sig_phi=1, sig_p=1, sig_f=1,
  psi=0.5,
  z = {
    zinit <- matrix(0, M, K)
    for(i in 1:n.obs) {
      pos <- which(y_aug[i,1:K]==1)
      zinit[i,pos[1]:pos[length(pos)]] <- 1
    }
    zinit
  },
  w = 1.0*(apply(y_aug, 1, sum)>0)
)

#' -----------------------------------------------------------------------------
#' Build, compile, run MCMC
#' -----------------------------------------------------------------------------
Rmodel <- nimbleModel(pxda_code,
                      constants = constants,
                      data = data,
                      inits = inits, check = TRUE)
conf <- configureMCMC(Rmodel,
                      monitors=c("p","mu_p","sig_p",
                                 "phi","mu_phi","sig_phi",
                                 "f", "mu_f","sig_f",
                                 "psi","Nsuper", "N")
                      )
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#' -----------------------------------------------------------------------------
#' Check initial convergence
#' -----------------------------------------------------------------------------
# st <- Sys.time()
# set.seed(8675309)
# samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 0, nchains = 3,
#                        thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
# gelman.diag(samples, autoburnin = FALSE)
# Sys.time()-st

#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------
set.seed(8675309)
samples <- runMCMC(Cmcmc, niter = 60000, nburnin = 10000, nchains = 1,
                       thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
samples_list <- as.list(Cmcmc$mvSamples)

#' -----------------------------------------------------------------------------
#' Summarize MCMC
#' -----------------------------------------------------------------------------

summary(samples)
HPDinterval(samples)
