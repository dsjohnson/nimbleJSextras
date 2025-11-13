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
n.aug <- 1000  # choose augmented size (tune as needed)
M <- n.obs + n.aug
y_aug <- matrix(0, nrow = M, ncol = K)
y_aug[1:n.obs, ] <- y.obs

constants <- list(M = M, K = K)
data <- list(y = y_aug)

#'------------------------------------------------------------------------------
#' PX-DA code
#' -----------------------------------------------------------------------------
pxda_code <- nimbleCode({
  # Priors
  for(t in 1:K){
    logit_p[t] ~ dnorm(mu_p, sd=sig_p)
    p[t] <- expit(logit_p[t])
    logit_rho[t] ~ dnorm(mu_rho, sd=sig_rho)
    rho[t] <- expit(logit_rho[t])
  }
  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mu_phi, sd=sig_phi)
    phi[t] <- expit(logit_phi[t])
  }
  sig_p ~ dexp(1)
  mu_p ~ dnorm(0,sd=1.5)
  sig_rho ~ dexp(1)
  mu_rho ~ dnorm(0,sd=10)
  mu_phi ~ dnorm(0,sd=1.5)
  sig_phi ~ dexp(1)

  # conditional recruitment
  xi[1] <- rho[1]
  for(t in 2:(K-1)){ xi[t] <- rho[t]/(1-prod(1-rho[t:K])) }
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
  logit_rho = rep(qlogis(0.1), K),
  psi = 0.5,
  mu_phi=0, mu_p=0, mu_rho=0,
  sig_phi=1, sig_p=1, sig_rho=1,
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
                      monitors = c("psi", "p", "phi", "rho", "N", "Nsuper","mu_phi","sig_phi","mu_p","sig_p","mu_rho","sig_rho")
                      )
Rmcmc <- buildMCMC(conf)
Cmodel <- compileNimble(Rmodel)
Cmcmc <- compileNimble(Rmcmc, project = Rmodel)

#' -----------------------------------------------------------------------------
#' Check initial convergence
#' -----------------------------------------------------------------------------
# set.seed(8675309)
# samples <- runMCMC(Cmcmc, niter = 5000, nburnin = 0, nchains = 3,
#                        thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
# gelman.diag(samples, autoburnin = FALSE)

#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------
samples <- runMCMC(Cmcmc, niter = 60000, nburnin = 10000, nchains = 1,
                       thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)

samples_list <- as.list(Cmcmc$mvSamples)

#' -----------------------------------------------------------------------------
#' Summarize MCMC
#' -----------------------------------------------------------------------------

summary(mcmc(samples_list$N))
summary(mcmc(samples_list$Nsuper))
summary(mcmc(samples_list$rho))
summary(mcmc(samples_list$phi))
summary(mcmc(samples_list$p))
summary(mcmc(samples_list$mu_rho))
summary(mcmc(samples_list$sig_rho))
summary(mcmc(samples_list$mu_phi))
summary(mcmc(samples_list$sig_phi))
summary(mcmc(samples_list$mu_p))
summary(mcmc(samples_list$sig_p))
