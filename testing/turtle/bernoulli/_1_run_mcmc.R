library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

setwd(here("testing","turtle","bernoulli"))

source("helper/make_hmm_mats.R")

honu_ch <- read_csv("honu_ch.csv") %>% as.matrix()

# 2013 was not surveyed; 2015 no new perm IDs put out
surv <- colnames(honu_ch)
surv <- ifelse(surv %in% c("X2013","X2020"), 0, 1)

ch_collapsed <- collapse_ch(honu_ch)
K <- ncol(honu_ch)
n <- nrow(honu_ch)

# Radial basis matrix
# B <- make_GRBF(1:K, 12)


source("honu_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, nd=7, n_states=6+2*7, #m=ncol(B),
    nobs_u=nrow(ch_collapsed$unique_ch)#, nobs=nrow(honu_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, surv=surv, #x=honu+ch,
    n=n#, B=B
  ),
  inits = list(
    #xi_raw = rep(1/K,K),
    alpha=0,
    rho_dot=1/(1+K),
    #mu_rho=qlogis(1/K), #beta = rep(0,ncol(B)), sig_beta=2,
    mu_phi=qlogis(0.9), #delta = rep(0,ncol(B)), sig_delta=2,
    mu_theta=log(2), #gamma = rep(0,ncol(B)), sig_gamma=2,
    # theta=2,
    mu_p1=qlogis(0.5), tau1=rep(0,K), sig_tau1=2,
    mu_p2=qlogis(0.9), tau2=rep(0,K), sig_tau2=2,
    lambda=2*n
  )
)
mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            # "alpha",
                            "rho", "xi", #"mu_rho", "beta", "sig_beta",
                            "theta",
                            "phi", "mu_phi", #"delta", "sig_delta",
                            "p",
                            "mu_p1", "tau1", "sig_tau1",
                            "mu_p2", "tau2", "sig_tau2",
                            "nu","Nsuper","lambda", "pstar"
                          )
)
# mcmcConf$removeSamplers('beta')
# mcmcConf$addSampler(target = 'beta', type = 'RW_block')
# mcmcConf$removeSamplers('delta')
# mcmcConf$addSampler(target = 'delta', type = 'RW_block')
# mcmcConf$removeSamplers('gamma')
# mcmcConf$addSampler(target = 'gamma', type = 'RW_block')
mcmcConf$removeSamplers('tau1')
mcmcConf$addSampler(target = 'tau1', type = 'RW_block')
mcmcConf$removeSamplers('tau2')
mcmcConf$addSampler(target = 'tau2', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)

c_js_model <- compileNimble(js_model)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 25000, nburnin = 5000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st

samples_list <- as.list(c_js_mcmc$mvSamples)
save(samples_list, honu_ch, ch_collapsed, surv, B, file="honu_mcmc_sample.RData")

