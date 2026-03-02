library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

setwd(here("testing","turtle","bernoulli"))

source("helper/make_hmm_mats.R")

honu_ch <- read_csv("honu_ch.csv") %>% as.matrix()
effort <- read_csv("raw/survey_summary.csv") %>% filter(year<=2016)

# 2013 was not surveyed; 2015 no new perm IDs put out
surv <- colnames(honu_ch)
surv <- ifelse(surv %in% c("X2013","X2020"), 0, 1)

ch_collapsed <- collapse_ch(honu_ch)
K <- ncol(honu_ch)
n <- nrow(honu_ch)

# design matrices
B <- cbind(1,make_GRBF(1:K, 5))
# X <- model.matrix(~0+total_days, data=effort); X <- X/120
X <- cbind(1,effort$total_days/120)

source("honu_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, nd=7, n_states=6+2*7, m=ncol(B), #mp=2,
    nobs_u=nrow(ch_collapsed$unique_ch)#, nobs=nrow(honu_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, surv=surv, #x=honu_ch,
    n=n, B=B, X=X
  ),
  inits = list(
    beta_rho_int=0, beta_rho=rep(0,ncol(B)), sig_rho=2,
    beta_phi_int=0, beta_phi=rep(0,ncol(B)), sig_phi=2,
    beta_theta_int=0, beta_theta=rep(0,ncol(B)), sig_theta=2,
    beta_p = rep(0,ncol(X)), sig_p=2,
    alpha=rep(5/6,K),
    lambda=2*n
  )
)
mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            "rho", "xi","beta_rho_int","beta_rho", "sig_rho",
                            "phi","beta_phi_int","beta_phi","sig_phi",
                            "theta","beta_theta_int","beta_theta","sig_theta",
                            "p", "alpha", "beta_p", "sig_p",
                            "nu","Nsuper","lambda", "pstar"
                          ))
mcmcConf$removeSamplers('beta_rho')
mcmcConf$addSampler(target = 'beta_rho', type = 'RW_block')
mcmcConf$removeSamplers('beta_phi')
mcmcConf$addSampler(target = 'beta_phi', type = 'RW_block')
mcmcConf$removeSamplers('beta_theta')
mcmcConf$addSampler(target = 'beta_theta', type = 'RW_block')
mcmcConf$removeSamplers('beta_p')
mcmcConf$addSampler(target = 'beta_p', type = 'RW_block')
# mcmcConf$removeSamplers('tau')
# mcmcConf$addSampler(target = 'tau', type = 'RW_block')
mcmcConf$removeSamplers('alpha')
mcmcConf$addSampler(target = 'alpha', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)

# Compile model and MCMC sampler
c_js_model <- compileNimble(js_model)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 65000, nburnin = 15000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st

samples_list <- as.list(c_js_mcmc$mvSamples)
save(samples_list, honu_ch, ch_collapsed, surv, file="honu_mcmc_sample.RData")

