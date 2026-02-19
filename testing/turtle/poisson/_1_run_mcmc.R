library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

setwd(here("testing","turtle","poisson"))

source("helper/make_hmm_mats.R")

honu_ch <- read_csv("honu_ch_pois.csv") %>% as.matrix()

# Years 2013 and 2020 were not surveyed
surv <- colnames(honu_ch)
surv <- ifelse(surv %in% c("X2013","X2020"), 0, 1)
surv <- cbind(surv,surv)
surv[colnames(honu_ch)=="X2015",1] <- 0


# 5 animals return to breeding without a foraging period so they are removed
skip_forage <- function(x){
  rrr <- rle(1-(rle(x)$values==0))
  out <- max(rrr$lengths)
  return(out>1)
}
honu_ch <- honu_ch[as.logical(1-apply(honu_ch, 1, skip_forage)),]
# honu_ch <- honu_ch + 1


ch_collapsed <- collapse_ch(honu_ch)
K <- ncol(honu_ch)
n <- nrow(honu_ch)

# Radial basis matrix
B <- make_GRBF(1:K, 12)


source("honu_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, nd=7, n_states=11, m=ncol(B),
    nobs_u=nrow(ch_collapsed$unique_ch)#, nobs=nrow(honu_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, surv=surv, #x=honu+ch,
    n=n, B=B
  ),
  inits = list(
    #xi_raw = rep(1/K,K),
    mu_rho=0, beta = rep(0,ncol(B)), sig_beta=2,
    mu_phi=qlogis(0.9), delta = rep(0,ncol(B)), sig_delta=2,
    mu_theta=0, gamma = rep(0,ncol(B)), sig_gamma=2,
    # theta=2,
    mu_r1=log(6), tau1=rep(0,K), sig_tau1=2,
    mu_r2=log(6), tau2=rep(0,K), sig_tau2=2,
    lambda=2*n
  )
)

c_js_model <- compileNimble(js_model)
mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            "rho", "xi", "mu_rho", "beta", "sig_beta",
                            "theta",
                            "phi", "mu_phi", "delta", "sig_delta",
                            "r",
                            "mu_r1", "tau1", "sig_tau1",
                            "mu_r2", "tau2", "sig_tau2",
                            "nu","Nsuper","lambda", "pstar"
                          )
)
mcmcConf$removeSamplers('beta')
mcmcConf$addSampler(target = 'beta', type = 'RW_block')
mcmcConf$removeSamplers('delta')
mcmcConf$addSampler(target = 'delta', type = 'RW_block')
mcmcConf$removeSamplers('gamma')
mcmcConf$addSampler(target = 'gamma', type = 'RW_block')
mcmcConf$removeSamplers('tau1')
mcmcConf$addSampler(target = 'tau1', type = 'RW_block')
mcmcConf$removeSamplers('tau2')
mcmcConf$addSampler(target = 'tau2', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 30000, nburnin = 10000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st

samples_list <- as.list(c_js_mcmc$mvSamples)
save(samples_list, honu_ch, ch_collapsed, surv, B, file="honu_mcmc_sample.RData")

