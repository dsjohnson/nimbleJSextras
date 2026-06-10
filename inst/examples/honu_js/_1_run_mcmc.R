library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

wd <- file.path(system.file(package="nimbleJSextras"), "examples", "honu_js")

setwd(here("isnt", "examples", "honu_js"))

source("helper/make_hmm_mats.R")

honu_ch <- read_csv("honu_ch.csv") %>% as.matrix()

# 2013 was not surveyed; 2015 no new perm IDs put out
surv <- colnames(honu_ch)
surv <- ifelse(surv == "X2013", 0, 1)

ch_collapsed <- collapse_ch(honu_ch)

K <- ncol(honu_ch)
n <- nrow(honu_ch)
m <- 5

source("honu_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, m=m, n_states=6+2*m,
    nobs_u=nrow(ch_collapsed$unique_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, surv=surv
  ),
  inits = list(
    eps_rho=qlogis(rep(1/K, K)), eps_phi=qlogis(0.92), eps_tau = rep(0, m), eps_p=rep(qlogis(0.9),K-1), eps_a=rep(0,K),
    lambda=2*n
  )
)
mcmcConf <- configureMCMC(js_model,
                          useConjugacy = FALSE,
                          monitors=c(
                            "xi", "rho", "phi","p", "alpha", "tau",
							"ft_pdf", "nu","Nsuper","lambda", "pstar"
                          ))
mcmcConf$removeSamplers('eps_rho')
mcmcConf$addSampler(target = 'eps_rho', type = 'RW_block')
mcmcConf$removeSamplers('eps_tau')
mcmcConf$addSampler(target = 'eps_tau', type = 'RW_block')
mcmcConf$removeSamplers('eps_a')
mcmcConf$addSampler(target = 'eps_a', type = 'RW_block')
mcmcConf$removeSamplers('eps_p')
mcmcConf$addSampler(target = 'eps_p', type = 'RW_block')
mcmcConf$removeSamplers('lambda')
mcmcConf$addSampler(target = 'lambda', type = 'conjugate')
js_mcmc <- buildMCMC(mcmcConf)

# Compile model and MCMC sampler
c_js_model <- compileNimble(js_model)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 70000, nburnin = 20000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st

samples_list <- as.list(c_js_mcmc$mvSamples)
save(samples_list, honu_ch, ch_collapsed, surv, file="honu_mcmc_sample.RData")

