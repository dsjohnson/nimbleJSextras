library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

setwd(here("testing","turtle","final_analysis"))

source("helper/make_hmm_mats2.R")

honu_ch <- read_csv("honu_ch.csv") %>% as.matrix()
# effort <- read_csv("raw/survey_summary.csv") %>% filter(year<=2016)

# 2013 was not surveyed; 2015 no new perm IDs put out
surv <- colnames(honu_ch)
surv <- ifelse(surv %in% c("X2013","X2020"), 0, 1)

ch_collapsed <- collapse_ch(honu_ch)
K <- ncol(honu_ch)
n <- nrow(honu_ch)

m <- 5

# design matrices
# B_rho <- make_GRBF(1:K, 20)$basis; k_rho <- ncol(B_rho)
# B_mu <- make_GRBF(1:(K-1), 5)$basis; k_mu <- ncol(B_mu)
# B_phi <- make_GRBF(1:(K-1), 5)$basis; k_phi <- ncol(B_phi)


source("honu_model2.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, m=m, n_states=6+2*m,
    #k_rho=k_rho, k_mu=k_mu, #k_phi=k_phi,
    nobs_u=nrow(ch_collapsed$unique_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, surv=surv, #x=honu_ch,
    n=n#, B_rho=B_rho, B_mu=B_mu#, B_phi=B_phi,
  ),
  inits = list(
    #beta_rho_int=0, sig_rho=5, beta_rho=rep(0,k_rho), #beta_rho=rep(0,k_rho),
    # rho=(1/K)/(1/K + 2-c(1:K)/K),
    eps_rho=qlogis(rep(1/K, K)),
    # beta_phi_int=0, beta_phi=rep(0,k_phi), sig_phi=5,
    eps_phi=qlogis(0.92),
    # beta_mu_int=0, sig_mu=5, beta_mu=rep(0,k_mu), #beta_mu=rep(0,k_mu),
    eps_tau = rep(0, m),
    # beta_p_int=0, eps_p=rep(0,K), sig_p=5,
    eps_p=rep(qlogis(0.9),K-1),
    eps_a=rep(0,K),
    lambda=2*n
  )
)
mcmcConf <- configureMCMC(js_model,
                          useConjugacy = FALSE,
                          monitors=c(
                            "xi","rho",#"beta_rho_int", "beta_rho", "sig_rho",
                            "phi",#"beta_phi_int","beta_phi","sig_phi",
                            #"mu","beta_mu_int","sig_mu", "beta_mu",
                            "tau","ft_pdf",
                            "p", "alpha", #"beta_p_int", "eps_p", "sig_p",
                            "nu","Nsuper","lambda", "pstar"
                          ))
# mcmcConf$removeSamplers('rho')
# mcmcConf$addSampler(target = 'rho', type = 'RW_block')
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

