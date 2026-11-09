library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

# Set working directory.

setwd(here("honu_js"))

source(file.path("R","make_hmm_mats.R"))
source(file.path("R","matern_rbf.R"))

honu_ch <- read_csv(file.path("raw_data","honu_ch.csv")) %>% as.matrix()

# 2013 was not surveyed; 2015 no new perm IDs put out
surv <- colnames(honu_ch)
surv <- ifelse(surv == "X2013", 0, 1)

ch_collapsed <- collapse_ch(honu_ch)

K <- ncol(honu_ch)
n <- nrow(honu_ch)
m <- 5 # number of breeding return time parameters

nu <- 2.5
knots <- make_padded_knots(1:K, n_interior = 10, nu = nu)$knots
B <-  matern_basis(1:K, knots, nu=nu)


source(file.path("R", "honu_model.R"))

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, m=m, n_states=6+2*m,
    nobs_u=nrow(ch_collapsed$unique_ch),
    df_xi = length(knots), B = B
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, n=n, surv=surv
  ),
  inits = list(
    sigma_beta_xi=10,
    beta_xi = rep(0,ncol(B)),
    logit_phi=qlogis(0.92),
    logit_tau = rep(0, m),
    logit_p=rep(qlogis(0.9),K), logit_eta=rep(2.5,K-2),
    lambda=2*n
  )
)
mcmcConf <- configureMCMC(js_model,
                          useConjugacy = FALSE,
                          monitors=c(
                            "sigma_beta_xi","beta_xi", "xi",
                            "phi", "p", "p0", "tau",
                            "ft_pdf", "Nu","N","lambda", "pstar"
                          ))
mcmcConf$removeSamplers('beta_xi')
mcmcConf$addSampler(target = 'beta_xi', type = 'RW_block')
mcmcConf$removeSamplers('logit_tau')
mcmcConf$addSampler(target = 'logit_tau', type = 'RW_block')
mcmcConf$removeSamplers('logit_eta')
mcmcConf$addSampler(target = 'logit_eta', type = 'RW_block')
mcmcConf$removeSamplers('logit_p')
mcmcConf$addSampler(target = 'logit_p', type = 'RW_block')

js_mcmc <- buildMCMC(mcmcConf)

# Compile model and MCMC sampler
c_js_model <- compileNimble(js_model)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(42)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 300000, nburnin = 50000, nchains = 1,
                   thin = 5, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st
# Time difference of 3.467576 hours

samples_list <- as.list(c_js_mcmc$mvSamples)
# save(samples_list, honu_ch, ch_collapsed, surv,
#      file=file.path("processed_data_samples", "honu_mcmc_sample.RData"))

