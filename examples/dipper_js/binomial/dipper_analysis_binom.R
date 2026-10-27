library(marked)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(modeest)
library(here)

setwd(here("dipper_js", "binomial"))

#'------------------------------------------------------------------------------
#' Get dipper data from R package 'marked'
#' -----------------------------------------------------------------------------
data("dipper")

#' Convert capture histories to binary matrix: 0 = not captured, 1 = captured
x <- strsplit(dipper$ch, "") %>% lapply(., as.numeric) %>% do.call(rbind, .)

K <- ncol(x)
n <- nrow(x)

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dipper_model_binom.R")

# Build nimble model and add constants/data/inits
js_model <- nimbleModel(
  code = js_code,
  constants = list(K = K, nobs = n, n_states = 3),
  data = list(
    x = x,
    n = n,
    size_vec = rep(1, K)  # Vector of ones for size
  ),
  inits = list(
    logit_p = rep(qlogis(0.5), K),
    logit_phi = rep(qlogis(0.7), K - 1),
    log_f = rep(0, K - 1),
    mu_phi = 0, mu_p = 0, mu_f = 0,
    sig_phi = 1, sig_p = 1, sig_f = 1,
    lambda = 2 * nrow(x)
  )
)

# Compile model
c_js_model <- compileNimble(js_model)

# Configure MCMC
mcmcConf <- configureMCMC(
  js_model, 
  monitors = c("p", "mu_p", "sig_p",
               "phi", "mu_phi", "sig_phi",
               "f", "mu_f", "sig_f",
               "lambda", "Nsuper", "N")
)

# Explicitly configure z and nu_t as posterior predictive samplers
mcmcConf$removeSamplers(c('z', 'nu_t'))
mcmcConf$addSampler(target = 'z', type = 'sampler_posterior_predictive')
mcmcConf$addSampler(target = 'nu_t', type = 'sampler_posterior_predictive')

# Build and compile MCMC
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)

#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------
set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 60000, nburnin = 10000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et - st

# Extract MCMC quantities as a named list
samples_list <- as.list(c_js_mcmc$mvSamples)

#' -----------------------------------------------------------------------------
#' Summarize MCMC
#' -----------------------------------------------------------------------------
summary(samples)
HPDinterval(samples)

# Yearly abundance estimates
Ndf <- data.frame(
  year = 1:ncol(x), 
  est = colMeans(samples_list$N), 
  se = apply(samples_list$N, 2, sd),
  hpd = HPDinterval(mcmc(samples_list$N))
)
print(Ndf)

# Plot abundance trajectory
ggplot(Ndf) + 
  geom_point(aes(x = year, y = est), size = 3) +
  geom_errorbar(aes(x = year, ymin = hpd.lower, ymax = hpd.upper), width = 0.2) +
  geom_path(aes(x = year, y = est)) + 
  ylab("Abundance") + xlab("Year") +
  theme_bw()
