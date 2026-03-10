library(marked)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(here)

local_wd <- file.path(here(), "inst", "examples", "dipper_js","parallel")
setwd(local_wd)

#'------------------------------------------------------------------------------
#' Get dipper data from R package 'marked'
#' -----------------------------------------------------------------------------
data("dipper")
#' Convert capture histories to a matrix
x <- strsplit(dipper$ch, "") %>% lapply(.,as.numeric) %>% do.call(rbind,.)
#' Observed states: 1 = not captured, 2 = captured
x <- x+1

#' Collapse to unique capture histories
ux <- collapse_ch(x)

K <- ncol(x)
n <- nrow(x)

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dipper_model_parallel.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, nobs_u=nrow(ux$unique_ch)
    ),
  data = list(
    ux = ux$unique_ch, w=ux$w, n=n
    ),
  inits = list(
    logit_p = rep(qlogis(0.5), K),
    logit_phi = rep(qlogis(0.7), K-1),
    log_f = rep(0, K-1),
    mu_phi=0, mu_p=0, mu_f=0,
    sig_phi=1, sig_p=1, sig_f=1,
    lambda=2*nrow(x)
    )
)
c_js_model <- compileNimble(js_model)
js_mcmc <- buildMCMC(js_model,
                     monitors=c("p","mu_p","sig_p",
                                "phi","mu_phi","sig_phi",
                                "f", "mu_f","sig_f", "xi",
                                "lambda","Nsuper", "nu"))
c_js_mcmc <- compileNimble(js_mcmc)


#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------

st <- Sys.time()
set.seed(8675309)
samples <- runMCMC(c_js_mcmc, niter = 60000, nburnin = 10000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et - st

samples_list <- as.list(c_js_mcmc$mvSamples)

#' -----------------------------------------------------------------------------
#' Summarize MCMC
#' -----------------------------------------------------------------------------

lambda <- mcmc(samples_list$lambda)
summary(lambda)
HPDinterval(lambda)

Nsuper <- mcmc(samples_list$Nsuper)
summary(Nsuper)
HPDinterval(Nsuper)

f <- mcmc(samples_list$f)
summary(f)
HPDinterval(f)

phi <- mcmc(samples_list$phi)
summary(phi)
HPDinterval(phi)

p <- mcmc(samples_list$p)
summary(p)
HPDinterval(p)

re_par <- samples[,c("mu_f","sig_f","mu_phi","sig_phi","mu_p","sig_p")]
summary(re_par)
HPDinterval(re_par)


save(samples_list, x, file=file.path(local_wd, "dipper_stage1_sample.RData"))






