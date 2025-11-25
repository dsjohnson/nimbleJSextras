library(marked)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(here)

local_wd <- file.path(here(), "inst", "examples", "dipper_js")
setwd(local_wd)

#'------------------------------------------------------------------------------
#' Get dipper data from R package 'marked'
#' -----------------------------------------------------------------------------
data("dipper")
#' Convert capture histories to a matrix
x <- strsplit(dipper$ch, "") %>% lapply(.,as.numeric) %>% do.call(rbind,.)
K <- ncol(x)

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dipper_model_bern.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(K=ncol(x), nobs=nrow(x)),
  data = list(x = x, n = nrow(x)),
  inits = list(
  logit_p = rep(qlogis(0.5), K),
  logit_phi = rep(qlogis(0.7), K-1),
  logit_rho = rep(qlogis(0.1), K),
  psi = 0.5,
  mu_phi=0, mu_p=0, mu_rho=0,
  sig_phi=1, sig_p=1, sig_rho=1,
  lambda=2*nrow(x)
  )
)
c_js_model <- compileNimble(js_model)
js_mcmc <- buildMCMC(js_model,
                     monitors=c("p","mu_p","sig_p",
                                "rho", "mu_rho","sig_rho",
                                "phi","mu_phi","sig_phi",
                                "lambda","Nsuper", "N"))
c_js_mcmc <- compileNimble(js_mcmc)


#' -----------------------------------------------------------------------------
#' Check initial convergence
#' -----------------------------------------------------------------------------
# set.seed(8675309)
# samples <- runMCMC(c_js_mcmc, niter = 5000, nburnin = 0, nchains = 3,
#                    thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
# gelman.diag(samples, autoburnin = FALSE)


#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------

set.seed(8675309)
samples <- runMCMC(c_js_mcmc, niter = 60000, nburnin = 10000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
samples_list <- as.list(c_js_mcmc$mvSamples)
samples_list <- lapply(samples_list, mcmc)

#' -----------------------------------------------------------------------------
#' Summarize MCMC
#' -----------------------------------------------------------------------------

summary(samples_list$lambda)
HPDinterval(samples_list$lambda)

summary(samples_list$Nsuper)
HPDinterval(samples_list$Nsuper)

summary(samples_list$rho)
HPDinterval(samples_list$rho)

summary(samples_list$phi)
HPDinterval(samples_list$phi)

summary(samples_list$p)
HPDinterval(samples_list$p)

summary(samples_list$mu_rho)
HPDinterval(samples_list$mu_rho)

summary(samples_list$sig_rho)
HPDinterval(samples_list$sig_rho)

summary(samples_list$mu_phi)
HPDinterval(samples_list$mu_phi)

summary(samples_list$sig_phi)
HPDinterval(samples_list$sig_phi)

summary(samples_list$mu_p)
HPDinterval(samples_list$mu_p)

summary(samples_list$sig_p)
HPDinterval(samples_list$sig_p)

Ndf <- data.frame(year = 1:ncol(x), est=colMeans(samples_list$N), hpd = HPDinterval(mcmc(samples_list$N)))
ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year") +
  theme_bw()
ggsave("dipper_Nt.pdf", width=6.5, height=4)






