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
#' Get dipper data from marked
#' -----------------------------------------------------------------------------
data("dipper")
#' Convert capture histories to a matrix
x <- strsplit(dipper$ch, "") %>% lapply(.,as.numeric) %>% do.call(rbind,.)


#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dipper_model_binom.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(K=ncol(x), mu_beta=rep(1, ncol(x)), nobs=nrow(x)),
  data = list(x = x, n = nrow(x), ones=rep(1,ncol(x))),
  inits = list(
    beta=rep(1/7,7), p=c(NA,rep(0.5,5),NA), phi = rep(0.5,6), lambda=1.1*nrow(x)
    )
)
c_js_model <- compileNimble(js_model)

js_mcmc <- buildMCMC(
  js_model,
  monitors=c("beta","phi","p","lambda","nu","pstar","Nsuper", "Nd", "N")
  )
c_js_mcmc <- compileNimble(js_mcmc)

samples <- runMCMC(c_js_mcmc, niter = 25000, nburnin = 5000, nchains = 1, thin2=1, thin=1)
samples_list <- as.list(c_js_mcmc$mvSamples)


summary(mcmc(samples_list$phi))
summary(mcmc(samples_list$beta))
summary(mcmc(samples_list$p))
summary(mcmc(samples_list$N))
summary(mcmc(samples_list$Nsuper))

Ndf <- data.frame(year = 1:ncol(x), est=apply(samples_list$N, 2, median), hpd = HPDinterval(mcmc(samples_list$N)))

ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year")




