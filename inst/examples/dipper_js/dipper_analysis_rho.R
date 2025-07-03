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
#' Observed states: 1 = not captured, 2 = captured
x <- x+1


#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("js_nimble_model_rho.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(K=ncol(x), nobs=nrow(x)),
  data = list(x = x, n = nrow(x)),
  inits = list(
    rho=c(rep(0.5,6),NA), p=rep(0.5,7), phi = rep(0.5,6), lambda=1.1*nrow(x)
    )
)
c_js_model <- compileNimble(js_model)

js_mcmc <- buildMCMC(
  js_model,
  monitors=c("rho","phi","p","lambda","nu","Nsuper","Nu", "Nd", "N", "nu_t", "xi","beta")
  )
c_js_mcmc <- compileNimble(js_mcmc)

samples <- runMCMC(c_js_mcmc, niter = 25000, nburnin = 5000, nchains = 1, thin2=1, thin=1)
samples_list <- as.list(c_js_mcmc$mvSamples)


summary(mcmc(samples_list$phi))
summary(mcmc(samples_list$rho))
summary(mcmc(samples_list$p))
summary(mcmc(samples_list$N))
summary(mcmc(samples_list$beta))
summary(mcmc(samples_list$Nd))

Ndf <- data.frame(year = 1:ncol(x), est=apply(samples_list$N,2,median), hpd = HPDinterval(mcmc(samples_list$N)))

ggplot(Ndf2) + geom_point(aes(x=year, y=est), size=3) +
  # geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year")




