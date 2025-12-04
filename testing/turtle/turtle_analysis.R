library(marked)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(here)

local_wd <- file.path(here(), "testing", "turtle")
setwd(local_wd)

#'------------------------------------------------------------------------------
#' Get weekly turtle detections
#' -----------------------------------------------------------------------------

honu <- read.table("weekly_detection.txt") #pull in detection txt
honu <- honu[, colnames(honu) != "TURTLEID"] #remove turtle ID column
honu <- unname(honu) #remove header
rownames(honu) <- NULL
honu <- as.matrix(honu)
honu <- honu + 1


#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("turtle_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(K=ncol(honu), nobs=nrow(honu)),
  data = list(x = honu, n = nrow(honu)),
  inits = list(
    p_dot = 0.5,
    phi_dot = 0.5,
    rho_dot = 0.5,
    psi23_dot = 0.5,
    psi32_dot = 0.5,
    lambda=2*nrow(honu)
    )
)
c_js_model <- compileNimble(js_model)
js_mcmc <- buildMCMC(js_model,
                     monitors=c("p_dot","phi_dot","psi23_dot","psi32_dot",
                                "Nsuper","N_nest", "N_internest","N"))
c_js_mcmc <- compileNimble(js_mcmc)


#' -----------------------------------------------------------------------------
#' Check initial convergence
#' -----------------------------------------------------------------------------
set.seed(8675309)
samples <- runMCMC(c_js_mcmc, niter = 10000, nburnin = 0, nchains = 3,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
gelman.diag(samples, autoburnin = FALSE)


#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------

set.seed(8675309)
samples <- runMCMC(c_js_mcmc, niter = 5000, nburnin = 0, nchains = 1,
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

summary(samples_list$rho_dot)
HPDinterval(samples_list$rho_dot)

summary(samples_list$phi_dot)
HPDinterval(samples_list$phi_dot)

summary(samples_list$p_dot)
HPDinterval(samples_list$p_dot)



Ndf <- data.frame(week = 1:ncol(honu), est=colMeans(samples_list$N), hpd = HPDinterval(mcmc(samples_list$N)))
ggplot(Ndf) + geom_point(aes(x=week, y=est), size=1) +
  geom_errorbar(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=week, y=est)) + ylab("Abundance") + xlab("Week") +
  theme_bw()
# ggsave("dipper_Nt.pdf", width=6.5, height=4)

Ndf <- data.frame(week = 1:ncol(honu), est=colMeans(samples_list$N_nest), hpd = HPDinterval(mcmc(samples_list$N_nest)))
ggplot(Ndf) + geom_point(aes(x=week, y=est), size=1) +
  geom_errorbar(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=week, y=est)) + ylab("Nesting individuals") + xlab("Week") +
  theme_bw()
ggsave("dipper_Nt.pdf", width=6.5, height=4)





