library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(here)

local_wd <- file.path(here(), "inst", "examples", "dolphin")
setwd(local_wd)

#'------------------------------------------------------------------------------
#' Get dipper data from marked
#' -----------------------------------------------------------------------------
dolphin_df <- read_csv("dolphin_captures_rd.csv") %>% filter(occasion_start!="2020-11-28")

n <- length(unique(dolphin_df$id))
K <- nrow(dolphin_df)/n
x <- matrix(dolphin_df$capture, nrow=n, ncol=K, byrow=T)
Rt <- dolphin_df$n_surv[1:K]
surveyed <- ifelse(Rt>0,1,0)
n_surveys <- sum(surveyed)
survey_num <- cumsum(surveyed) * surveyed
survey_num <- ifelse(survey_num=="0", NA, survey_num)
Midx <- factor(approx(x=1:K, y=survey_num, xout=1:K, method="constant", f=1)$y)
Mxi <- model.matrix(~0+Midx)
Mp <- sweep(Mxi, 1, surveyed, "*")

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dolphin_js_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(K=ncol(x), nobs=nrow(x), n_surveys=n_surveys),
  data = list(x = x, n = nrow(x), Rt=Rt, surveyed=surveyed, Mxi=Mxi, Mp=Mp),
  inits = list(
    # beta_xi_R=c(NA,rep(0,n_surveys-1)), beta_xi_P=c(NA,rep(0,n_surveys-1)),beta_xi_T=c(NA,rep(0,n_surveys-1)),
    # sig_xi_R = 1, sig_xi_P=1, sig_xi_T=1,
    rho_R = c(0.5,0.5), rho_P = c(0.5,0.5), rho_T = c(0.5,0.5),
    phi_RP=0.85, delta_phi = 0.5,
    mlogit_alpha=c(NA,0.35,1.12),
    mu_p=0, logit_p_RT = rep(0,n_surveys), delta_p=0.5, lambda=2*n
  )
)
c_js_model <- compileNimble(js_model)

mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            "alpha", "rho_R", "rho_P", "rho_T",
                            "phi_RP","delta_phi",
                            "p_RT","delta_p", "eps", "mu_p", "sig_p",
                            "lambda",
                            "N_R", "N_P", "N_T","Nsuper"
                          ))
mcmcConf$removeSamplers('eps')
mcmcConf$addSampler(target = 'eps', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)

samples <- runMCMC(c_js_mcmc, niter = 25000, nburnin = 5000, nchains = 1, thin2=1, thin=1)
samples_list <- as.list(c_js_mcmc$mvSamples)
samples_list <- lapply(samples_list, mcmc)


summary(samples_list$Nsuper)

summary(samples_list$rho_R)
summary(samples_list$rho_P)
summary(samples_list$rho_T)

summary(samples_list$phi_RP)
summary(samples_list$delta_phi * samples_list$phi_RP)

summary(samples_list$p_RT)
summary(samples_list$delta_p)

summary(samples_list$alpha)

Ndf <- data.frame(
  occ = 1:K,
  est_R=colMeans(samples_list$N_R),
  hpd_R = HPDinterval(samples_list$N_R),
  est_P=colMeans(samples_list$N_P),
  hpd_P = HPDinterval(samples_list$N_P),
  est_T=colMeans(samples_list$N_T),
  hpd_T = HPDinterval(samples_list$N_T)
)

ggplot(Ndf) + #geom_point(aes(x=year, y=est), size=3) +
  geom_ribbon(aes(x=occ, ymin=hpd_R.lower, ymax=hpd_R.upper), fill="blue", alpha=0.2) +
  geom_path(aes(x=occ, y=est_R), color="blue") +
  geom_ribbon(aes(x=occ, ymin=hpd_P.lower, ymax=hpd_P.upper), fill="darkred", alpha=0.2) +
  geom_path(aes(x=occ, y=est_P), color="darkred") +
  geom_ribbon(aes(x=occ, ymin=hpd_R.lower, ymax=hpd_R.upper), fill="gold", alpha=0.2) +
  geom_path(aes(x=occ, y=est_T), color="gold") +
  ylab("Abundance") + xlab("Occasion")




