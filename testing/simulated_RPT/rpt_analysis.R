library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)
library(here)

local_wd <- file.path(here(), "inst", "examples", "simulated_RPT")
setwd(local_wd)

#'------------------------------------------------------------------------------
#' Get dipper data from marked
#' -----------------------------------------------------------------------------
load("rpt_data.RData")
x <- sim$ch + 1
x_collapsed <- collapse_ch(x)
K <- ncol(x)
n <- nrow(x)

yidx <- rbind(
  c(1,12),
  c(13,24),
  c(25,36)
)

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("rpt_model_code.R")

js_model <- nimbleModel(
  code = rpt_code,
  constants = list(
    K=K,  nobs_u=nrow(x_collapsed$unique_ch)#, nobs=nrow(x), yidx=yidx
    ),
  data = list(
    ux=x_collapsed$unique_ch, w=x_collapsed$w, #x=x,
    n=n#, xi_mean=rep(1/K,K)
    ),
  inits = list(
    # xi_R=rep(1/K,K), xi_P=rep(1/K,K), xi_T=rep(1/K,K),
    rho = matrix(rep((1/K)/(1/K + 2-c(1:K)/K),3),K,3),
    phi_RP=0.75, phi_T=0.01,
    alpha_raw=c(.2, .25, .55),
    mu_p=-1.75, tau=rep(0,K), sig_p=0.01,
    delta_p=0.7,
    lambda=2*n
  )
)
c_js_model <- compileNimble(js_model)

mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            #"xi_R","xi_P","xi_T",
                            "rho","alpha",
                            "phi_RP","phi_T", "phi_tilde",
                            "p", "mu_p","tau", "sig_p","delta_p",
                            "Nsuper", "nu_raw","lambda",
                            "pstar"#,"loglik"#,"N_R_yr", "N_P_yr","N_T_yr"
                          )
)
mcmcConf$removeSamplers('tau')
mcmcConf$addSampler(target = 'tau', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)
samples <- runMCMC(c_js_mcmc, niter = 10000, nburnin = 5000, nchains = 1, thin2=1, thin=1)



samples_list <- as.list(c_js_mcmc$mvSamples) %>% lapply(mcmc)

rho_R <- as.list(c_js_mcmc$mvSamples)$rho[,,1] %>% mcmc
rho_P <- as.list(c_js_mcmc$mvSamples)$rho[,,2] %>% mcmc
rho_T <- as.list(c_js_mcmc$mvSamples)$rho[,,3] %>% mcmc


summary(samples_list$Nsuper)
plot(samples_list$Nsuper)
plot(-samples_list$loglik)
summary(samples_list$alpha)
plot(x=1:36, colMeans(rho_R))
plot(x=1:36, colMeans(rho_P))
plot(x=1:36, colMeans(rho_T))
summary(samples_list$phi_RP)
summary(samples_list$phi_T)
summary(samples_list$p_exp)
summary(samples_list$mu_p)
plot(samples_list$delta_p)
plot(samples_list$phi_T)


Ndf <- data.frame(
  year = 2018:2020,
  est_R=colMeans(samples_list$N_R_yr),
  hpd_R = HPDinterval(samples_list$N_R_yr),
  est_P=colMeans(samples_list$N_P_yr),
  hpd_P = HPDinterval(samples_list$N_P_yr),
  est_T=colMeans(samples_list$N_T_yr),
  hpd_T = HPDinterval(samples_list$N_T_yr)
)

ggplot(Ndf) + #geom_point(aes(x=year, y=est), size=3) +
  geom_ribbon(aes(x=year, ymin=hpd_R.lower, ymax=hpd_R.upper), fill="blue", alpha=0.2) +
  geom_path(aes(x=year, y=est_R), color="blue") +
  geom_ribbon(aes(x=year, ymin=hpd_P.lower, ymax=hpd_P.upper), fill="darkred", alpha=0.2) +
  geom_path(aes(x=year, y=est_P), color="darkred") +
  geom_ribbon(aes(x=year, ymin=hpd_R.lower, ymax=hpd_R.upper), fill="gold", alpha=0.2) +
  geom_path(aes(x=year, y=est_T), color="gold") +
  ylab("Abundance") + xlab("Year")




