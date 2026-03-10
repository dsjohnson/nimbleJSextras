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
dolphin_df <- read_csv("bottlenoseDolph.csv")
x <- as.matrix(dolphin_df) + 1
x_collapsed <- collapse_ch(x)
dates <- as.Date(colnames(dolphin_df))
year <- year(dates)-2017
h <- as.integer(diff(dates))/365
K <- ncol(x)
n <- nrow(x)

yidx <- matrix(NA,3,2)
for(j in 1:3){
  yidx[j,1] <- min(which(year==j))
  yidx[j,2] <- max(which(year==j))
}

#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("dolphin_js_model_caruso.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K,  nobs_u=nrow(x_collapsed$unique_ch), nobs=nrow(x)#, yidx=yidx
  ),
  data = list(
    ux=x_collapsed$unique_ch, w=x_collapsed$w, #x=x,
    n=n, h=h
  ),
  inits = list(
    rho = matrix(rep((1/K)/(1/K + 2-c(1:K)/K),3),K,3),
    phi_RP=0.71, phi_T=1.06e-8,
    alpha_raw=c(57/312, 81/312, 174/312),
    mu_p=qlogis(0.19), tau_raw=rep(0,K), #sig_p=0.01,
    delta_p=0.74,
    lambda=2*n
  )
)
c_js_model <- compileNimble(js_model)

mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            "rho","alpha",
                            "phi_RP","phi_T",
                            "p", "mu_p","tau","delta_p", #sig_p,
                            "Nsuper", "nu_raw","lambda",
                            "pstar","n2ll"#,"N_R_yr", "N_P_yr","N_T_yr"
                          )
)
mcmcConf$removeSamplers('tau_raw')
mcmcConf$addSampler(target = 'tau_raw', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)

samples <- runMCMC(c_js_mcmc, niter = 2000, nburnin = 0, nchains = 1, thin2=1, thin=1)

samples_list <- as.list(c_js_mcmc$mvSamples) %>% lapply(mcmc)

rho_R <- as.list(c_js_mcmc$mvSamples)$rho[,,1] %>% mcmc
rho_P <- as.list(c_js_mcmc$mvSamples)$rho[,,2] %>% mcmc
rho_T <- as.list(c_js_mcmc$mvSamples)$rho[,,3] %>% mcmc
p <- as.list(c_js_mcmc$mvSamples)$p[,,1] %>% mcmc


plot(samples_list$Nsuper)
summary(samples_list$alpha)
plot(x=dates, colMeans(rho_R))
plot(x=dates, colMeans(rho_P))
plot(x=dates, colMeans(rho_T))
plot(x=dates, colMeans(p))
plot(samples_list$phi_RP)
summary(p)
summary(samples_list$mu_p)
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




