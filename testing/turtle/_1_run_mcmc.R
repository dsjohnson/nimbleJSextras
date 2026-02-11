library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)

setwd("~/research/projects/r_packages/nimbleJSextras/testing/turtle")

honu_ch <- read.delim("honu_ch.txt",sep=" ")[,-1]
honu_ch <- honu_ch[,16:46]

# 2013 was not survyed
skip_surv <- which(colnames(honu_ch)=="X2013")
honu_ch[,skip_surv] <- 0

honu_ch <- honu_ch[rowSums(honu_ch)>0,]

# one animal return to breeding without a foraging period so it was removed
min_forage <- function(x){
  rrr <- rle(x)
  out <- max(rrr$lengths[rrr$values==1])
  return(out>1)
}
honu_ch <- honu_ch[as.logical(1-apply(honu_ch, 1, min_forage)),]

honu_ch <- honu_ch + 1


ch_collapsed <- collapse_ch(honu_ch)
K <- ncol(honu_ch)
n <- nrow(honu_ch)

source("helper/make_hmm_mats.R")
source("honu_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=K, nd=7, n_states=10, nobs_u=nrow(ch_collapsed$unique_ch)#, nobs=nrow(honu_ch)
  ),
  data = list(
    ux=ch_collapsed$unique_ch, w=ch_collapsed$w, skip_surv=skip_surv, #x=honu+ch,
    n=n
  ),
  inits = list(
    xi_raw = rep(1/K,K), theta=2,
    mu_phi=qlogis(0.9), eps=rep(0,K-1), sig_eps=2,
    mu_p=qlogis(0.9), tau=rep(0,K), sig_tau=2,
    lambda=2*n
  )
)

c_js_model <- compileNimble(js_model)
mcmcConf <- configureMCMC(js_model,
                          monitors=c(
                            "xi","theta",
                            "phi", "mu_phi", "eps", "sig_eps",
                            "p", "mu_p","tau", "sig_tau",
                            "nu","Nsuper","lambda", "pstar"#,"N_avail", "N_nest"
                          )
)
mcmcConf$removeSamplers('eps')
mcmcConf$addSampler(target = 'eps', type = 'RW_block')
mcmcConf$removeSamplers('tau')
mcmcConf$addSampler(target = 'tau', type = 'RW_block')
js_mcmc <- buildMCMC(mcmcConf)
c_js_mcmc <- compileNimble(js_mcmc)


set.seed(8675309)
st <- Sys.time()
samples <- runMCMC(c_js_mcmc, niter = 25000, nburnin = 5000, nchains = 1,
                   thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
et <- Sys.time()
et-st

samples_list <- as.list(c_js_mcmc$mvSamples)
save(samples_list, honu_ch, ch_collapsed, file="honu_mcmc_sample.RData")

