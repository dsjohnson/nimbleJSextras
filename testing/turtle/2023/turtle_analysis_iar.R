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
K <- ncol(honu)


#' -----------------------------------------------------------------------------
#' Load and compile model code
#' -----------------------------------------------------------------------------
source("turtle_model_pc.R")

iar_list <- make_nimble_icar(K, 4)

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K = K,
    nobs = nrow(honu),
    L = length(iar_list$adj),
    adj = iar_list$adj, wts=iar_list$weights, num=iar_list$num, c=iar_list$c
    ),
  data = list(
    x = honu, n = nrow(honu)
    ),
  inits = list(
    p_dot = 0.9,
    sig_xi = 1,
    mlogit_xi = rep(0,K),
    sig_delta = 1,
    mlogit_delta = rep(0,K),
    psi23_dot = 0.5,
    psi32_dot = 0.5,
    lambda = 2*nrow(honu)
    )
)
c_js_model <- compileNimble(js_model)
js_mcmc <- buildMCMC(js_model,
                     monitors=c("p_dot","delta","psi23_dot","psi32_dot", "xi",
                                "Nsuper","N_nest", "N_internest","N", "pstar"))
c_js_mcmc <- compileNimble(js_mcmc)


#' -----------------------------------------------------------------------------
#' Check initial convergence
#' -----------------------------------------------------------------------------
# set.seed(8675309)
# samples <- runMCMC(c_js_mcmc, niter = 20000, nburnin = 0, nchains = 3,
#                    thin = 1, samplesAsCodaMCMC = TRUE, progress = TRUE)
# c_js_mcmc$run(niter = 10000, reset = FALSE)
#
# gelman.diag(samples)


#' -----------------------------------------------------------------------------
#' Run full MCMC
#' -----------------------------------------------------------------------------

set.seed(8675309)
samples <- runMCMC(c_js_mcmc, niter = 70000, nburnin = 20000, nchains = 1,
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

summary(samples_list$delta)
HPDinterval(samples_list$delta)

summary(samples_list$p_dot)
HPDinterval(samples_list$p_dot)

summary(samples_list$beta)
HPDinterval(samples_list$beta)

xi <- sweep(samples_list$xi,1,rowSums(samples_list$xi),FUN="/")
xi_df <- data.frame(week = 1:ncol(honu), est=colMeans(xi), hpd = HPDinterval(xi))
ggplot(xi_df) + #geom_point(aes(x=week, y=est), size=1) +
  geom_ribbon(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), fill=gray(0.7)) +
  geom_path(aes(x=week, y=est)) + ylab("Prob. first nest") + xlab("Week") +
  theme_bw()

delta <- sweep(samples_list$delta,1,rowSums(samples_list$delta),FUN="/")
delta_df <- data.frame(week = 2:(ncol(honu)+1), est=colMeans(delta), hpd = HPDinterval(delta))
ggplot(delta_df) + #geom_point(aes(x=week, y=est), size=1) +
  geom_ribbon(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), fill=gray(0.7)) +
  geom_path(aes(x=week, y=est)) + ylab("Prob. repro. complete") + xlab("Week") +
  theme_bw()


Ndf <- data.frame(week = 1:ncol(honu), est=colMeans(samples_list$N), hpd = HPDinterval(mcmc(samples_list$N)))
ggplot(Ndf) +
  geom_ribbon(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), fill=gray(0.7)) +
  geom_path(aes(x=week, y=est)) +
  geom_point(aes(x=week, y=est), size=1) +
  ylab("Abundance") + xlab("Week") + theme_bw()
# ggsave("dipper_Nt.pdf", width=6.5, height=4)

Ndf <- data.frame(week = 1:ncol(honu), est=colMeans(samples_list$N_nest), hpd = HPDinterval(mcmc(samples_list$N_nest)))
ggplot(Ndf) +
  geom_ribbon(aes(x=week, ymin=hpd.lower, ymax=hpd.upper), fill=gray(0.7)) +
  geom_path(aes(x=week, y=est)) +
  geom_point(aes(x=week, y=est), size=1) +
  ylab("Nesting individuals") + xlab("Week") + theme_bw()
ggsave("dipper_Nt.pdf", width=6.5, height=4)





