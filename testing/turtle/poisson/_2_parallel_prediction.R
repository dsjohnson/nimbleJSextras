library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(foreach)
library(future)
library(doFuture)
library(here)

### Read in MCMC sample from stage I

setwd(here("testing","turtle","poisson"))

load("honu_mcmc_sample.RData")
par_mcmc_list <- samples_list[c("xi","phi","theta","r","nu")]

# Define number of cores and chunk size
workers <- 10 #parallel::detectCores() - 2
total_iters <- dim(par_mcmc_list$xi)[1]
chunk_indices <- split(1:total_iters, cut(1:total_iters, workers))

# Raise the limit for the transfer
# options(future.globals.maxSize = 1500 * 1024^2)

# Create parallel sessions
plan(multisession, workers=workers) # Use multicore on Mac/Linux for faster shared memory
seeds <- 41 + 1:length(chunk_indices)

st <- Sys.time()

results <- foreach(j = 1:length(chunk_indices),
  .options.future = list(
    seed = TRUE,
    packages = c("nimble","nimbleEcology","nimbleJSextras"),
    globals=c("par_mcmc_list","honu_ch","chunk_indices","seeds")
    )
) %dofuture% {
  # Process a 'chunk' of iterations inside the worker
  # We subset the input matrices inside the worker to save memory
  set.seed(seeds[j])
  source("helper/make_hmm_mats.R")
  source("helper/predict_abundance.R")
  idx <- chunk_indices[[j]]
  worker_pred <- compileNimble(predict_abundance)
  pred <- worker_pred(
    xi = par_mcmc_list$xi[idx, , drop = FALSE],
    phi = par_mcmc_list$phi[idx, , drop = FALSE],
    theta = par_mcmc_list$theta[idx, , drop=FALSE],
    r = par_mcmc_list$r[idx, , , drop = FALSE],
    nd = 7,
    nu = par_mcmc_list$nu[idx],
    x = honu_ch
  )
  pred$N_avail
}

et <- Sys.time()
et-st
plan("sequential")

# Extract and stack N_avail
N_avail <- do.call(rbind, results) %>% mcmc()

###

# Yearly abundance of reproductive females
yr <- colnames(honu_ch) %>% str_remove("X") %>% as.numeric()
Ndf <- data.frame(year = yr, est=round(colMeans(N_avail)), hpd = HPDinterval(N_avail,0.9))
ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year") +
  theme_bw()
# ggsave("honu_Nt.pdf", width=6.5, height=4)
# get_trend <- function(x){
#   L <- length(x)
#   fit <- lm(log(x) ~ c(1:L))
#   return((exp(fit$coefficients[2])-1)*100)
# }
# ttt <- mcmc(apply(N_avail[,1:15], 1, get_trend))


# capture probs
r1 <- mcmc(par_mcmc_list$r[,,1])
p1 <- 1-exp(-r1)
r2 <- mcmc(par_mcmc_list$r[,,2])
p2 <- 1-exp(-r2)
pdf <- data.frame(year = yr,
                         est1=colMeans(p1), hpd1 = HPDinterval(p1, prob=0.9),
                         est2=colMeans(p2), hpd2 = HPDinterval(p2,prob=0.9))
ggplot(pdf) + geom_point(aes(x=year, y=est1), size=3, color="red") +
  geom_errorbar(aes(x=year, ymin=hpd1.lower, ymax=hpd1.upper), width=0, color="red") +
  # geom_path(aes(x=year, y=est1), color="red") +
  geom_point(aes(x=year, y=est2), size=3, color="blue") +
  geom_errorbar(aes(x=year, ymin=hpd2.lower, ymax=hpd2.upper), width=0, color="blue") +
  # geom_path(aes(x=year, y=est2), color="blue") +
  ylab("Capture probability") + xlab("Year") +
  theme_bw()

# Phi
phi <- mcmc(par_mcmc_list$phi)
phidf <- data.frame(year = yr[-length(yr)], est=colMeans(phi), hpd = HPDinterval(phi,0.9))
ggplot(phidf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0) +
  geom_path(aes(x=year, y=est)) + ylab("Survival probability") + xlab("Year") +
  theme_bw()

# theta
theta <- mcmc(par_mcmc_list$theta+2)
thetadf <- data.frame(year = yr[-length(yr)], est=colMeans(theta), hpd = HPDinterval(theta,0.9))
ggplot(thetadf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0) +
  geom_path(aes(x=year, y=est)) + ylab("Expected renesting interval") + xlab("Year") +
  theme_bw()


rho <- mcmc(samples_list$xi)
rhodf <- data.frame(year = yr, est=colMeans(rho), hpd = HPDinterval(rho,0.9))
ggplot(rhodf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0) +
  geom_path(aes(x=year, y=est)) + ylab("Entry probability") + xlab("Year") +
  theme_bw()


# Renest dist
theta <- samples_list$theta
rn_dist <- NULL
for(i in 1:length(theta)){
  rn_dist <- rbind(rn_dist, dpois(c(0:9)-1, theta[i]))
}
rn_dist <- mcmc(rn_dist)
rndf <- data.frame(year=c(1:10), est=colMeans(rn_dist), hpd = HPDinterval(rn_dist))
ggplot(rndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Renest probability") + xlab("Year") +
  theme_bw()

mean(theta+2)
HPDinterval(mcmc(theta+2))



