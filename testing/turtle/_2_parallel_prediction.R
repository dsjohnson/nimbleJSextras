library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(foreach)
library(doFuture)

### Read in MCMC sample from stage I

load("honu_mcmc_sample.RData")
par_mcmc_list <- samples_list[c("xi","phi","theta","p","nu")]

# Define number of cores and chunk size
n_cores <- 10 #parallel::detectCores() - 2
total_iters <- dim(par_mcmc_list$xi)[1]
chunk_indices <- split(1:total_iters, cut(1:total_iters, n_cores))

# Raise the limit for the transfer
# options(future.globals.maxSize = 1500 * 1024^2)

# Create parallel sessions
plan(multisession) # Use multicore on Mac/Linux for faster shared memory
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
    theta = par_mcmc_list$theta[idx],
    p = par_mcmc_list$p[idx, , drop = FALSE],
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
N_avail <- do.call(rbind, results)[,-c(1:10)] %>% mcmc()

###

# Yearly abundance of reproductive females
Ndf <- data.frame(year = 1998:2018, est=round(colMeans(N_avail)), hpd = HPDinterval(N_avail))
ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year") +
  theme_bw()
# ggsave("honu_Nt.pdf", width=6.5, height=4)
get_trend <- function(x){
  L <- length(x)
  fit <- lm(log(x) ~ c(1:L))
  return((exp(fit$coefficients[2])-1)*100)
}
ttt <- mcmc(apply(N_avail[,1:15], 1, get_trend))

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


plot(samples_list$Nsuper)
summary(samples_list$Nsuper)
plot(x=c(1994:2023), colMeans(samples_list$N))
plot(x=c(1994:2023), colMeans(samples_list$N_nest))

summary(samples_list$N)

xi <- samples_list$xi
xi <- sweep(xi,1,rowSums(xi),"/")
plot(colMeans(xi))

plot(samples_list$pstar)

summary(samples_list$phi)

summary(samples_list$p)
plot(colMeans(samples_list$p))
plot(samples_list$sig_tau)

plot(samples_list$theta+2)


