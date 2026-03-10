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

local_wd <- file.path(here(), "inst", "examples", "dipper_js","parallel")
setwd(local_wd)

# source("helper/make_hmm_mats.R")
# source("helper/predict_abundance.R")
load("dipper_stage1_sample.RData")
par_mcmc_list <- samples_list[c("xi","phi","p","nu")]

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
    globals=c("par_mcmc_list","x","chunk_indices","seeds")
    )
) %dofuture% {
  # Process a 'chunk' of iterations inside the worker
  # We subset the input matrices inside the worker to save memory
  set.seed(seeds[j])
  source("helper/predict_abundance.R")
  idx <- chunk_indices[[j]]
  worker_pred <- compileNimble(predict_abundance)
  N_avail <- worker_pred(
    xi = par_mcmc_list$xi[idx, , drop = FALSE],
    phi = par_mcmc_list$phi[idx, , drop = FALSE],
    p = par_mcmc_list$p[idx, , drop = FALSE],
    nu = par_mcmc_list$nu[idx],
    x = x
  )
  N_avail
}

et <- Sys.time()
et-st
plan("sequential")

# Extract and stack N_avail
N_avail <- do.call(rbind, results) %>% mcmc()

###

# Yearly abundance of reproductive females
Ndf <- data.frame(year = 1:ncol(x), est=colMeans(N_avail), hpd = HPDinterval(mcmc(N_avail)))
ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year") +
  theme_bw()

summary(mcmc(N_avail))
HPDinterval(mcmc(N_avail))


