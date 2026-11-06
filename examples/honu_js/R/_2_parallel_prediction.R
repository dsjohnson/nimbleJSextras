library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(tidyverse)
library(coda)
library(foreach)
library(future)
library(doFuture)
library(here)
library(latex2exp)
library(cowplot)

setwd(here("honu_js"))

# Read in MCMC sample from stage I
load(file.path("processed_data_samples", "honu_mcmc_sample.RData"))
par_mcmc_list <- samples_list[c("xi","phi","tau","p","p0","Nu")]

# Define number of cores and chunk size
workers <- 10 #parallel::detectCores() - 2
total_iters <- dim(par_mcmc_list$xi)[1]
chunk_indices <- split(1:total_iters, cut(1:total_iters, workers))

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
  # Source utility functions
  source(file.path("R", "make_hmm_mats.R"))
  source(file.path("R", "predict_abundance.R"))
  
  # compile predictive sampler on the local worker
  worker_pred <- compileNimble(predict_abundance)
  
  idx <- chunk_indices[[j]]
  N_t <- worker_pred(
    xi = par_mcmc_list$xi[idx, , drop = FALSE],
    phi = par_mcmc_list$phi[idx],
    tau = par_mcmc_list$tau[idx, , drop=FALSE],
    p = par_mcmc_list$p[idx, ,drop = FALSE],
    p0 = par_mcmc_list$p0[idx, ,drop = FALSE],
    Nu = par_mcmc_list$Nu[idx],
    x = honu_ch
  )
  N_t
}

et <- Sys.time()
et-st
# Time difference of 7.657337 mins
plan("sequential")



# Extract and stack Yearly N_t sample
N_t <- do.call(rbind, results) %>% mcmc()

# Estimate super population abundance
summary(mcmc(samples_list$N))
# 
# Iterations = 1:50000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 50000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 9158.657        405.951          1.815         24.599 
# 
# 2. Quantiles for each variable:
#   
#   2.5%   25%   50%   75% 97.5% 
# 8478  8877  9115  9401 10056 
plot(mcmc(samples_list$N))
HPDinterval(mcmc(samples_list$N))
#       lower upper
# var1  8416  9964
quantile(samples_list$N, 0.99)
#      99% 
# 10339.01 

# Yearly abundance of reproductive females
yr <- colnames(honu_ch) %>% str_remove("X") %>% as.numeric()
obs <- colSums((honu_ch==2) + (honu_ch==3))
Ndf <- data.frame(
  year = yr,
  obs = obs,
  est=round(colMeans(N_t)),
  hpd = HPDinterval(mcmc(N_t,0.9))) %>%
  mutate(
    aaa=ifelse(year<1979, 0.25, 1)
  )
ggplot(Ndf) + geom_point(aes(x=year, y=est), alpha=Ndf$aaa, size=2) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2, alpha=Ndf$aaa) +
  geom_point(aes(x=year, y=obs), size=2, pch=3) +
  # geom_path(aes(x=year, y=est, color=col)) +
  ylab(TeX(r"(Yearly abundance ($N_t$))")) + xlab("Year") +
  theme_bw()
# ggsave(file.path("figures","abund.png"), dpi=600, width=6.5, height=4)
# ggsave(file.path("figures","abund.pdf"), width=6.5, height=4)

# Plot capture probabilities
p1 <- mcmc(par_mcmc_list$p)
p0 <- mcmc(par_mcmc_list$p0)
# p0 <- cbind(1,p)*eta
pdf <- data.frame(
  year = yr,
  est1=colMeans(p0), hpd1 = HPDinterval(p0, prob=0.9),
  est2=colMeans(p1), hpd2 = HPDinterval(p1,prob=0.9)
)
ggplot(pdf) + geom_point(aes(x=year, y=est1), size=3, color="red") +
  geom_errorbar(aes(x=year, ymin=hpd1.lower, ymax=hpd1.upper), width=0, color="red") +
  # geom_path(aes(x=year, y=est1), color="red") +
  geom_point(aes(x=year, y=est2), size=3, color="blue") +
  geom_errorbar(aes(x=year, ymin=hpd2.lower, ymax=hpd2.upper), width=0, color="blue") +
  # geom_path(aes(x=year, y=est2), color="blue") +
  ylab(TeX(r"(Capture probability ($p_t$))")) + xlab("Year") +
  theme_bw()
# ggsave(file.path("figures","capt_probs.png"), dpi=600, width=6.5, height=4)
# ggsave(file.path("figures","capt_probs.pdf"), width=6.5, height=4)


# Phi estimate
phi <- mcmc(par_mcmc_list$phi)
summary(phi)
# 
# Iterations = 1:50000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 50000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 9.128e-01      1.671e-03      7.471e-06      1.481e-05 
# 
# 2. Quantiles for each variable:
#   
#   2.5%    25%    50%    75%  97.5% 
# 0.9095 0.9117 0.9128 0.9140 0.9161 
HPDinterval(phi)
#           lower     upper
# var1 0.9095183 0.9160704
plot(phi)

# Estimate entry probability
xi <- mcmc(samples_list$xi)
xidf <- data.frame(year = yr, est=colMeans(xi), hpd = HPDinterval(xi,0.9))
ggplot(xidf) + geom_line(aes(x=year, y=est), linewidth=1) +
  geom_ribbon(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), alpha=0.25) +
  geom_hline(yintercept=1/ncol(honu_ch), color="red", lty=2) + 
  # geom_path(aes(x=year, y=est))
  ylab(TeX(r"(Entry probability ($\xi_t$))")) + xlab("Year") +
  theme_bw()
# ggsave(file.path("figures","entry_probs.png"), dpi=600, width=6.5, height=4)
# ggsave(file.path("figures","entry_probs.pdf"), width=6.5, height=4)

# Examine sigma_beta_xi
mean(samples_list$sigma_beta_xi)
# [1] 0.1830536
HPDinterval(mcmc(samples_list$sigma_beta_xi))
#              lower     upper
# var1 4.827946e-06 0.5601561
plot(mcmc(samples_list$sigma_beta_xi))


# estimate and plot remigration interval
ft <- mcmc(samples_list$ft_pdf)
ftdf <- data.frame(
  `Remig interval` = 2:11,
  est=colMeans(ft),
  hpd = HPDinterval(ft,0.9))
ggplot(ftdf) +
  geom_ribbon(aes(x=Remig.interval, ymin=hpd.lower, ymax=hpd.upper), fill="gray") +
  geom_path(aes(x=Remig.interval, y=est)) + ylab("PDF") + xlab("Remigration interval") +
  theme_bw()
# ggsave(file.path("figures","remig_int.png"), dpi=600, width=6.5, height=4)
# ggsave(file.path("figures","remig_int.pdf"), width=6.5, height=4)


# CDF of remigration interval
cdf_ft <- t(apply(ft, 1, cumsum)) |> mcmc()
HPDinterval(mcmc(cdf_ft[,6])) # nesters prior to study return on or before yr 7
#           lower upper
# var1 0.9955336     1

# Estimate average remigration interval
ft_e <- mcmc(apply(ft, 1, \(x) x%*%c(2:11)))
plot(ft_e)
summary(ft_e)
# 
# Iterations = 1:50000
# Thinning interval = 1 
# Number of chains = 1 
# Sample size per chain = 50000 
# 
# 1. Empirical mean and standard deviation for each variable,
# plus standard error of the mean:
#   
#   Mean             SD       Naive SE Time-series SE 
# 3.7447320      0.0350196      0.0001566      0.0009403 
# 
# 2. Quantiles for each variable:
#   
#   2.5%   25%   50%   75% 97.5% 
# 3.679 3.721 3.744 3.767 3.817 
HPDinterval(ft_e)
#          lower    upper
# var1 3.677621 3.814687
plot(ft_e)



# save(N_t, file=file.path("processed_data_samples","honu_abund_mcmc_sample.RData"))




