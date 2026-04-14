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


### Read in MCMC sample from stage I

setwd(here("testing","turtle","final_analysis"))

# source("helper/make_hmm_mats.R")
# source("helper/predict_abundance.R")
load("honu_mcmc_sample.RData")
par_mcmc_list <- samples_list[c("xi","phi","tau","p","alpha","nu")]

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
  source("helper/make_hmm_mats2.R")
  source("helper/predict_abundance2.R")
  idx <- chunk_indices[[j]]
  worker_pred <- compileNimble(predict_abundance)
  pred <- worker_pred(
    xi = par_mcmc_list$xi[idx, , drop = FALSE],
    phi = par_mcmc_list$phi[idx],
    tau = par_mcmc_list$tau[idx, , drop=FALSE],
    p = par_mcmc_list$p[idx, ,drop = FALSE],
    alpha = par_mcmc_list$alpha[idx, ,drop = FALSE],
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
obs <- colSums((honu_ch==2) + (honu_ch==3))
Ndf <- data.frame(
  year = yr,
  obs = obs,
  est=round(colMeans(N_avail)),
  hpd = HPDinterval(N_avail,0.9)) %>%
  mutate(
    aaa=ifelse(year<1979, 0.25, 1)
  )
ggplot(Ndf) + geom_point(aes(x=year, y=est), alpha=Ndf$aaa, size=2) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2, alpha=Ndf$aaa) +
  geom_point(aes(x=year, y=obs), size=2, pch=3) +
  # geom_path(aes(x=year, y=est, color=col)) +
  ylab(TeX(r"(Yearly abundance ($N_t$))")) + xlab("Year") +
  theme_bw()
ggsave("abund.png", dpi="retina", width=6.5, height=4)

#
get_growth <- function(N){
  mean(diff(N)/N[length(N)]) * 100
}
gr <- mcmc(apply(N_avail[,-c(1:6)], 1, get_growth))
summary(gr)
HPDinterval(gr)

# capture probs
p <- mcmc(par_mcmc_list$p)[,-1]
alpha <- mcmc(par_mcmc_list$alpha)
p0 <- cbind(1,p)*alpha
pdf <- data.frame(
  year = yr,
  est1=colMeans(p0), hpd1 = HPDinterval(p0, prob=0.9),
  est2=c(NA,colMeans(p)), hpd2 = rbind(c(NA,NA),HPDinterval(p,prob=0.9))
  )
ggplot(pdf) + geom_point(aes(x=year, y=est1), size=3, color="red") +
  geom_errorbar(aes(x=year, ymin=hpd1.lower, ymax=hpd1.upper), width=0, color="red") +
  # geom_path(aes(x=year, y=est1), color="red") +
  geom_point(aes(x=year, y=est2), size=3, color="blue") +
  geom_errorbar(aes(x=year, ymin=hpd2.lower, ymax=hpd2.upper), width=0, color="blue") +
  # geom_path(aes(x=year, y=est2), color="blue") +
  ylab(TeX(r"(Capture probability ($p_t$))")) + xlab("Year") +
  theme_bw()
ggsave("capt_probs.png", dpi="retina", width=6.5, height=4)


# Phi
phi <- mcmc(par_mcmc_list$phi)
summary(phi)
HPDinterval(phi)
plot(phi)

# rho and xi
rho <- mcmc(samples_list$rho)
rhodf <- data.frame(year = yr, est=colMeans(rho), hpd = HPDinterval(rho,0.9))
ppp1 <- ggplot(rhodf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0) +
  # geom_path(aes(x=year, y=est))
  ylab(TeX(r"(Recruitment probability ($\rho_t$))")) + xlab("Year") +
  theme_bw()

xi <- mcmc(samples_list$xi)
xidf <- data.frame(year = yr, est=colMeans(xi), hpd = HPDinterval(xi,0.9))
ppp2 <- ggplot(xidf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0) +
  # geom_path(aes(x=year, y=est))
  ylab(TeX(r"(Entry probability ($\xi_t$))")) + xlab("Year") +
  theme_bw()

ppp <- plot_grid(ppp1, ppp2, ncol=1, labels = "AUTO")
ggsave(ppp, filename="rho_entry_probs.png", dpi="retina", width=6.5, height=6.5)


ft <- mcmc(samples_list$ft_pdf)
#xxx <- apply(ft, 1, \(x){min(which(cumsum(x)>0.99))})
ftdf <- data.frame(
  `Renest interval` = 2:11,
  est=colMeans(ft),
  hpd = HPDinterval(ft,0.9))
ggplot(ftdf) +
  #geom_point(aes(x=Renest.interval, y=est), size=2) +
  geom_ribbon(aes(x=Renest.interval, ymin=hpd.lower, ymax=hpd.upper), fill="gray") +
  geom_path(aes(x=Renest.interval, y=est)) + ylab("PDF") + xlab("Renesting interval") +
  theme_bw()
ggsave("renest_int.png", dpi="retina", width=6.5, height=4)



ft_e <- mcmc(apply(ft, 1, \(x) x%*%c(2:11)))
plot(ft_e)
summary(ft_e)
HPDinterval(ft_e)
hist(ft_e)

save(N_avail, file="honu_abund_mcmc_sample.RData")




