library(tidyverse)
library(nimble)
library(nimbleEcology)
library(nimbleJSextras)
library(coda)

setwd("~/research/projects/r_packages/nimbleJSextras/testing/turtle")

ch <- read.delim("yearly_detection.txt",sep=" ")[,-1]
ch <- ch[,16:46]
ch[,26] <- 0
x <- ch[rowSums(ch)>0,] + 1
# one animal return to breeding without a foraging period so it was removed
min_forage <- function(x){
  rrr <- rle(x)
  out <- max(rrr$lengths[rrr$values==2])
  return(out>1)
}
x <- x[as.logical(1-apply(x, 1, min_forage)),]

x_collapsed <- collapse_ch(x)
K <- ncol(x)
n <- nrow(x)

skip_surv <- which(colnames(x)=="X2013")

source("load_turtle_mats.R")
source("yearly_turtle_model.R")

js_model <- nimbleModel(
  code = js_code,
  constants = list(
    K=ncol(x), nd=7, n_states=10, nobs_u=nrow(x_collapsed$unique_ch)#, nobs=nrow(x)
  ),
  data = list(
    ux=x_collapsed$unique_ch, w=x_collapsed$w, skip_surv=skip_surv, #x=x,
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
                            "Nsuper","lambda", "pstar"#,"N", "N_nest","nu"
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

samples_list <- as.list(c_js_mcmc$mvSamples) %>% lapply(mcmc)

# Yearly abundance of reproductive females
Ndf <- data.frame(year = 1994:2023, est=colMeans(samples_list$N), hpd = HPDinterval(mcmc(samples_list$N)))
ggplot(Ndf) + geom_point(aes(x=year, y=est), size=3) +
  geom_errorbar(aes(x=year, ymin=hpd.lower, ymax=hpd.upper), width=0.2) +
  geom_path(aes(x=year, y=est)) + ylab("Abundance") + xlab("Year") +
  theme_bw()
ggsave("honu_Nt.pdf", width=6.5, height=4)


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
HPDinterval(theta+2)


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


