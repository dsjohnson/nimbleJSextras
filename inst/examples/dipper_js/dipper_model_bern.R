##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleJSextras)

js_code <- nimbleCode({

  for(t in 1:K){
    logit_p[t] ~ dnorm(mu_p, sd=sig_p)
    p[t] <- expit(logit_p[t])
    logit_rho[t] ~ dnorm(mu_rho, sd=sig_rho)
    rho[t] <- expit(logit_rho[t])
  }
  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mu_phi, sd=sig_phi)
    phi[t] <- expit(logit_phi[t])
  }
  sig_p ~ dexp(1)
  mu_p ~ dnorm(0,sd=1.5)
  sig_rho ~ dexp(1)
  mu_rho ~ dnorm(0,sd=10)
  mu_phi ~ dnorm(0,sd=1.5)
  sig_phi ~ dexp(1)

  xi[1] <- rho[1]
  for(t in 2:(K-1)){ xi[t] <- rho[t]/(1-prod(1-rho[t:K])) }
  xi[K] <- 1

  #' ---------------------------------------------------------------------------
  #' Initial entry probability
  #' ---------------------------------------------------------------------------
  pi[1] <- 1-xi[1]
  pi[2] <- xi[1]
  pi[3] <- 0

  #' ---------------------------------------------------------------------------
  #' Detection probabilities and matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:K){
    pmat[t,1] <- 0
    pmat[t,2] <- p[t]
    pmat[t,3] <- 0
  }

  #' ---------------------------------------------------------------------------
  #' Transition probability matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    Gamma[1,1,t] <- 1-xi[t+1]
    Gamma[1,2,t] <- xi[t+1]
    Gamma[1,3,t] <- 0
    Gamma[2,1,t] <- 0
    Gamma[2,2,t] <- phi[t]
    Gamma[2,3,t] <- 1 - phi[t]
    Gamma[3,1,t] <- 0
    Gamma[3,2,t] <- 0
    Gamma[3,3,t] <- 1
  }


  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- calc_pstar(
    piVector = pi,
    PArray = pmat,
    GammaArray = Gamma,
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------

  for(i in 1:nobs){
    x[i, 1:K] ~ dJS(
      pstar,
      pi,
      pmat,
      Gamma
    )
  }
  n ~ dpois(lambda*pstar)



  # pmat[1:K,1] <- 0
  # pmat[1:K,2] <- p
  # pmat[1:K,3] <- 0



  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu

  nu_t[1:3, 1:K] ~ dnu(nUndet=nu, pi=pi, PArray=pmat, GammaArray=Gamma, len=K)

  for(i in 1:nobs){
    z[i,1:K] ~ dState(pi, pmat, Gamma, x[i,1:K])
    }
  avail[1:nobs,1:K] <- z[1:nobs,1:K]==2

  for(t in 1:K){
    Nd[t] <- sum(avail[1:nobs, t])
    N[t] <- Nd[t] + nu_t[2,t]
  }
})
