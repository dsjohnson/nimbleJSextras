##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  #' -----------------------------------------------------------------------------
  #' Capture probability
  #' -----------------------------------------------------------------------------
  for(t in 1:K){
    logit_p[t] ~ dnorm(mu_p, sd=sig_p)
    p[t] <- expit(logit_p[t])
  }
  mu_p ~ dnorm(0,sd=1.5)
  sig_p ~ dexp(1)

  #' ---------------------------------------------------------------------------
  #' Detection probabilities and matrix
  #' ---------------------------------------------------------------------------
  Pmats[1:K,1] <- 0
  Pmats[1:K,2] <- p[1:K]
  Pmats[1:K,3] <- 0

  #' ---------------------------------------------------------------------------
  #' Entry and exit probabilities
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    logit_phi[t] ~ dnorm(mu_phi, sd=sig_phi)
    phi[t] <- expit(logit_phi[t])
  }
  mu_phi ~ dnorm(0,sd=1.5)
  sig_phi ~ dexp(1)

  d[1] <- 1
  xi[1] <- 1
  for(t in 2:K){
    log_f[t-1] ~ dnorm(mu_f, sd=sig_f)
    f[t-1] <- exp(log_f[t-1])
    xi[t] <- d[t-1] * f[t-1]
    d[t] <- d[t-1] * (phi[t-1]+f[t-1])
  }
  mu_f ~ dnorm(0,sd=1.5)
  sig_f ~ dexp(1)

  for(t in 1:(K-1)){
    xi_tilde[t] <- xi[t]/sum(xi[t:K])
  }
  xi_tilde[K] <- 1

  #' ---------------------------------------------------------------------------
  #' Initial entry probability
  #' ---------------------------------------------------------------------------
  pi[1] <- 1-xi_tilde[1]
  pi[2] <- xi_tilde[1]
  pi[3] <- 0

  #' ---------------------------------------------------------------------------
  #' Transition probability matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    Gamma[1,1,t] <- 1-xi_tilde[t+1]
    Gamma[1,2,t] <- xi_tilde[t+1]
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

  pstar <- pstar_binom(
    init = pi[1:3],
    prob = Pmats[1:K,1:3],
    size = ones[1:K],
    probTrans = Gamma[1:3, 1:3, 1:(K-1)],
    len=K
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs){
    x[i, 1:K] ~ dJS_binom(
      init = pi[1:3],
      prob = Pmats[1:K,1:3],
      size = ones[1:K],
      probTrans = Gamma[1:3, 1:3, 1:(K-1)],
      pstar=pstar, len = K
    )
  }
  n ~ dpois(lambda*pstar)

  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n + nu

  nu_t[1:3, 1:K] <- sample_undet_binom(n=nu, init=pi[1:3], prob= Pmats[1:K,1:3],
                                       size=ones[1:K],
                                       probTrans=Gamma[1:3, 1:3, 1:(K-1)])

  for(i in 1:nobs){
    det_state[i,1:K] <- sample_det_binom(x=x[i,1:K], init=pi[1:3],
                                      prob=Pmats[1:K,1:3],
                                      size=ones[1:K],
                                      probTrans=Gamma[1:3, 1:3, 1:(K-1)])
  }
  available[1:nobs,1:K] <- det_state[1:nobs,1:K]==2

  for(t in 1:K){
    Nd[t] <- sum(available[1:nobs, t])
    N[t] <- Nd[t] + nu_t[2,t]
  }
})
