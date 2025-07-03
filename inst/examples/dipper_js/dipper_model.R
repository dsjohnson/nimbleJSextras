##
## CODE FOR JOLLY-SEBER CAPTURE-RECAPTURE MODEL
##

require(nimble)
require(nimbleEcology)
require(nimbleJSextras)

js_code <- nimbleCode({

  #' ---------------------------------------------------------------------------
  #' Unconditional detection probability
  #' ---------------------------------------------------------------------------
  pstar <- pstar_Do(
    init = pi[1:3],
    probObs = pmat[1:3, 1:2, 1:K],
    probTrans = Gamma[1:3, 1:3, 1:(K-1)],
    len=K,
    checkRowSums = 1
  )

  #' ---------------------------------------------------------------------------
  #' HMM likelihood for observed individuals
  #' ---------------------------------------------------------------------------
  for(i in 1:nobs){
    x[i, 1:K] ~ dJS_Do(
      init = pi[1:3],
      probObs = pmat[1:3, 1:2, 1:K],
      probTrans = Gamma[1:3, 1:3, 1:(K-1)],
      len = K, pstar=pstar, checkRowSums = 1
    )
  }
  n ~ dpois(lambda*pstar)

  #' ---------------------------------------------------------------------------
  #' Initial entry probability
  #' ---------------------------------------------------------------------------
  pi[1] <- 1-beta[1]
  pi[2] <- beta[1]
  pi[3] <- 0

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:K){
    pmat[1,1,t] <- 1
    pmat[1,2,t] <- 0
    pmat[2,1,t] <- 1-p[t]
    pmat[2,2,t] <- p[t]
    pmat[3,1,t] <- 1
    pmat[3,2,t] <- 0
  }

  for(t in 1:K){
    p[t] ~ dunif(0,1)
  }

  #' ---------------------------------------------------------------------------
  #' Use the code below for Royle & Dorazio (2008) parameterization
  #' ---------------------------------------------------------------------------
  # for(t in 2:(K-1)){
  #   p[t] ~ dunif(0,1)
  # }
  # p[1] <- 1
  # p[K] <- 1

  #' ---------------------------------------------------------------------------
  #' Transition probability matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:(K-1)){
    Gamma[1,1,t] <- 1-xi[t]
    Gamma[1,2,t] <- xi[t]
    Gamma[1,3,t] <- 0
    Gamma[2,1,t] <- 0
    Gamma[2,2,t] <- phi[t]
    Gamma[2,3,t] <- 1 - phi[t]
    Gamma[3,1,t] <- 0
    Gamma[3,2,t] <- 0
    Gamma[3,3,t] <- 1

    xi[t] <- beta[t+1]/(1-sum(beta[1:t]))
    phi[t] ~ dunif(0,1)
  }

  beta[1:K] ~ ddirch(mu_beta[1:K])
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu

  nu_t[1:3, 1:K] <- sample_state_undet_D(nu, pi[1:3], Gamma[1:3, 1:3, 1:(K-1)])

  # Bu[1:K] ~ dmulti(beta[1:K], nu)
  # Su[1] <- 0
  # Nu[1] <- Bu[1]
  # for(t in 2:K){
  #   Su[t] ~ dbinom(phi[t-1], Nu[t-1])
  #   Nu[t] <- Bu[t] + Su[t]
  # }

  for(i in 1:nobs){
    det_state[i,1:K] <- sample_state_det_Do(x[i,1:K], init = pi[1:3],
                                    probObs = pmat[1:3, 1:2, 1:K],
                                    probTrans = Gamma[1:3, 1:3, 1:(K-1)])
  }
  alive[1:nobs, 1:K] <- det_state[1:nobs,1:K]==2

  for(t in 1:K){
    Nd[t] <- sum(alive[1:nobs, t])
    N[t] <- Nd[t] + nu_t[2,t]
    # N[t] <- Nd[t] + Nu[t]
  }
})
