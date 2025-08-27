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

  pi[1] <- 1-xi[1]
  pi[2] <- xi[1]
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

    phi[t] ~ dunif(0,1)
  }

  #' recruitment model
  for(t in 1:(K)){
    rho[t] ~ dunif(0,1)
    xi[t] <- rho[t]/(1-prod(1-rho[t:K]))
  }
  beta[1] <- xi[1]
  for(t in 2:K){
    beta[t] <- (rho[t] * prod(1-rho[1:(t-1)]))/((1-prod(1-rho[1:K])))
  }
  # rho[K] <- rho[K-1]
  # xi[K] <- 1

  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu
  nu_t[1:3,1] ~ dmulti(pi[1:3], nu)
  Nu[1] <- nu_t[2,1]
  for(t in 2:K){
    for(l in 1:3){
      Tu[l,1:3,t-1] ~ dmulti(Gamma[l,1:3,t-1], nu_t[l,t-1])
      nu_t[l,t] <- sum(Tu[1:3,l,t-1])
    }
    Nu[t] <-  nu_t[2,t]
  }

  for(i in 1:nobs){
    state[i,1:K] <- sample_state_Do(x[i,1:K], init = pi[1:3],
                                    probObs = pmat[1:3, 1:2, 1:K],
                                    probTrans = Gamma[1:3, 1:3, 1:(K-1)])
    alive[i,1:K] <- state[i,1:K]==2
  }

  for(t in 1:K){
    Nd[t] <- sum(alive[1:nobs, t])
    N[t] <- Nd[t] + Nu[t]
  }
})
