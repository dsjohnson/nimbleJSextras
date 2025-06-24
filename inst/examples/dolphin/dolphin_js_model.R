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
  # Resident
  pi[1] <- 1-beta_R[1]
  pi[2] <- beta_R[1]
  pi[3] <- 0
  pi[4] <- 0
  # Part-time
  pi[5] <- 1-beta_P[1]
  pi[6] <- beta_R[1]
  pi[7] <- 0
  pi[8] <- 0
  # Transient
  pi[9] <- 1-beta_T[1]
  pi[10] <- beta_T[1]
  pi[11] <- 0
  pi[12] <- 0

  #' ---------------------------------------------------------------------------
  #' Detection matrix
  #' ---------------------------------------------------------------------------
  for(t in 1:K){
    # Resident
    pmat[1,1,t] <- 1
    pmat[1,2,t] <- 0
    pmat[2,1,t] <- 1 - p_[t]
    pmat[2,2,t] <- p_[t]
    pmat[3,1,t] <- 1
    pmat[3,2,t] <- 0
    pmat[4,1,t] <- 1
    pmat[4,2,t] <- 0
    # Part-time resident
    pmat[5,1,t] <- 1
    pmat[5,2,t] <- 0
    pmat[6,1,t] <- 1 - p[t]
    pmat[6,2,t] <- p[t]
    pmat[7,1,t] <- 1
    pmat[7,2,t] <- 0
    pmat[8,1,t] <- 1
    pmat[8,2,t] <- 0
    # Transient
    pmat[9,1,t] <- 1
    pmat[9,2,t] <- 0
    pmat[10,1,t] <- 1 - p_[t]
    pmat[10,2,t] <- p[t]
    pmat[11,1,t] <- 1
    pmat[11,2,t] <- 0
    pmat[12,1,t] <- 1
    pmat[12,2,t] <- 0
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

  #' Constant per capita recruitment model
  f ~ dunif(0,2)
  d[1] <- 1
  prop_beta[1] <- 1
  for(t in 2:K){
    d[t] <- d[t-1] * (phi[t-1]+f)
    prop_beta[t] <- f * d[t]
  }
  norm <- sum(prop_beta[1:K])
  beta[1:K] <- prop_beta[1:K]/norm

  #' lambda prior
  lambda ~ dgamma(1.0e-6, 1.0e-6)

  #' ---------------------------------------------------------------------------
  #' Code below is for the posterior predicted abundance
  #' ---------------------------------------------------------------------------

  nu ~ dpois(lambda*(1-pstar))
  Nsuper <- n+nu
  Bu[1:K] ~ dmulti(beta[1:K], nu)

  Su[1] <- 0
  Nu[1] <- Bu[1]
  for(t in 2:K){
    Su[t] ~ dbinom(phi[t-1], Nu[t-1])
    Nu[t] <- Bu[t] + Su[t]
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

  Bp[1:K] ~ dmulti(beta[1:K], Nsuper)
  Sp[1] <- 0
  Np[1] <- Bp[1]
  for(t in 2:K){
    Sp[t] ~ dbinom(phi[t-1], Np[t-1])
    Np[t] <- Bp[t] + Sp[t]
  }



})
